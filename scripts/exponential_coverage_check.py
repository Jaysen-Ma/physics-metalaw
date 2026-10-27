#!/usr/bin/env python3
"""Coverage checks for untrimmed exponential fitting.

Compare seven equal-rate branches under homogeneous equations, heterogeneous
sizes/supports with within-equation dependence, and chapter dependence with
chapter resampling. These synthetic count vectors are not expression trees.
Results check specified simulations, not coverage on the actual corpora.
"""
from __future__ import annotations
import argparse
from concurrent.futures import ProcessPoolExecutor
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
from scipy.stats import beta as beta_dist

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from physics_metalaw.exponential_fit import bootstrap_type_matrix, simultaneous_intervals, summarise_bootstrap

SEED = 20261026
BETA = 0.35
SCENARIOS = ('homogeneous_equations', 'heterogeneous_equations', 'chapter_dependence')


def exponential_probabilities(n_types, beta=BETA):
    p = np.exp(-beta * np.arange(1, n_types + 1))
    return p / p.sum()


def generate_group(rng, scenario, group):
    sizes = [89, 221, 312, 126, 105, 73, 42]
    supports = [9, 10, 11, 12, 13, 14, 15]
    p = exponential_probabilities(12 if scenario == 'homogeneous_equations' else supports[group])
    if scenario == 'homogeneous_equations':
        return rng.multinomial(12, p, size=80).astype(float), None
    if scenario == 'heterogeneous_equations':
        # Marginal probabilities retain beta=0.35; tokens within an equation
        # share a random probability vector, producing extra dependence.
        rows = []
        for _ in range(sizes[group]):
            q = rng.dirichlet((10 + 5 * group) * p)
            rows.append(rng.multinomial(int(rng.integers(8, 25)), q))
        return np.asarray(rows, dtype=float), None
    n_chapters = [9, 23, 28, 16, 16, 13, 4][group]
    rows, clusters = [], []
    for _ in range(n_chapters):
        q = rng.dirichlet(30 * p)
        n = int(rng.integers(3, 10))
        start = len(rows)
        rows.extend(rng.multinomial(12, q, size=n))
        clusters.append(np.arange(start, start + n))
    return np.asarray(rows, dtype=float), clusters


def one_dataset(task):
    scenario_index, dataset, n_boot = task
    scenario = SCENARIOS[scenario_index]
    rng = np.random.default_rng(np.random.SeedSequence(SEED, spawn_key=(scenario_index, dataset)))
    hats, draws, covered, bias = [], [], [], []
    failed = boundary = 0
    for group in range(7):
        matrix, clusters = generate_group(rng, scenario, group)
        fit, sample, status = bootstrap_type_matrix(matrix, threshold=None, n_boot=n_boot, rng=rng, clusters=clusters)
        summary = summarise_bootstrap(fit, sample, status)
        valid = summary['interval_validity'] == 'approximate_percentile'
        lo, hi = summary['ci95']
        covered.append(bool(lo <= BETA <= hi) if valid and lo is not None else None)
        bias.append(None if fit.beta is None else fit.beta - BETA)
        failed += summary['n_failed'] + summary['n_degenerate']
        boundary += summary['n_boundary']
        hats.append(np.nan if fit.beta is None else fit.beta)
        draws.append(sample)
    pairs = [(a,b) for a in range(7) for b in range(a+1,7)]
    delta = np.array([hats[a]-hats[b] for a,b in pairs])
    samples = np.column_stack([draws[a]-draws[b] for a,b in pairs])
    result = simultaneous_intervals(delta, samples)
    valid = result['interval_validity'] == 'approximate_simultaneous'
    false_positive = any(x['excludes_zero'] for x in result['intervals']) if valid else None
    return {'scenario': scenario, 'dataset': dataset, 'coverage_by_branch': covered,
            'estimation_error_by_branch': bias, 'any_false_positive': false_positive,
            'interval_validity': result['interval_validity'], 'failed_draws': failed,
            'boundary_draws': boundary}


def binomial_summary(flags):
    values = [bool(x) for x in flags if x is not None]
    n, k = len(values), sum(values)
    lo = 0.0 if k == 0 else float(beta_dist.ppf(.025, k, n-k+1))
    hi = 1.0 if k == n else float(beta_dist.ppf(.975, k+1, n-k))
    return {'n': n, 'count': k, 'n_missing': len(flags)-n,
            'rate': k/n if n else None, 'monte_carlo_ci95': [lo,hi] if n else [None,None]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--datasets', type=int, default=200)
    parser.add_argument('--boot', type=int, default=400)
    parser.add_argument('--workers', type=int, default=3)
    args = parser.parse_args()
    if args.datasets < 1 or args.boot < 2 or args.workers < 1:
        parser.error('datasets/workers must be positive; boot must be at least 2')
    tasks = [(s,d,args.boot) for s in range(len(SCENARIOS)) for d in range(args.datasets)]
    records = []
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        for i, result in enumerate(pool.map(one_dataset, tasks), 1):
            records.append(result)
            if i % 20 == 0:
                print(f'{i}/{len(tasks)} simulation datasets completed', flush=True)
    scenarios = {}
    for name in SCENARIOS:
        subset = [x for x in records if x['scenario'] == name]
        scenarios[name] = {
            'simultaneous_false_positive_rate': binomial_summary([x['any_false_positive'] for x in subset]),
            'coverage_by_branch': [binomial_summary([x['coverage_by_branch'][g] for x in subset]) for g in range(7)],
            'mean_estimation_error_by_branch': [float(np.mean([x['estimation_error_by_branch'][g] for x in subset if x['estimation_error_by_branch'][g] is not None])) for g in range(7)],
            'failed_draws': sum(x['failed_draws'] for x in subset),
            'boundary_draws': sum(x['boundary_draws'] for x in subset)}
    output = ROOT/'results/exponential_fit'
    if (args.datasets,args.boot) != (200,400):
        output = output/'diagnostics'/f'coverage_{args.datasets}_{args.boot}'
    output.mkdir(parents=True,exist_ok=True)
    payload = {'seed': SEED, 'n_datasets_per_scenario':args.datasets, 'n_boot':args.boot,
               'beta':BETA, 'threshold':None, 'n_branches':7, 'n_comparisons':21,
               'resampling':{'homogeneous_equations':'equations','heterogeneous_equations':'equations','chapter_dependence':'chapters'},
               'limitations':'Synthetic count vectors, not expression trees. Finite simulation precision. Chapter scenario includes a branch with only four clusters. This does not establish coverage on real corpora.',
               'code_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__).resolve(),ROOT/'physics_metalaw/exponential_fit.py']},
               'scenarios':scenarios}
    (output/'coverage_check.json').write_text(json.dumps(payload,indent=2)+'\n')
    (output/'coverage_datasets.json').write_text(json.dumps(records,indent=2)+'\n')
    print(json.dumps(scenarios,indent=2))

if __name__ == '__main__':
    main()
