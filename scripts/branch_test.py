#!/usr/bin/env python3
"""Is β the same in every branch of physics? Single-source test on the Feynman corpus.

Corpus: the 100 AI Feynman benchmark expressions (data/catalog/feynman.json) plus the
Feynman Lectures Vol. I-III catalogs (data/catalog/feynman_lectures_v{1,2,3}.json),
deduplicated against the benchmark and across volumes.

Model: p_r = exp(-β r) / Σ_{s=1..R} exp(-β s) over the R ranks observed in the group
(operator types ranked by their total count). β minimises the sum of squared
log-frequency residuals, one residual per rank.
Errors: bootstrap over expressions (2000 replicates, re-ranked each time) and the
least-squares residual standard error.
Tests: profile likelihood ratio of those residuals (shared β vs one β per group),
χ² and permutation p-values (expressions shuffled between groups), ΔAIC/ΔBIC,
Bonferroni for families of tests.

Run from the repository root (about 80 s on one CPU core):

    .venv/bin/python scripts/branch_test.py            # writes results/branch_test.json
    .venv/bin/python scripts/branch_test.py --quick    # 200 permutations, for a smoke test

Outputs: results/branch_test.json (every number in the README and writeup),
results/branch_beta.csv, results/branch_tests.csv, writeup/generated/numbers.tex.
"""
from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import re
import sys
import time
from collections import Counter
from pathlib import Path

import numpy as np
import sympy
from scipy.stats import chi2, norm

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from physics_metalaw import metalaw as M  # noqa: E402
from physics_metalaw.catalog import CATALOG_DIR, record_as_entry, scatter_records  # noqa: E402

SEED = 20261013
BOOT = 2000
RESULTS = ROOT / "results"
GENERATED = ROOT / "writeup" / "generated"

_spec = importlib.util.spec_from_file_location("bflc", ROOT / "scripts" / "build_feynman_lectures_catalog.py")
BF = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(BF)

# Fine branch labels (as written in the catalogs) -> analysis groups.
FINE = {"Classical mechanics": "Classical", "Special relativity": "SR", "Electromagnetism": "EM",
        "Waves and optics": "Waves+optics", "Thermodynamics and kinetic theory": "Thermo",
        "Quantum mechanics": "QM", "Fluid mechanics": "Fluid+continuum", "Continuum mechanics": "Fluid+continuum"}
ORDER = ["SR", "QM", "EM", "Waves+optics", "Classical", "Thermo", "Fluid+continuum"]
LONG = {"SR": "Special relativity", "QM": "Quantum mechanics", "EM": "Electromagnetism",
        "Waves+optics": "Waves and optics", "Classical": "Classical mechanics",
        "Thermo": "Thermodynamics and kinetic theory", "Fluid+continuum": "Fluid and continuum mechanics"}
COARSE = {"Classical": "Classical+SR", "SR": "Classical+SR", "EM": "EM+waves", "Waves+optics": "EM+waves",
          "QM": "QM", "Thermo": "Thermo", "Fluid+continuum": "Fluid+continuum"}
# The quoted uncertainty is the bootstrap. The least-squares residual standard
# error treats ranks as independent and is smaller than the bootstrap for the
# smaller branches, so it is reported but not quoted.
QUOTE_FIT = set()

# Sensitivity check: the alternative branch for each chapter flagged ambiguous in the
# chapter-to-branch table (None = no single alternative; the chapter is dropped).
# Relativistic electrodynamics (I.34, II.25, II.26, II.28) is SR in the main mapping
# (BF.BRANCH_OVERRIDE); its alternative here is electromagnetism.
ALT = {
    ("I", 5): None, ("I", 7): None, ("I", 12): None, ("I", 23): "Electromagnetism",
    ("I", 24): "Electromagnetism", ("I", 25): "Electromagnetism", ("I", 28): "Waves and optics",
    ("I", 31): "Electromagnetism", ("I", 32): "Waves and optics", ("I", 34): "Electromagnetism",
    ("I", 35): None, ("I", 36): None, ("I", 37): "Waves and optics", ("I", 41): "Quantum mechanics",
    ("I", 50): None, ("I", 51): "Fluid mechanics",
    ("II", 1): None, ("II", 10): None, ("II", 11): "Thermodynamics and kinetic theory", ("II", 12): None,
    ("II", 19): None, ("II", 22): None, ("II", 25): "Electromagnetism", ("II", 26): "Electromagnetism",
    ("II", 28): "Electromagnetism", ("II", 32): "Electromagnetism", ("II", 33): "Electromagnetism",
    ("II", 34): "Quantum mechanics", ("II", 35): "Quantum mechanics", ("II", 36): None, ("II", 37): None,
    ("III", 4): "Thermodynamics and kinetic theory", ("III", 13): None, ("III", 14): None,
    ("III", 15): None, ("III", 20): None, ("III", 21): None,
}


# ---------------------------------------------------------------- corpus
def _vc(url):
    m = re.search(r"/(I{1,3})_(\d+)\.html", url or "")
    return (m.group(1), int(m.group(2))) if m else None


def feynman_corpus():
    """[(record, branch label, meta)] for the benchmark plus deduplicated lecture records."""
    bench = scatter_records("feynman")
    labels = {r["id"]: r.get("branch") for r in json.loads((CATALOG_DIR / "feynman.json").read_text())["records"]}
    def bench_label(r):  # same main-mapping rule for benchmark expressions from relativistic-EM chapters
        b = labels.get(r.id)
        return "Special relativity" if b == "Electromagnetism" and _vc(r.source_url) in BF.BRANCH_OVERRIDE else b
    rows = [(r, bench_label(r), dict(source="benchmark", volume=None, chapter=None, ambiguous=False,
                                     notation=False)) for r in bench]
    bench_shape = {(_vc(r.source_url), BF.signature(r.formula), BF.consts(r.formula)) for r in bench}
    seen = set()
    for r in bench:
        try:
            seen.add(sympy.srepr(sympy.sympify(r.formula)))
        except Exception:  # noqa: BLE001
            pass
    dedupe, added = Counter(), Counter()
    for v in (1, 2, 3):
        key = f"feynman_lectures_v{v}"
        raw = {d["id"]: d for d in json.loads((CATALOG_DIR / f"{key}.json").read_text())["records"]}
        for r in scatter_records(key):
            d = raw[r.id]
            if d.get("benchmark_overlap"):
                dedupe[f"v{v}: one of the 100 benchmark expressions"] += 1
                continue
            if (_vc(r.source_url), BF.signature(r.formula), BF.consts(r.formula)) in bench_shape:
                dedupe[f"v{v}: same chapter, shape and constants as a benchmark expression"] += 1
                continue
            try:
                ex = sympy.srepr(sympy.sympify(r.formula))
            except Exception:  # noqa: BLE001
                ex = r.formula
            if ex in seen:
                dedupe[f"v{v}: exact duplicate of an earlier expression"] += 1
                continue
            seen.add(ex)
            added[f"v{v}"] += 1
            names = set(d.get("variable_names") or []) | set(re.findall(r"[A-Za-z_]\w*", r.formula))
            notation = any(n.startswith("avg") for n in names) or "i" in names or "I" in names
            rows.append((r, d.get("branch"), dict(source=f"lectures_v{v}", volume=d["volume"], chapter=d["chapter"],
                                                  ambiguous=bool(d.get("ambiguous_chapter")), notation=notation)))
    return rows, dict(dedupe), dict(added)


# ---------------------------------------------------------------- log least squares
class Counts:
    """Per-expression operator counts on one common set of operator types."""

    def __init__(self, records):
        per = [M.rank_counts([record_as_entry(r)]) for r in records]
        keys = sorted({k for c in per for k in c})
        self.mat = np.array([[c.get(k, 0) for k in keys] for c in per], float)

    def obs(self, idx):
        t = self.mat[idx].sum(axis=0)
        return np.sort(t[t > 0])[::-1]


def compare(C, groups):
    """Log-least-squares comparison. Λ is the first element."""
    return M.logls_compare([C.obs(g) for g in groups])


def test(C, name, groups, labels, reps, rng, bonf=1, family=None):
    G = len(groups)
    stat, b0, bs, _rss0, _rss1, n_ranks, N = compare(C, groups)
    pool = np.concatenate(groups)
    sizes = [len(g) for g in groups]
    cnt = 0
    for _ in range(reps):
        perm = rng.permutation(pool)
        parts, s = [], 0
        for n in sizes:
            parts.append(perm[s:s + n])
            s += n
        cnt += compare(C, parts)[0] >= stat - 1e-9
    p_perm = (cnt + 1) / (reps + 1)
    p_chi2 = float(chi2.sf(stat, G - 1))
    # Gaussian profile likelihood on the log-frequency residuals. n is the
    # number of ranks, the observations in the least-squares fit.
    daic = -stat + 2 * (G - 1)
    dbic = -stat + (G - 1) * np.log(n_ranks)
    row = dict(name=name, family=family, labels=list(labels), n=sizes, beta_sep=bs, beta_shared=b0,
               N=int(N), n_ranks=n_ranks,
               LR=stat, df=G - 1, p_chi2=p_chi2, p_perm=p_perm, perm_reps=reps, p_perm_floor=1 / (reps + 1),
               bonf=bonf, p_perm_bonf=min(1.0, p_perm * bonf), p_chi2_bonf=min(1.0, p_chi2 * bonf),
               dAIC=daic, dBIC=dbic)
    print(f"  {name}: n={dict(zip(labels, sizes))} β=" + ", ".join(f"{b:.3f}" for b in bs)
          + f" shared={b0:.3f} LR={stat:.2f} df={G - 1} p_chi2={p_chi2:.3g} p_perm={p_perm:.4f}"
          + (f" (Bonf ×{bonf}: {row['p_perm_bonf']:.4f})" if bonf > 1 else "")
          + f" ΔAIC={daic:+.1f} ΔBIC={dbic:+.1f}", flush=True)
    return row


def fit(records, seed=SEED):
    d = M.beta_with_errors([record_as_entry(r) for r in records], reps=BOOT, seed=seed)
    return {k: (float(v) if isinstance(v, (float, np.floating)) else v) for k, v in d.items()}


# ---------------------------------------------------------------- analysis of one labelling
def groups_of(rows, relabel=None, drop=None):
    """Fine groups -> list of row indices (into rows)."""
    out = {}
    for i, (r, b, meta) in enumerate(rows):
        if drop and drop(meta, b):
            continue
        if relabel:
            b = relabel(meta, b)
        g = FINE.get(b)
        if g:
            out.setdefault(g, []).append(i)
    return {k: np.array(out[k]) for k in ORDER if k in out and len(out[k]) >= 5}


def analyse(rows, C, groups, reps, rng, label, full=True):
    print(f"\n== {label} ==", flush=True)
    recs = [r for r, _, _ in rows]
    table = []
    for k in ORDER:
        if k not in groups:
            continue
        d = fit([recs[i] for i in groups[k]])
        d.update(branch=k, quoted_sigma=d["fit_sigma"] if k in QUOTE_FIT else d["boot_sigma"],
                 quoted_sigma_kind="least squares" if k in QUOTE_FIT else "bootstrap")
        table.append(d)
        print(f"  {k:16s} n={d['n']:4d} N_ops={d['n_instances']:5d} R={d['n_ranks']:2d} β={d['beta']:.3f} "
              f"boot={d['boot_sigma']:.3f} [{d['boot_p16']:.3f},{d['boot_p84']:.3f}] fit={d['fit_sigma']:.3f}",
              flush=True)
    keys = [k for k in ORDER if k in groups]
    tests = [test(C, "all branches", [groups[k] for k in keys], keys, reps, rng, family="global")]
    for k in keys:
        rest = np.concatenate([groups[u] for u in keys if u != k])
        tests.append(test(C, f"{k} vs rest", [groups[k], rest], [k, "rest"], reps, rng, bonf=len(keys),
                          family="each_vs_rest"))
    if full:
        npairs = len(keys) * (len(keys) - 1) // 2
        for i, a in enumerate(keys):
            for b in keys[i + 1:]:
                tests.append(test(C, f"{a} vs {b}", [groups[a], groups[b]], [a, b], reps, rng, bonf=npairs,
                                  family="pairwise"))
    return table, tests


def power(rows, C, groups, table, rng):
    print("\n== power ==", flush=True)
    sig = {d["branch"]: d["quoted_sigma"] for d in table}
    keys = [d["branch"] for d in table]
    z = norm.isf(0.025) + norm.isf(0.2)
    mdd = [dict(a=a, b=b, mdd80=float(z * np.hypot(sig[a], sig[b]))) for i, a in enumerate(keys) for b in keys[i + 1:]]
    sims = []
    for k in keys:
        dk = next(d for d in table if d["branch"] == k)
        rest = np.concatenate([groups[u] for u in keys if u != k])
        orest = C.obs(rest)
        brest = M.fit_beta_logls(orest)[0]
        infl = max(1.0, dk["boot_sigma"] / dk["fit_sigma"])
        for delta in (0.02, 0.04, 0.06, 0.08):
            hits = 0
            for _ in range(1000):
                oa = []
                for n, R, b in ((dk["n_instances"], dk["n_ranks"], brest + delta), (orest.sum(), orest.size, brest)):
                    rr = np.arange(1, R + 1)
                    p = np.exp(-b * rr)
                    p /= p.sum()
                    o = np.sort(rng.multinomial(int(round(n / infl ** 2)), p))[::-1].astype(float)
                    o = o[o > 0]
                    if o.size >= 2:
                        oa.append(o)
                hits += len(oa) == 2 and chi2.sf(M.logls_compare(oa)[0], 1) < 0.05
            sims.append(dict(group=k, delta=delta, power=hits / 1000, overdispersion=infl))
            print(f"  {k:16s} vs rest Δβ=+{delta:.2f} power={hits / 1000:.2f}", flush=True)
    return mdd, sims


def rank_frequency(C, idx, beta):
    o = C.obs(idx)
    r = np.arange(1, o.size + 1)
    lp = -beta * r
    return dict(rank=r.tolist(), freq=(o / o.sum()).tolist(), model=np.exp(lp - np.logaddexp.reduce(lp)).tolist())


def corpora(rows):
    print("\n== corpora ==", flush=True)
    out = []
    for key, label in (("feynman", "Feynman (AI Feynman benchmark)"), ("wiki", "Wikipedia"),
                       ("inflation", "Encyclopaedia Inflationaris"), ("leemis", "Leemis"), ("fr", "f(R)")):
        d = fit(scatter_records(key))
        d.update(key=key, label=label)
        out.append(d)
    d = fit([r for r, _, _ in rows])
    d.update(key="feynman_all", label="Feynman (benchmark + Lectures I-III)")
    out.insert(1, d)
    for d in out:
        print(f"  {d['label']:40s} n={d['n']:4d} β={d['beta']:.3f} boot={d['boot_sigma']:.3f} fit={d['fit_sigma']:.3f}")
    return out


def corpus_test(reps, rng):
    """EI vs the Feynman benchmark (the two released corpora of the upstream paper)."""
    a, b = scatter_records("feynman"), scatter_records("inflation")
    C = Counts(a + b)
    return test(C, "EI vs Feynman benchmark", [np.arange(len(a)), np.arange(len(a), len(a) + len(b))],
                ["Feynman benchmark", "EI"], reps, rng, family="corpus")


# ---------------------------------------------------------------- outputs
def fmt(x, nd=3):
    return f"{x:.{nd}f}"


def write_numbers_tex(res):
    """LaTeX macros and table rows used by writeup/main.tex."""
    GENERATED.mkdir(parents=True, exist_ok=True)
    main = res["main"]
    lines = ["% Generated by scripts/branch_test.py. Do not edit.", ]
    m = lambda name, val: lines.append(f"\\newcommand{{\\{name}}}{{{val}}}")
    alpha = {"SR": "SR", "QM": "QM", "EM": "EM", "Waves+optics": "Waves", "Classical": "Classical",
             "Thermo": "Thermo", "Fluid+continuum": "Fluid"}
    for d in main["table"]:
        a = alpha[d["branch"]]
        m(f"beta{a}", fmt(d["beta"]))
        m(f"sig{a}", fmt(d["quoted_sigma"]))
        m(f"n{a}", d["n"])
    g = main["tests"][0]
    m("LRall", f"{g['LR']:.1f}")
    m("dfall", g["df"])
    m("dAICall", f"{g['dAIC']:.0f}")
    m("dBICall", f"{g['dBIC']:.0f}")
    m("pchiall", tex_sci(g["p_chi2"]))
    m("permreps", main["tests"][0]["perm_reps"])
    m("permfloor", f"{1 / (main['tests'][0]['perm_reps'] + 1):.5f}")
    for c in res["corpora"]:
        a = {"feynman": "FeynBench", "feynman_all": "FeynAll", "wiki": "Wiki", "inflation": "EI", "leemis": "Leemis",
             "fr": "fR"}[c["key"]]
        m(f"beta{a}", fmt(c["beta"]))
        m(f"sig{a}", fmt(c["boot_sigma"]))
        m(f"n{a}", c["n"])
    pooled = res["pooled_all"]
    m("betaPooled", fmt(pooled["beta"]))
    m("sigPooled", fmt(pooled["boot_sigma"]))
    m("lsigPooled", fmt(pooled["fit_sigma"]))
    m("nPooled", pooled["n"])
    m("nLabelled", res["pooled_labelled"]["n"])
    for k, v in res["counts"].items():
        m(k, v)
    (GENERATED / "numbers.tex").write_text("\n".join(lines) + "\n")

    rows = []
    for d in main["table"]:
        evr = next(t for t in main["tests"] if t["name"] == f"{d['branch']} vs rest")
        rows.append(f"{LONG[d['branch']]} & {d['n']} & {d['n_instances']} & {d['n_ranks']} & "
                    f"{fmt(d['beta'])} & {fmt(d['boot_sigma'])} & {fmt(d['fit_sigma'])} & "
                    f"{fmt(d['quoted_sigma'])}{'$^{a}$' if d['quoted_sigma_kind'] == 'least squares' else ''} & "
                    f"{tex_p(evr['p_perm_bonf'], evr['p_perm_floor'] * evr['bonf'])} \\\\")
    (GENERATED / "branch_table.tex").write_text("\n".join(rows) + "\n")

    keys = [d["branch"] for d in main["table"]]
    pw = {(t["labels"][0], t["labels"][1]): t for t in main["tests"] if t["family"] == "pairwise"}
    rows = []
    for a in keys:
        cells = []
        for b in keys:
            t = pw.get((a, b)) or pw.get((b, a))
            cells.append("--" if a == b else tex_p(t["p_perm"], t["p_perm_floor"]))
        rows.append(f"{alpha[a] if a != 'Waves+optics' else 'Waves'} & " + " & ".join(cells) + " \\\\")
    (GENERATED / "pairwise_table.tex").write_text("\n".join(rows) + "\n")

    rows = []
    for s in res["sensitivity"]:
        bmap = {d["branch"]: d for d in s["table"]}
        g = s["tests"][0]
        cells = [f"{fmt(bmap[k]['beta'])}" if k in bmap else "--" for k in keys]
        rows.append(f"{s['label']} & {s['n_total']} & " + " & ".join(cells) +
                    f" & {g['LR']:.1f} & {tex_p(g['p_perm'], g['p_perm_floor'])} \\\\")
    (GENERATED / "sensitivity_table.tex").write_text("\n".join(rows) + "\n")

    rows = []
    for c in res["corpora"]:
        rows.append(f"{c['label']} & {c['n']} & {c['n_instances']} & {fmt(c['beta'])} & {fmt(c['boot_sigma'])} & "
                    f"{fmt(c['fit_sigma'])} \\\\")
    (GENERATED / "corpus_table.tex").write_text("\n".join(rows) + "\n")


def tex_sci(p):
    e = int(np.floor(np.log10(p)))
    return f"{p / 10 ** e:.1f}\\times 10^{{{e}}}"


def tex_p(p, floor):
    return f"$<{floor:.4f}$".replace("<0.0", "<0.0") if p <= floor + 1e-12 else (f"{p:.4f}" if p < 0.01 else f"{p:.2f}")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true", help="200 permutations, no power simulation")
    args = ap.parse_args()
    reps = 200 if args.quick else 4000
    t0 = time.time()
    rng = np.random.default_rng(SEED)
    rows, dedupe, added = feynman_corpus()
    recs = [r for r, _, _ in rows]
    C = Counts(recs)
    print(f"Feynman corpus: {len(rows)} expressions (benchmark 100 + lectures {added}); removed {dedupe}")
    raw = Counter(b or "(unlabelled)" for _, b, _ in rows)
    print("labels:", dict(raw.most_common()))

    groups = groups_of(rows)
    table, tests = analyse(rows, C, groups, reps, rng, "main: fine branches")
    labelled = np.concatenate(list(groups.values()))
    pooled_all, pooled_lab = fit(recs), fit([recs[i] for i in labelled])
    print(f"  pooled all n={pooled_all['n']} β={pooled_all['beta']:.3f}±{pooled_all['boot_sigma']:.3f}; "
          f"branch-labelled n={pooled_lab['n']} β={pooled_lab['beta']:.3f}±{pooled_lab['boot_sigma']:.3f}")
    coarse = {}
    for k, idx in groups.items():
        coarse.setdefault(COARSE[k], []).append(idx)
    coarse = {k: np.concatenate(v) for k, v in coarse.items()}
    print("\n== coarse split ==")
    ck = list(coarse)
    coarse_test = test(C, "coarse: all", [coarse[k] for k in ck], ck, reps, rng, family="coarse")

    sens = []
    scenarios = [
        ("Main assignment", None, None),
        ("Ambiguous chapters dropped", None, lambda m, b: m["ambiguous"]),
        ("Ambiguous chapters reassigned", lambda m, b: (ALT.get((m["volume"], m["chapter"]), b)
                                                        if m["ambiguous"] else b),
         lambda m, b: m["ambiguous"] and (m["volume"], m["chapter"]) in ALT and ALT[(m["volume"], m["chapter"])] is None),
        ("Bracket/imaginary-unit expressions dropped", None, lambda m, b: m["notation"]),
        ("Vol. III dropped", None, lambda m, b: m["volume"] == "III"),
    ]
    for label, relabel, drop in scenarios[1:]:
        g = groups_of(rows, relabel, drop)
        t, ts = analyse(rows, C, g, reps // 2, rng, f"sensitivity: {label}", full=False)
        sens.append(dict(label=label, n_total=int(sum(len(v) for v in g.values())), table=t, tests=ts))
    main_sens = dict(label=scenarios[0][0], n_total=int(len(labelled)), table=table, tests=tests)
    sens.insert(0, main_sens)

    mdd, sims = ([], []) if args.quick else power(rows, C, groups, table, rng)
    corp = corpora(rows)
    ct = corpus_test(reps, rng)
    rf = {d["branch"]: rank_frequency(C, groups[d["branch"]], d["beta"]) for d in table}
    nlect = sum(1 for _, _, m in rows if m["source"] != "benchmark")
    counts = dict(nFeynAllExpr=len(rows), nLectAdded=nlect, nAmbiguousKept=sum(
        1 for i in labelled if rows[i][2]["ambiguous"]), nNotation=sum(1 for i in labelled if rows[i][2]["notation"]))
    res = dict(seed=SEED, bootstrap_reps=BOOT, perm_reps=reps, quick=args.quick,
               corpus=dict(n=len(rows), dedupe=dedupe, added=added, labels=dict(raw)),
               main=dict(table=table, tests=tests), coarse_test=coarse_test,
               pooled_all=pooled_all, pooled_labelled=pooled_lab, sensitivity=sens,
               power=dict(mdd80=mdd, simulation=sims), corpora=corp, corpus_test=ct,
               rank_frequency=rf, counts=counts, alt_assignment={f"{v}.{c}": b for (v, c), b in ALT.items()},
               runtime_s=time.time() - t0)
    RESULTS.mkdir(exist_ok=True)
    suffix = "_quick" if args.quick else ""
    (RESULTS / f"branch_test{suffix}.json").write_text(json.dumps(res, indent=1, default=float) + "\n")
    if not args.quick:
        with open(RESULTS / "branch_beta.csv", "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["branch", "n_expr", "n_ops", "n_ranks", "beta_logls", "boot_sigma", "boot_p16", "boot_p84",
                        "fit_sigma", "quoted_sigma", "quoted_sigma_kind"])
            for d in table:
                w.writerow([d["branch"], d["n"], d["n_instances"], d["n_ranks"], *(f"{d[k]:.4f}" for k in (
                    "beta", "boot_sigma", "boot_p16", "boot_p84", "fit_sigma", "quoted_sigma")), d["quoted_sigma_kind"]])
        with open(RESULTS / "branch_tests.csv", "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["family", "test", "n", "beta_sep", "beta_shared", "N_ops", "LR", "df", "p_chi2", "p_perm",
                        "perm_reps", "bonferroni", "p_perm_bonf", "dAIC", "dBIC"])
            for t in tests + [coarse_test, ct]:
                w.writerow([t["family"], t["name"], "|".join(map(str, t["n"])), "|".join(f"{b:.4f}" for b in t["beta_sep"]),
                            f"{t['beta_shared']:.4f}", t["N"], f"{t['LR']:.2f}", t["df"], f"{t['p_chi2']:.3g}",
                            f"{t['p_perm']:.5f}", t["perm_reps"], t["bonf"], f"{t['p_perm_bonf']:.5f}",
                            f"{t['dAIC']:.1f}", f"{t['dBIC']:.1f}"])
        write_numbers_tex(res)
    print(f"\ndone in {time.time() - t0:.0f} s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
