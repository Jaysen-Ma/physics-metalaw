"""Complexity-stratified operator statistics for the three released corpora.

Answers: is beta a property of a corpus, or does it drift with the complexity
of the expressions in it? Run before generating anything -- if beta already
moves with complexity inside the human corpora, a single aggregate beta is not
the right object to compare a machine-generated catalogue against.
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
import metalaw as M

DATA = Path(__file__).parent / ".scratch/physics-formulae-patterns/data"
CORPORA = [("inflation", "Inflation.csv"), ("feynman", "Feynman.csv"),
           ("wiki", "Wiki.csv")]
RNG = np.random.default_rng(20260827)
N_BOOT = 2000


def bootstrap_beta(entries, n_boot=N_BOOT):
    """Resample formulas (not operator instances) with replacement."""
    betas = []
    idx = np.arange(len(entries))
    for _ in range(n_boot):
        sample = [entries[i] for i in RNG.choice(idx, size=len(idx))]
        try:
            betas.append(M.fit_beta(M.rank_counts(sample).values())[0])
        except RuntimeError:
            continue
    b = np.array(betas)
    return np.median(b), np.percentile(b, 16), np.percentile(b, 84)


def report(name, entries):
    counts = M.rank_counts(entries)
    beta, sigma = M.fit_beta(counts.values())
    lo_med, lo, hi = bootstrap_beta(entries)
    shares = M.arity_shares(entries)
    ks = np.array([e.complexity_declared for e in entries])
    print(f"  {name:14s} n={len(entries):3d}  k=[{ks.min():2d},{ks.max():3d}] "
          f"med {np.median(ks):5.1f}  ranks={len(counts):2d}  "
          f"beta={beta:.3f} (boot {lo:.3f}-{hi:.3f})  "
          f"unary={shares['unary']:.3f}  N={int(shares['n_instances'])}")
    return beta, lo, hi


for corpus, fname in CORPORA:
    entries = M.load_corpus(DATA / fname, corpus)
    print(f"\n=== {corpus} ===")
    ks = np.array([e.complexity_declared for e in entries])
    print(f"  checksum: sum(complexity)={sum(e.complexity_declared for e in entries)} "
          f"recomputed={sum(e.complexity_recomputed for e in entries)} "
          f"exact={sum(e.complexity_recomputed == e.complexity_declared for e in entries)}/{len(entries)}")
    report("ALL", entries)
    # terciles of complexity
    edges = np.quantile(ks, [0, 1 / 3, 2 / 3, 1.0])
    for j in range(3):
        lo_k, hi_k = edges[j], edges[j + 1]
        sel = [e for e in entries
               if (e.complexity_declared >= lo_k and
                   (e.complexity_declared < hi_k or j == 2))]
        if len(sel) < 8:
            continue
        report(f"k in [{lo_k:.0f},{hi_k:.0f}]", sel)


# ---------------------------------------------------------------------------
# Matched absolute complexity windows. This is the comparison that matters:
# terciles above are per-corpus, so they do NOT compare like with like.
# ---------------------------------------------------------------------------
WINDOWS = [(3, 8), (8, 12), (12, 18), (18, 26), (26, 80)]
MIN_PER_WINDOW = 8

print("\n\n=== beta in matched complexity windows (68% bootstrap intervals) ===")
loaded = {c: M.load_corpus(DATA / f, c) for c, f in CORPORA}
header = " | ".join(f"{c:^28s}" for c in loaded)
print(f"{'window':>10s} | {header}")
print(f"{'':>10s} | " + " | ".join(f"{'n':>3s} {'beta':>6s} {'68% CI':>15s}"
                                   for _ in loaded))
for lo, hi in WINDOWS:
    cells = []
    for corpus, entries in loaded.items():
        sel = [e for e in entries if lo <= e.complexity_declared < hi]
        if len(sel) < MIN_PER_WINDOW:
            cells.append(f"{len(sel):3d} {'--':>6s} {'':>15s}")
            continue
        beta, _ = M.fit_beta(M.rank_counts(sel).values())
        _, b_lo, b_hi = bootstrap_beta(sel, n_boot=3000)
        cells.append(f"{len(sel):3d} {beta:6.3f} {f'[{b_lo:.3f},{b_hi:.3f}]':>15s}")
    print(f"{f'[{lo},{hi})':>10s} | " + " | ".join(cells))

print("\nComplexity distributions (what an aggregate beta silently averages over):")
for corpus, entries in loaded.items():
    ks = np.array([e.complexity_declared for e in entries])
    q = np.percentile(ks, [25, 75]).round(1)
    print(f"  {corpus:10s} n={len(ks):3d}  median={np.median(ks):5.1f}  "
          f"mean={ks.mean():5.1f}  IQR={list(q)}")
