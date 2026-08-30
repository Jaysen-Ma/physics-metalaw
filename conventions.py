"""Alternative operator conventions, and how much of `beta` depends on them.

`metalaw.py` implements the meta-law paper's own convention and is validated
against their published per-formula columns; it must not drift. This module
recounts the same corpora under alternative conventions, so that the part of
`beta` which is a fact about the *counting rule* rather than about physics can
be measured instead of assumed.

The convention under test folds the whole power family into one binary operator
with its exponent counted as a node:

    authors'                        folded
    sqrt(x)   sqrt + x        (2)   Pow + x + 0.5   (3)
    x**2      Pow2 + x        (2)   Pow + x + 2     (3)
    x**5      Pow + x + 5     (3)   Pow + x + 5     (3)

Taking the exponent of `sqrt` as the single literal `0.5` rather than as `1/2`
(which would add a Div and two literals, making it 5 nodes) is a choice; it is
the one that keeps every power costing exactly 3.

`ln` is also merged into `log`. Feynman writes `ln` in I.44.4 while
Inflationaris and Wikipedia write `log` for the same function, so under the
authors' convention the natural logarithm occupies a different rank in
different corpora.

Run directly for the comparison report:  .venv/bin/python conventions.py
"""
from __future__ import annotations

from collections import Counter
from pathlib import Path

import numpy as np

import metalaw as M

DATA = Path(__file__).parent / ".scratch/physics-formulae-patterns/data"
CORPORA = [("inflation", "Inflation.csv"), ("feynman", "Feynman.csv"),
           ("wiki", "Wiki.csv")]

# Unary operators that are a power in disguise. Each costs one extra node once
# its exponent is written out.
UNARY_POWERS = ("Pow2", "Pow3", "Pow4", "sqrt")

CONVENTIONS = ("base", "fold", "fold_exp")


def recount(entry: M.EIEntry, convention: str = "base") -> tuple[Counter, int]:
    """Rank counts and complexity for one formula under `convention`.

    base      the authors' convention, unchanged (see `metalaw`)
    fold      every unary power becomes a binary Pow with its exponent counted
    fold_exp  additionally exp(x) -> Pow(e, x), with `e` as a named constant

    Nullary tallies come from the authors' own columns, as in `metalaw.load_*`;
    folding only ever *adds* to them.
    """
    if convention not in CONVENTIONS:
        raise ValueError(f"unknown convention {convention!r}; expected one of "
                         f"{CONVENTIONS}")
    ops = Counter(entry.operator_counts)
    n_numbers, n_constants = entry.n_numbers, entry.n_constants
    added = 0

    if convention != "base":
        folded = sum(ops.pop(k, 0) for k in UNARY_POWERS)
        ops["Pow"] += folded
        n_numbers += folded          # the exponent literal becomes a node
        added += folded
        if convention == "fold_exp":
            n_exp = ops.pop("exp", 0)
            ops["Pow"] += n_exp
            n_constants += n_exp     # the base `e` becomes a named constant
            added += n_exp
        ops["log"] += ops.pop("ln", 0)

    counts = Counter({k: v for k, v in ops.items() if v})
    counts[M.NULLARY_VARIABLE] += entry.n_variables
    counts[M.NULLARY_NUMBER] += n_numbers
    counts[M.NULLARY_CONSTANT] += n_constants
    return counts, entry.complexity_declared + added


def rank_counts(entries: list[M.EIEntry], convention: str = "base") -> Counter:
    total: Counter = Counter()
    for e in entries:
        total += recount(e, convention)[0]
    return total


def fit_beta(entries: list[M.EIEntry], convention: str = "base") -> float:
    return M.fit_beta(rank_counts(entries, convention).values())[0]


def bootstrap_beta(entries, convention="base", n_boot=2000, seed=20260830):
    """68% interval, resampling formulas (not operator instances)."""
    rng = np.random.default_rng(seed)
    idx = np.arange(len(entries))
    betas = []
    for _ in range(n_boot):
        sample = [entries[i] for i in rng.choice(idx, size=len(idx))]
        try:
            betas.append(fit_beta(sample, convention))
        except RuntimeError:
            continue
    b = np.array(betas)
    return np.percentile(b, 16), np.percentile(b, 84)


def _report() -> None:
    loaded = {c: M.load_corpus(DATA / f, c) for c, f in CORPORA}
    labels = {"base": "authors' set", "fold": "folded Pow",
              "fold_exp": "+ exp folded"}

    print("=== Operator sets are not the same across corpora ===\n")
    tot = {c: rank_counts(e) for c, e in loaded.items()}
    shared = set.intersection(*[set(t) for t in tot.values()])
    for c, t in tot.items():
        only = sorted(set(t) - shared)
        print(f"  {c:10s} {len(t):2d} ranks; unique to it: {', '.join(only) or '-'}")
    print(f"  shared by all three: {len(shared)}")
    print("  note: Feynman writes `ln`, the others `log`, for the same function.")

    print("\n\n=== Aggregate beta ===\n")
    print(f"  {'convention':14s}" + "".join(f"{c:>26s}" for c in loaded))
    beta = {}
    for conv in CONVENTIONS:
        cells = []
        for c, ents in loaded.items():
            b = fit_beta(ents, conv)
            lo, hi = bootstrap_beta(ents, conv)
            beta[(c, conv)] = b
            cells.append(f"{b:.3f} [{lo:.3f},{hi:.3f}] {len(rank_counts(ents, conv)):2d}rk")
        print(f"  {labels[conv]:14s}" + "".join(f"{x:>26s}" for x in cells))

    print("\n  Encyclopaedia Inflationaris barely moves; Feynman moves a lot, so")
    print("  the gap the meta-law paper is built on is convention-dependent:\n")
    for conv in CONVENTIONS:
        ei = beta[("inflation", conv)]
        print(f"    {labels[conv]:14s} EI={ei:.3f} ({ei - beta[('inflation','base')]:+.3f})"
              f"   Feynman-minus-EI gap = {beta[('feynman', conv)] - ei:+.3f}")

    print("\n\n=== Complexity-stratified, each binned by its own complexity ===")
    for conv in ("base", "fold"):
        print(f"\n  --- {labels[conv]} ---")
        print(f"  {'window':>9s} |" + "|".join(f"{c:^26s}" for c in loaded))
        for lo_k, hi_k in [(3, 8), (8, 12), (12, 18), (18, 26), (26, 90)]:
            cells = []
            for c, ents in loaded.items():
                sel = [e for e in ents if lo_k <= recount(e, conv)[1] < hi_k]
                if len(sel) < 8:
                    cells.append(f"  n={len(sel):<3d}{'--':>19s}")
                    continue
                b = fit_beta(sel, conv)
                lo, hi = bootstrap_beta(sel, conv, n_boot=1500)
                cells.append(f"  n={len(sel):<3d}{b:6.3f} [{lo:.3f},{hi:.3f}]")
            print(f"  {f'[{lo_k},{hi_k})':>9s} |" + "|".join(cells))


if __name__ == "__main__":
    _report()
