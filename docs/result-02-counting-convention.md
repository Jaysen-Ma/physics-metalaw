# Result 2 — part of `beta` is a fact about the counting rule, not about physics

Date: 2026-08-31. Code: `conventions.py`, `export_corpus.py`, `build_viewer.py`.
Data: the authors' released corpus, `github.com/andrei-const/physics-formulae-patterns`
(CC0, cloned to `.scratch/`, gitignored).

Result 1 established that `beta` depends on expression complexity. This one
establishes that it also depends on the *counting convention* — which is a human
choice the meta-law paper never states, and which had to be reverse-engineered
from their released data before any of Result 1 could be computed.

## The convention, re-derived and tested

Two rules in `metalaw.py` were inferred rather than read off the paper. Both are
now confirmed against the authors' own per-formula columns:

* **Symbolic names are counted as instances, not distinct names.** In
  `exp(-(theta/sigma)**2/2)/(sqrt(2*pi)*sigma)` the `sigma` appearing twice is
  charged twice: `variables = 3`, not 2. Holds on 100/100 Feynman, 65/70
  Inflationaris, 40/41 Wikipedia.
* **An integer exponent of 2, 3 or 4 collapses into the operator and the literal
  vanishes.** The `numbers` column matches the literals present only *after* the
  collapse: 68/70 Inflationaris against 0/70 without it, 97/100 Feynman against
  56/100, 40/41 Wikipedia against 22/41. The rule is unambiguously theirs.

There is no node for "the number `sqrt(2)`" and none for "2*pi". Every numeric
idea is built from the same leaves and operators as everything else:
`2*pi` is 3 nodes, `sqrt(2)` is 2, `sqrt(2/3)` is 4, `x**5` is 3, `x**2` is 2.

## What the ranked list actually contains

Every variable in a corpus collapses into a single rank entry, and so does every
plain number and every named constant. The Inflationaris ladder `beta` is fitted
to begins: `variable` 362, `Mult` 220, `number` 168, `Div` 167, `Pow4` 81,
`constant` 70. Three of the top six ranks are bucket categories, not operators.
Nullary nodes are 44-49% of every corpus.

## Data quality of the released corpora

| corpus | rows | won't parse | complexity matches theirs | published rank list |
|---|---|---|---|---|
| Feynman | 100 | 0 | 97/100 | 19 ranks, reproduced (1173 vs 1170 nodes) |
| Wikipedia | 41 | 0 | 38/41 | 15 ranks, reproduced **exactly** (550 nodes) |
| Inflationaris | 70 | **6** | 59/70 | 23 published vs **25** recounted |

Two things worth carrying forward:

* **Wikipedia is complete relative to the published analysis.** Result 1 left
  open whether the 41 released formulas are the corpus behind the published
  Wikipedia `beta`. They reproduce `wiki_data.txt` entry for entry — 15 ranks,
  550 instances. Nothing is missing.
* **Inflationaris has an unexplained rank discrepancy.** Six entries have
  unbalanced brackets and do not parse (repairs in `metalaw.EI_BRACKET_REPAIRS`).
  After repair the recount yields 25 distinct operator types against their 23,
  with rank 7 at 57 against their 59. Not yet diagnosed. Since `beta` is fitted
  to exactly this ladder, it should be diagnosed before leaning harder on the
  Inflationaris numbers.

## The variable/constant split is per formula, not per symbol

The CSVs publish only *counts* of variables and constants, never which symbol
carries which role — and the split is a physics judgement made per row, not a
property of the symbol. Eight Feynman formulas declare a constant although every
symbol in them appears as a variable elsewhere in the same corpus; in
`I.15.3x`, `(x-u*t)/sqrt(1-u**2/c**2)`, the `c` is the constant.

The split is therefore **not recoverable** from the released data. `metalaw`
takes the authors' three columns as given, and `export_corpus.py` marks such
symbols `unknown` rather than guessing.

## The operator sets are not consistent across corpora

| corpus | ranks | only in this corpus | shared with one other |
|---|---|---|---|
| Inflationaris | 25 | `abs arctan cosh coth sech sinh tan` | `cos log sin tanh` |
| Feynman | 19 | `arcsin ln` | `cos sin tanh` |
| Wikipedia | 15 | — | `log` |

Only 11 operators plus the 3 nullary buckets are shared by all three. And
**Feynman writes `ln` where the other two write `log`** for the same natural
logarithm (`I.44.4`, `n*kb*T*ln(V2/V1)`), so that function occupies a different
rank depending on the corpus. It is one instance and moves no published number,
but the vocabulary was never normalised across corpora.

This matters beyond tidiness: the normalised exponential runs over *the ranks
present*, so a corpus with a richer operator vocabulary has a longer tail to fit
through. Inflationaris has 25 ranks against Feynman's 19, most of the excess
being singletons. Part of its low `beta` is vocabulary breadth, which sits
alongside the complexity confound of Result 1 rather than replacing it.

## Folding the power family

`sqrt(x)` is `Pow(x, 0.5)`. Folding the whole family into one binary operator
with its exponent counted gives every power a uniform cost of 3:

| | authors' | folded |
|---|---|---|
| `sqrt(x)` | `sqrt + x` (2) | `Pow + x + 0.5` (3) |
| `x**2` | `Pow2 + x` (2) | `Pow + x + 2` (3) |
| `x**5` | `Pow + x + 5` (3) | `Pow + x + 5` (3) |

Taking the exponent of `sqrt` as the literal `0.5` rather than as `1/2` (which
would add a `Div` and two literals, making it 5 nodes) is a choice; it is the one
that keeps every power at exactly 3. `ln` is merged into `log` at the same time.

| convention | Inflationaris | Feynman | Wikipedia | gap (F − EI) |
|---|---|---|---|---|
| authors' set | 0.282 [0.270, 0.296] | 0.352 [0.336, 0.368] | 0.338 [0.316, 0.368] | **+0.071** |
| folded `Pow` | 0.289 [0.282, 0.296] | 0.324 [0.308, 0.337] | 0.319 [0.307, 0.334] | **+0.035** |
| + `exp` folded | 0.284 [0.276, 0.291] | 0.319 [0.303, 0.331] | 0.305 [0.291, 0.320] | +0.035 |

68% intervals from 2000 bootstrap resamples over formulas, seeded.

**Re-bucketing alone barely moves `beta`; counting the exponent moves it a lot,
and unevenly.** Merging five power ranks into one shifts every corpus by less
than its error bar. But restoring the exponent as a node shifts Inflationaris by
only +0.007 while Feynman falls 0.352 → 0.324, so the Feynman-minus-Inflationaris
gap the meta-law paper builds its discussion on **halves**. The intervals still
separate, so the gap is real — but half of its size is a consequence of the
`PowN` rule.

The asymmetry has a mechanical cause: merging only redistributes instances
between ranks and leaves the totals alone, whereas counting the exponent creates
86 new `number` nodes in Feynman, which pushes `number` up the ladder and
changes the totals. Folding also lengthens every formula that uses a power —
median complexity goes 17 → 19 (Inflationaris), 11 → 12.5 (Feynman),
10 → 11 (Wikipedia).

## Result 1 survives, with one change

Binning each convention by its own complexity:

| window | EI (authors') | EI (folded) | Feynman (authors') | Feynman (folded) |
|---|---|---|---|---|
| 3-7 | (n=2) | (n=0) | 0.545 | 0.519 |
| 8-11 | 0.329 | 0.357 | 0.413 | 0.469 |
| 12-17 | 0.321 | 0.298 | 0.303 | 0.286 |
| 18-25 | 0.296 | 0.293 | (n=7) | 0.268 |
| 26+ | 0.259 | 0.291 | (n=4) | (n=4) |

`beta` still falls steeply with complexity, and the corpora still agree where
they overlap — at 12-17 under folding the three read 0.298 / 0.286 / 0.317. What
changes is that **Inflationaris goes nearly flat above complexity 12** (0.298,
0.293, 0.291) where under the authors' convention it kept declining. So the one
claim Result 1 left as genuinely Inflationaris-specific — that its `beta(k)` is
flatter than the others — is itself convention-dependent, and gets *stronger*
here rather than weaker.

## Consequences for this project

1. **Never report a `beta` without its convention.** A number that moves by 0.028
   depending on whether an exponent is charged is not a bare constant.
2. **The null model must use the same convention as the measurement.** The
   paper's generator is a random walk over `+1` binary / `0` unary / `-1` nullary
   steps. Whether `square` is a unary primitive (one node, one step) or a binary
   `Pow` (three nodes) changes which step it occupies, so it changes the
   achievable arity mix at fixed `k` and therefore `beta_null(k)` itself. A
   unary `square` is cheap: it lets a random tree buy a nonlinearity for one
   node. Running the null model under one convention and comparing against
   `beta(k)` measured under the other would be an unforced error.
3. **Keep the authors' convention as primary** so comparison against their
   published numbers stays valid, and carry the folded convention as a standing
   robustness axis. `conventions.py` reports both.
4. **The Inflationaris 25-vs-23 rank gap is open.** Diagnose before relying on
   Inflationaris-specific conclusions.
