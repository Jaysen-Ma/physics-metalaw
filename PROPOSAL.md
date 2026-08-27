# Is the "meta-law of nature" a fact about nature, or a fact about physicists?

Research proposal. Status: design stage, no experiments run yet.
Session handoff and current status live in `.notes/STATE.md` (gitignored).

---

## 1. Background

Count how often each mathematical operator appears across a corpus of physics
equations, rank operators from most to least frequent, and the usage decays
exponentially with rank: `f(r) ~ exp(-beta * r)`. Constantin, Bartlett, Desmond
and Ferreira measured `beta` across four corpora and found the Feynman Lectures,
Wikipedia's named equations and the Cambridge Handbook all decay at roughly the
same rate, about `exp(-r/3)`. They proposed this as a candidate "meta-law of
nature": a statistical regularity in how physical laws are *written*.

One corpus came out different. **Encyclopaedia Inflationaris** (EI), the
catalogue of ~75 inflationary potentials proposed over forty years, decays more
slowly, about `exp(-r/4)`. The authors listed several candidate explanations but
could not determine which was responsible.

Published values, **verified 2026-08-27** against arXiv:2408.11065 and
independently reproduced from the authors' released data (see
`docs/result-01-beta-vs-complexity.md`):

| corpus | n | published beta | reproduced |
|---|---|---|---|
| Feynman Lectures | 100 | 0.342 +/- 0.010 | 0.352 +/- 0.011 |
| Wikipedia named equations | 41 released | 0.35 +/- 0.03 | 0.338 +/- 0.016 |
| Cambridge Handbook | 852 | 0.33 +/- 0.04 | data not released |
| Encyclopaedia Inflationaris | 70 | 0.266 +/- 0.014 | 0.282 +/- 0.011 |

Unary shares are also confirmed: 0.177 in Inflationaris (paper: "no more than
1/5"), 0.126 in Feynman and 0.096 in Wikipedia (paper: "under 1/7").

The five explanations the authors offer for Inflationaris being different, in
their own terms: all 71 expressions describe the same physical quantity; they
involve a single variable; natural units suppress multiplication and division;
that in turn raises the unary share; and the corpus is speculative, which "may
allow for a broader variety of operator combinations". Only the last is a
property of the *people* rather than of the *problem* — the earlier framing of
this proposal attributed richer sociological claims to them than they make.
They do **not** compare against any machine-generated catalogue, and they do
**not** investigate complexity dependence.

## 2. Research question

**If inflationary potentials are generated bottom-up — by searching expression
space and keeping only what is physically viable, with no human input about what
a potential "should" look like — do the operator statistics match the
human-written catalogue, or differ from it?**

This is a clean experiment because the original authors' candidate explanations
split into two disjoint groups:

- **Properties of the problem.** Single field variable, natural units, an
  energy-scale prefactor. A machine-generated catalogue shares all of these.
- **Properties of the people.** Familiar functional forms, expressions that
  communicate well, new models built by modifying old ones. A machine-generated
  catalogue shares none of these.

Comparing the two catalogues therefore isolates the human contribution.

## 3. Hypothesis

Bottom-up viable potentials will be dominated by **nested saturating unary
functions** (`tanh`, `arctan`, `cosh`, ...), giving a unary share far above the
<20% seen in EI, and a correspondingly different `beta`.

**Mechanism.** Inflation needs a flat region for the field to roll slowly.
Nesting saturating functions produces a plateau essentially for free — that is
what those functions do. The human-style alternative, polynomials, produces
steep power-law potentials, which is exactly what the bound on the
tensor-to-scalar ratio (`r < 0.061`) has been ruling out. So humans and blind
search solve the *same* problem — build a plateau — by structurally opposite
means. Theorists write compact structured forms like Starobinsky's
`(1 - exp(-sqrt(2/3) phi))^2`; search stacks saturating functions.

### 3.1 Status of the supporting evidence — READ THIS

The "300+ viable potentials at complexity 5, ~80% unary-nested, e.g.
`cosh(arctan(cosh(tanh(phi))))`" result comes from **subsequent unpublished work
(c. 2023), not from the thesis** — confirmed by the author on 2026-08-15.

**The code and data for it no longer exist.** Searched on 2026-08-15: no
inflation-related source anywhere under `/home/zotar` outside this repo, no
archives, and on GitHub the account cited in the thesis
(`Coding-appreciator/ESR`) now 404s while the current account holds no
predecessor repo. The result is therefore a **recollection, not a citable
dataset**. It cannot be reproduced, error-barred, or checked for the
selection artefacts listed in `docs/point-cloud-sr-plan.md`.

**Consequences.** (a) Treat the hypothesis as a plausible untested prior with a
good mechanism — do not write it as prior data, and do not cite a number.
(b) Step 0 is the decisive gate: the belief that viable potentials are common
now rests on nothing measured. (c) Regenerating that result is itself a
deliverable, and doing it from scratch is the right call — it needs the error
bars, controls and reproducibility the original lacked.

For the record, what the **thesis** contains — checked against
`../ISR_latex/main.tex`, and none of it is the above:

| Thesis result | What it is |
|---|---|
| `ats(phi) = c1*arctan(c2*sin(c3*phi))` | One function class, from GA crossover of natural + arctan inflation |
| 1768 solutions | Viable **coefficient tuples** inside `ats`, in one log-region of its 3D parameter space |
| 599 viable potentials | Viable **coefficient tuples** inside `c1*|arctan(c2*phi)|^c3`, found by GA over Table 5 ranges |
| ~4000 function classes | ESR enumeration at max complexity 6 over basis `{x,a ; inv,sin,arctan ; +,-,*,/,pow}`, fitted to *data points* from one viable curve by MDL — **never screened for viability** |
| Sec 4.5 "Function Relatives" | Symbolic degeneracies: `sin`->`cos` by axis relabel, absorbable additive constants, parity flips |

The thesis never enumerated viable *function classes* and never reported an
operator mix. `tanh` and `cosh` appear nowhere in it. Its large counts are
parameter-space densities *within an already-known-good class*, which says
nothing about the density of viable classes in expression space. What it does
contribute is the fitness criterion (reimplemented in `inflation_fitness.py`)
and the degeneracy taxonomy needed for Step 4.

### 3.2 Prior evidence that viability is RARE

The ACHL paper (source in `../ciga_tex/InflationGA.tex`, our fitness's origin)
reports the opposite of the "viability is common" assumption, and it is the best
quantitative evidence available before Step 0 runs:

- "a typical random selection of 5000 models would contain **no viable
  potentials at all**" (`InflationGA.tex:523`); a 50-generation GA run visiting
  ~5000 states yields only `O(10)` viable models.
- The sextic-polynomial dataset — 309,097 viable models — required **70,000 GA
  runs** visiting `~1e8` states, about 1% of the search space
  (`InflationGA.tex:525`). Acceptance is on the order of `1e-4` *with* GA
  guidance, lower without.

**This is about coefficient space, not expression space** — the two questions
differ, and the distinction drives the cost model. Our viability predicate is
"does there **exist** a parameter setting making this expression viable", which
is far more permissive than a random coefficient draw. But it means the per
expression cost is an **inner parameter optimisation**, and ACHL demonstrates
that random coefficient sampling is inadequate as that inner search. Step 0 must
therefore measure two things, not one: the fraction of expressions that are
viable, *and* the cost of deciding it. The second sets the whole budget.

Note also that the `~300,000` in ACHL's abstract is a plausible corruption path
for the remembered "300+" figure. Another reason to re-derive rather than trust.

## 4. Methodology

Complexity `k` = node count of the expression tree (ESR convention; matches the
thesis definition at `main.tex:112`).

### Step 0 — Viability density vs complexity

Generate random expressions at several `k` and measure the fraction that are
viable. Produce a **curve**, not a single number: it sizes the sampling budget
per bin downstream, and "density of viable inflation models in expression space
as a function of complexity" is a standalone result.

Decision gate: if viable potentials are reasonably common, plain rejection
sampling works and is statistically clean. If they are rare, directed search is
needed, which injects bias that then has to be corrected.

### Step 1 — Low-complexity catalogue, exhaustive

Use the ESR code to enumerate every expression up to `k = 8` over a chosen
operator basis. For each, optimise over free parameters to decide whether any
setting is viable: enough e-folds, and `n_s`, `r`, `m_star` within observational
bounds. Complete and unbiased at low complexity.

### Step 2 — High-complexity catalogue, by sampling

Exhaustive enumeration dies above `k ~ 8`. Randomly generate expressions up to
`k ~ 30` and filter for viability. **Use the same random generator the original
paper used for its null model**, so results are directly comparable.

Sample at *fixed* `k` per bin (the paper's generator samples uniform trees of
exact size, so this is free): exhaustive for `k <= 8`, then bins at
`k ~ 12, 16, 20, 25, 30`. Rejection-sample each bin until enough survivors for
stable statistics.

### Step 3 — The matched control

Run the same generator twice: once keeping everything, once keeping only viable
potentials. Two catalogues identical in operator basis, complexity distribution
and generator, differing *only* by the physics filter. Any statistical
difference is caused by the physics and nothing else. This is what makes the
result defensible, and it is also what neutralises the operator-vocabulary
caveat (Section 6) — so it is not optional.

### Step 4 — Deduplication

Many expressions differ symbolically but describe the same physical model. The
thesis degeneracy taxonomy (Sec 4.5) covers axis relabelling, absorbable
additive constants and parity flips. The robust approach is to fingerprint each
candidate by its **inflationary behaviour** (`n_s`, `r`, `m_star`, e-fold count,
plateau geometry) rather than symbolic form, then cluster.

### Step 5 — Measure and compare

Report the **operator-arity mix (nullary/unary/binary) first, `beta` second**.
When 80% of operators are unary, `beta` summarises a distribution so
structurally unlike the human ones that the single number hides the story.

Report everything **per complexity bin against the matched control**. Bootstrap
CIs per bin are mandatory — a ~10-operator basis gives only ~10 ranks to fit an
exponential through, so per-bin fits will be noisy. Do not over-interpret
wiggles.

## 5. Complexity dependence of beta

This is a first-class question, not an add-on: **is `beta` a constant, or a
function of complexity?**

### 5.1 Why beta almost certainly moves with k even with no physics

Tree combinatorics constrains the achievable arity mix at fixed `k`. With `b`
binary nodes, `u` unary nodes and `l` leaves: `l = b + 1` and `k = 2b + 1 + u`.
So the reachable nullary/unary/binary shares shift mechanically as `k` grows,
and the *unfiltered* generator's `beta` already drifts with `k`.

The quantity of interest is therefore not `beta(viable, k)` but the **gap**
`beta(viable, k) - beta(control, k)` per bin. Any publishable claim lives in
that gap.

### 5.2 The aggregation confound in EI's published beta — MEASURED, and it is large

**This is no longer a conjecture. It was measured on 2026-08-27 and it accounts
for most of EI's apparent anomaly.** Full write-up and caveats in
`docs/result-01-beta-vs-complexity.md`; the short version:

`beta` falls from about 0.54 at complexity 3-7 to about 0.26 at complexity 26+,
in every corpus. That span is roughly four times the EI-versus-Feynman gap the
paper is built on. In the complexity windows where the corpora overlap they
agree — at complexity 12-17 the 68% intervals are EI 0.321 [0.296, 0.345],
Feynman 0.303 [0.283, 0.320], Wikipedia 0.343 [0.306, 0.386]. Median
complexities are EI 17, Feynman 11, Wikipedia 10, and only EI populates
complexity 18 and above, which is where the low values sit.

So EI does not decay more slowly because it is cosmology. It decays more slowly
because its expressions are longer. What survives as genuinely EI-specific is
the *slope*: `beta(k)` is much flatter in EI than in Feynman or Wikipedia.

Consequence for this project: never compare a generated catalogue against the
aggregate 0.266. Compare window by window.

The original argument, which still holds:

An aggregate `beta` conflates two things: per-complexity operator usage, and the
catalogue's complexity *distribution*. EI's `beta` is a weighted average over
wherever its ~75 models happen to sit in complexity. Part of "EI decays slower"
could simply be that human models occupy a different complexity range than the
other corpora.

Stratifying fixes this and buys an extra comparison: compute the complexities of
the EI potentials under the same counting convention, then compare EI's `beta`
against the viable catalogue **reweighted to EI's complexity histogram**. If
they still disagree there, the human-fingerprint claim is much stronger.

### 5.3 Binning must use simplified complexity

A generated complexity-25 expression frequently simplifies to complexity 7
wearing dead weight. Binning by *generated* `k` pollutes high-`k` bins with
low-`k` physics plus inert decoration, flattening `beta(k)` artificially. **Bin
by the complexity of the simplified, deduplicated representative.** Keep raw
binning as a robustness check — the discrepancy is itself a result about
redundancy in expression space.

### 5.4 Two competing mechanisms — which wins is the finding

1. **Nesting demand.** Plateaus require a minimum depth of saturating-function
   nesting, so at low `k` the filter is maximally selective and the unary excess
   is largest.
2. **Dilution.** At high `k`, most of a viable expression's tree is physically
   inert over the observable field range (a few e-folds), so survivor statistics
   drift back toward the generator prior:
   `beta(viable, k) - beta(control, k) -> 0` as `k -> infinity`.

If dilution wins at large `k`, then "the beta of viable inflationary potentials"
is not a number at all — it is a function of complexity with the physics signal
concentrated at low `k`. That is itself an argument against reading the meta-law
as a universal constant, and it says *where* to compare against humans: at the
`k`-range EI actually occupies. If the gap survives there, that is Outcome B
with teeth.

**Cost note.** Stratification multiplies compute by the number of bins, and
high-`k` bins are the expensive ones (likely lower viability rate, plus costlier
parameter fits per candidate). Decide after Step 0 whether `k = 30` is
affordable or whether the curve should stop at 20-25.

## 6. Outcomes and their meaning

| | Result | Interpretation |
|---|---|---|
| **A** | Bottom-up statistics match EI | EI's deviation is forced by the structure of the problem, not human habit. Supports an intrinsic reading of the meta-law, but means `beta` is domain-dependent, mildly undercutting universality. |
| **B** | Bottom-up differs from EI, matches general physics | The domain does not force EI's signature, so that signature is a measurable fingerprint of human theory-building. |
| **C** | Bottom-up differs from everything | The exponential regularity is a property of *published, curated* equations rather than of physically-valid expressions in general. Most direct challenge to the original claim. |

**Expected:** a version of B, but stronger than a shifted decay rate — a
wholesale difference in the *kind* of expression that survives, with plateau
construction as the concrete mechanism. Held as a prior with a mechanism, not as
a replication (see Section 3.1).

## 7. Caveats to build in from the start

1. **Operator vocabulary partly determines the answer.** Hand the generator many
   saturating functions and it returns many saturating functions. The unfiltered
   control corpus (Step 3) is what neutralises this.
2. **The amplitude prefactor — RESOLVED for the primary convention.** Every
   viable potential carries a prefactor of order `1e-10`. The authors **do**
   count it: every EI entry in their corpus is written `(M**4)*(...)`, and the
   per-formula complexity they publish includes `M` + `Pow4` + `Mult` = 3 nodes
   for it (verified 2026-08-27 on 59/70 formulas). So the primary convention
   counts the prefactor as 3 nodes. Still **report both ways**, since it is a
   constant additive contamination of every expression in the catalogue.
3. **The fitness function is one specific viability criterion**, from Abel,
   Constantin, Harvey & Lukas (arXiv:2208.13804), with documented deliberate
   deviations (see the `inflation_fitness.py` module docstring). Results are
   conditional on it; a sensitivity check on the cut values is worth budgeting.

## 8. Existing assets

| Path | What it is |
|---|---|
| `inflation_fitness.py` | Clean reimplementation of the ACHL fitness. `evaluate()` / `fitness()` on a `Potential`; `from_sympy()` bridges symbolic input; all cuts in `FitnessConfig`. Tracked, no known issues. |
| `fitness.py` | Thesis-era original script. Untracked, superseded, kept only for reference. |
| `../ISR_latex/main.tex` | The MSc thesis. Sec 4.5 = degeneracy taxonomy for Step 4. |
| `../ciga_tex/InflationGA.tex` | **Full LaTeX source of the ACHL paper** (arXiv:2208.13804), the origin of our fitness function — a pre-publication draft with author comments still in it. Authoritative for exact cut definitions, state-space sizes, GA settings and runtimes. Use it instead of the arXiv PDF when reproducing. |
| `docs/point-cloud-sr-plan.md` | A prior full research plan for an **adjacent project** (see Section 10). Rescued from `~/.claude/plans/` on 2026-08-15. |
| `metalaw.py` | Operator counting under the paper's convention, reverse-engineered and validated against their published complexity column and rank list; plus the normalised-exponential `beta` fit. |
| `analyse_corpora.py` | Complexity-stratified `beta` for the three released corpora, with bootstrap intervals and two artifact checks. Produces `docs/result-01-beta-vs-complexity.md`'s numbers. |
| `docs/result-01-beta-vs-complexity.md` | **Result 1.** `beta` is not a corpus constant; it falls steeply with complexity, and the corpora agree once matched. |
| `.scratch/physics-formulae-patterns/` | Clone of the authors' released code and data (gitignored). Contains all three corpora with per-formula complexity, the published rank-frequency lists, and their Mathematica random-expression generator. |

## 10. Adjacent plan: free-form point clouds -> symbolic regression

`docs/point-cloud-sr-plan.md` (written 2026-08-05, then lost outside the repo)
proposes a different attack on the same field: represent the potential
non-parametrically, evolve it until it inflates, and *then* symbolically regress
the surviving curves — so functional form is an output rather than an input.

Three things in it bear directly on this proposal and should be used regardless
of whether that project is pursued:

1. **The `w = V'/V` change of variable.** Every slow-roll quantity becomes at
   most one derivative from the sampled object (`eps = w^2/2`, `n_s = 1 - w^2 + 2w'`,
   `r = 8w^2`), and the inflating region is nearly geometric: `|w| < 1`. This is
   the single best lead for making the inner viability decision cheap — the cost
   driver identified in Section 3.2. Caveat for *this* project: EI catalogues
   `V(phi)`, so operator statistics must be measured on `V`, not `w`. Sampling
   expressions in `w`-space to raise acceptance would change what is being
   counted, and is not a free win.
2. **The amplitude direction is analytic, not searched.** `n_s` and `r` are
   invariant under `V -> lambda V`, so the overall scale is fixed exactly by
   `A_s = V/(24 pi^2 eps)` rather than searched. This removes a whole dimension
   from any inner optimisation, and it is a genuine criticism of ACHL's `f_m`
   term. It also interacts with Caveat 7.2: the prefactor is a *derived
   normalisation*, which strengthens the case for reporting operator counts both
   with and without it.
3. **The "must not be repeated" list** (end of that document) — no baseline,
   distributions reported over quantities that were cut on, a figure
   contradicting its caption, no error bars, no seeds. That list was written
   about the thesis and applies verbatim to this project. In particular:
   *never report an operator-frequency distribution without the matched
   unfiltered control*, which is exactly why Step 3 is non-negotiable.

## 11. References

- **arXiv:2408.11065** — Constantin, Bartlett, Desmond & Ferreira, *Statistical
  Patterns in the Equations of Physics and the Emergence of a Meta-Law of
  Nature* (2024). The paper under test. Companion code and data:
  https://github.com/andrei-const/physics-formulae-patterns (cloned to
  `.scratch/`). Releases all three corpora with per-formula complexity, the
  ranked operator counts, the plotting notebook, and a Mathematica notebook
  `RandomExpressions.nb` implementing the null generator.
- **arXiv:2310.16786** — Sousa, Bartlett, Desmond & Ferreira, *Optimal
  Inflationary Potentials* (2023). **Directly overlaps Step 1**: they use ESR to
  exhaustively enumerate scalar-field potentials to complexity ~8-9 over two
  operator bases — Set A `{x, a, inv, exp, log, sin, sqrt|.|, cube, square, +,
  *, -, /}` and Set B `{x, a, inv, exp, log, +, *, -, /, power}` — treat each as
  single-field slow roll, and score by minimum description length rather than a
  hard viability cut. Parameters fitted with BFGS from multiple starts; end of
  inflation at `epsilon = 1`; `N = 60` e-folds; `n_S = 0.9649 +/- 0.0042`,
  `r < 0.028` (95%), `A_S`. Data "on reasonable request" — not released.
  Three authors overlap with the meta-law paper, so the bottom-up catalogue
  existed in the same group and was never compared against Inflationaris.
- **arXiv:2208.13804** — Abel, Constantin, Harvey & Lukas, *Cosmic Inflation and
  Genetic Algorithms*. Source of `inflation_fitness.py`. Local LaTeX source at
  `../ciga_tex/InflationGA.tex`.
- **ESR** (Exhaustive Symbolic Regression), Bartlett, Desmond & Ferreira:
  https://github.com/DeaglanBartlett/ESR — reachable, not yet installed.
- **arXiv:2601.14288** — *identified 2026-08-27*: Peng, Yuan, Lai, Jiang, Ye,
  Zhang & Piao, *DeepInflation: an AI agent for research and model discovery of
  inflation* (Jan 2026, rev. Jun 2026). A multi-agent LLM system with a symbolic
  regression engine and a retrieval-augmented knowledge base that proposes
  viable single-field slow-roll potentials. Relevant as *contrast*, not as
  method: it is retrieval-augmented over the human literature, so its output
  inherits exactly the human prior this project is trying to measure. If it ever
  becomes a comparison arm, that is the point to make about it.
