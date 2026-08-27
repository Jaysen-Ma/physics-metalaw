# Free-form point clouds → symbolic regression for inflationary potentials

---

# Part 1 — The plan in plain English

## What the old work did, and why it stalled

Inflation needs a potential — some curve `V(φ)` describing the energy of the inflaton field.
Nobody knows what that curve is. What we know is what it has to *produce*: about 60 e-folds of
expansion, a spectral tilt near the measured value, a small tensor-to-scalar ratio, and the right
overall amplitude.

The paper's approach was: **guess the shape of the curve, then hunt for numbers.** They assumed
the potential was a sextic polynomial (or a quartic plus a cosine, or a quartic plus an
exponential), chopped each coefficient into 64 bins, encoded the whole thing as a bit string, and
let a genetic algorithm breed bit strings until it found combinations that worked. It worked in
the sense that it found a lot of them — a few hundred thousand.

Your thesis pushed one step further by hand-building a *new* shape (`arctan` of a `sin`) out of two
old ones, then hunting for numbers inside that.

But both are doing the same thing. **The shape of the curve is always something you decide before
you start.** The search only ever adjusts numbers inside a shape you already committed to. So the
answer to "what functional forms allow inflation?" can never be anything other than the list you
typed in at the beginning. The paper knows this — its section on genetic programming says so
explicitly, then says they'll deal with it in a future paper. Your thesis says the same thing in
its future-work section and sketches the fix: replace part of the function with "point clouds".

That's the idea we're going to actually build.

## The new idea, in one paragraph

Stop assuming a shape. Represent the potential as **a bag of loose points** — just a set of heights
at a set of positions, joined up into a smooth curve, with no formula behind it at all. Wiggle
those points around until the curve produces real inflation. Do this many times, and you end up
with a large collection of curves that all work, none of which was ever forced into a polynomial
or a cosine or anything else.

*Then* bring in symbolic regression — but for a completely different job than before. Previously,
symbolic regression was a fancy way to propose a formula to go hunt inside. Now it's the last step:
you have a pile of curves that demonstrably work, and you ask the machine **"what formula describes
these?"** The functional form stops being an assumption you feed in and becomes a result you read out.

That inversion is the whole project.

## The one thing that makes this hard, and how we get around it

There's a catch, and it's the reason nobody just does this.

Two of the quantities you need are built from the *curvature* of the potential — how sharply it
bends. Curvature is a second derivative, and second derivatives are brutal on a curve drawn through
loose points. Nudge one point slightly and the curvature nearby swings wildly. The physics you're
trying to measure gets drowned in noise from your own interpolation. On top of that, if you scatter
points at random, almost none of the resulting curves will inflate at all — you'd spend all your
compute throwing candidates away.

The fix is to change *what the points represent*. Instead of the points being the height of the
potential, let them be its **steepness** — specifically, how fast the potential is changing relative
to how big it is. Call that `w`.

This sounds like a small bookkeeping change and it isn't. Three things happen at once:

- **Everything becomes one derivative shallower.** Every quantity we need is now either `w` itself or
  `w` differentiated once. The noise problem largely evaporates.
- **You can build working curves instead of searching for them.** Inflation happens exactly where the
  potential is shallow enough — which in this language is a simple, direct condition on `w`. So
  instead of generating random candidates and discarding 9,999 out of 10,000, you *construct*
  candidates that inflate by design. This is the difference between the project being feasible on
  your machine and not.
- **One whole dimension of the search disappears.** Scaling the entire potential up or down doesn't
  change the tilt or the tensor ratio at all — it only changes the overall amplitude, which is
  measured very precisely. So instead of searching for the right overall size, we just solve for it
  in one line. Free speedup, and it removes a source of bias.

And you still get an ordinary potential `V(φ)` out at the end — the thing you asked for, and the
thing a reader wants to look at. `w` is the workspace; `V(φ)` is the deliverable.

That third point is also a genuine mistake in the original paper worth writing up. They spent a
sixth of their search space, and one of their four scoring terms, chasing a quantity that could
have been calculated directly. Worse, they allowed it a 10% tolerance — using a normalisation from
1996 — for something now measured to about half a percent. That loose constraint sits on exactly
the axis their headline result about the tensor-to-scalar ratio depends on.

## What we do about statistics, which is where the old work is weakest

The paper's most quotable claim is that the tensor-to-scalar ratio prefers to land between 0.0001
and 0.0004. It presents this as a hint about nature.

It isn't. It's a histogram of whatever the genetic algorithm happened to visit. Genetic algorithms
are optimisers — they're built to converge onto good regions, not to sample them fairly. Run one
and count what it found, and you're measuring the algorithm's habits at least as much as the
physics. Layer on top the log-spaced coefficient bins, the arbitrary field-range box, and the fact
they never de-duplicated their results (they admit a "seven-fold degeneracy" and that the true
count may be ten times smaller than advertised), and the number doesn't support the interpretation.

So we won't use a genetic algorithm to build the collection. Because our reformulation makes almost
every candidate usable, we can sample the space of working potentials *fairly* instead — which
turns "a pile of stuff an optimiser found" into a well-defined object.

There's still a prior, and there always will be: we have to decide how wiggly a curve is allowed to
be. That choice genuinely determines which formulas can come out the other end. Rather than bury it
in a smoothing constant, we make it an explicit dial and **report every result as a function of that
dial.** If the tensor-to-scalar preference survives turning the dial, it's real and you can say so.
If it moves, then it was always a statement about the prior, and demonstrating *that* is itself a
worthwhile result — it just happens to be a partly deflationary one about your supervisor's paper.

## Why the symbolic regression will work this time, when it didn't in your thesis

In the thesis, ESR's best answer scored a description length of about 49,000 and wasn't even viable
until you manually stuck a prefactor on it. That reads like the method failing. I don't think it
was.

Potentials here have values around 10⁻¹¹. If you hand a fitting routine numbers that small without
telling it what counts as a meaningful error, the goodness-of-fit term blows up and swamps
everything. That's almost certainly what those numbers are. (There's a figure in your `ISR_latex`
folder, `esrfit.png`, that shows the fit visibly coming apart — it's sitting in the directory but
never referenced anywhere in `main.tex`, so it never made it into the document.)

The fix is unglamorous: strip out the parts that don't matter before fitting. Slide every curve so
it sits at the same place, scale every curve to the same size, and fit only the *shape*. As a bonus,
symbolic regression then stops wasting its complexity budget inventing constants to describe an
offset and a prefactor that we already know are meaningless.

The second improvement is bigger. Fitting curves one at a time gives you one formula per curve and
no way to tell signal from noise. Instead, we fit a single formula **skeleton** to hundreds of curves
at once — same structure for all of them, different numbers for each — and rank skeletons by how
well they compress the *entire collection*. This asks a much better question: not "what fits this
one potential" but "what family is this whole region of working potentials?"

That is also, without needing any of their proposed machinery, the thing your thesis defined but
never computed — it proposed scoring a whole function class by the "volume of viable potentials"
inside it. This is the rigorous version of that. And it quietly dissolves the elaborate two-layer
scheme the parent paper proposed for handling constants.

Before believing anything it discovers, we feed it potentials we already know — Starobinsky,
α-attractors, natural inflation, hilltop — and check it recovers them. Ground truth first. Neither
the paper nor the thesis ever did this.

## The order of work

1. **Get the machine ready.** Your environment is bare — jupyter and nothing else, no numpy. It's
   also ARM, so check early that both symbolic regression packages install; that's the main thing
   that could bite. Also, before anything else: check your former supervisor hasn't already
   published the follow-up his paper promised. Colliding with that would be the worst outcome here.
2. **Make the scoring code trustworthy, and reproduce the paper.** The cursor-generated
   `inflation_fitness.py` is genuinely good — better than you'd expect, and honest in its comments
   about where it deliberately departs from the paper. It needs a handful of specific bug fixes and
   an actual test suite. Then reproduce the paper's own published numbers. For an independent
   preprint this isn't optional; it's the thing that earns you the right to disagree with them later.
3. **Build the point-cloud engine.** The new core.
4. **Run symbolic regression on the results.** Discovery.
5. **Physics.** The honest version of the tensor-to-scalar statement, plus re-running under current
   data — the measured tilt appears to have moved significantly since Planck 2018, which on its own
   would relocate the entire viable region and is a publishable result. (I'm going from memory on
   that; confirm the numbers before relying on it.)
6. **Write it up.**

Roughly: steps 1–2 are engineering with a known answer, and you should not skip them because
they're the boring part. Step 3 is where the new idea either works or doesn't. Steps 4–5 are where
the paper comes from.

## One thing I'd push back on, briefly

You picked `V(φ)` over the alternative Hubble-flow space, and that's a defensible call — it keeps
you directly comparable to the paper, which matters for a preprint. The `w` trick above recovers
most of what you'd have gained numerically, so I've built the plan around your choice.

But one argument for the other space doesn't go away: Starobinsky and the α-attractors are *simple*
there and *ugly* in `V(φ)`. Since this entire project ends with ranking formulas by simplicity,
it's worth knowing that your notion of "simple" is partly an artifact of the variable you chose.
Cheap way to settle it: run the final ranking in both spaces and see which known models come out
simple in each. That's a sharp test of whether the complexity ranking means anything, and it costs
very little once everything else exists. It's listed under the optional ideas below.

---

# Part 2 — Technical plan

## Context

Three years ago you did a master's project (`ISR_latex/main.tex`) under Andrei Constantin, extending
Abel–Constantin–Harvey–Lukas, *"Cosmic Inflation and Genetic Algorithms"* (arXiv:2208.13804,
source in `ciga_tex/`). The paper searched **coefficients inside hand-chosen functional forms**
(sextic polynomial; quartic+cos; quartic+exp×quartic) with a binary-encoded GA, and reported
~300k "viable" sextic models plus a claimed preference for `0.0001 ≤ r ≤ 0.0004`. Your thesis
carried this one step further with a single hand-derived crossover
`ats = c₁arctan(c₂ sin(c₃φ))` plus one ESR call.

Both stop at the same place: **the functional form is always an input, never an output.**
The paper's own §IV.B admits this and proposes genetic programming as the fix, then defers it
("we shall return to these more ambitious approaches in a future publication"). Your thesis §6
sketches the better idea and defers it too: *"instead of fitting the chosen function to ESR, we
could replace part of the function with a series of 'point clouds'."*

**This project executes that idea in its strong form.** Represent the potential as a free-form
point cloud with no functional form at all, sample the space of *viable* potentials directly,
then use symbolic regression as a **discovery/compression step on the resulting ensemble** to
read off what functional forms viable inflation actually contains.

```
PAPER / THESIS:   pick form → search coefficients → "is it viable?"
THIS:             search free control points → ensemble of viable curves → SR → "what form IS this?"
```

Target: **independent preprint**. That sets the bar — the reproduction of the paper's own numbers
is not optional, it is the credibility anchor.

**Decisions already taken (yours):** free-form over the whole potential; point clouds live in
V(φ) space (not ε(N) Hubble-flow space); keep the paper's Heaviside/log fitness for now and add a
real likelihood later.

---

## The core technical device: search in w(φ) = V′/V

You chose V(φ) space. The obstruction there is real and must be engineered around: η = V″/V needs
a clean **second** derivative of an interpolated curve, so free knot values produce η noise that
swamps n_s, and most random clouds have no usable slow-roll region at all.

The fix that keeps you in V(φ) space: **put the point cloud on the logarithmic slope**

```
w(φ) ≡ V′/V = d(ln V)/dφ
```

Then every slow-roll quantity is at most **one** derivative away from the cloud:

| quantity | in terms of w | derivatives of the cloud |
|---|---|---|
| ε   | w²/2                | 0 |
| η   | w² + w′             | 1 |
| n_s | 1 − w² + 2w′        | 1 |
| r   | 8w²                 | 0 |
| N   | ∫ dφ / w            | 0 (quadrature) |
| V   | V₀ exp(∫ w dφ)      | 0 (quadrature) |

Three consequences that make the whole programme feasible:

1. **Conditioning.** You get ε(N)-grade numerical behaviour while the axis stays φ and the output
   stays a genuine potential V(φ) — exactly what you asked for.
2. **Viability becomes near-geometric.** End of inflation ε = 0.5 ⟺ |w| = 1. The inflating region
   is just `{|w| < 1}` intersected with `|w² + w′| ≤ 0.5`. You can *construct* clouds that inflate
   for 60 e-folds instead of searching for them. Acceptance goes from ~10⁻⁴ (the paper's regime,
   where one 50-generation run of ~5000 models yields O(10) hits) to ~1.
3. **The amplitude is analytic, not searched.** n_s and r are invariant under V → λV (ε and η both
   are), so the overall scale is a pure normalisation: V₀ is fixed exactly by A_s (or by the paper's
   m_* = 0.027), via `A_s = V/(24π²ε) = V/(12π²w²)`.

That last point is a concrete correction to the paper worth stating in the preprint: their `f_m`
term searches a direction that should have been projected out analytically, spends a whole
dimension of the state space on it (64× on the sextic class), and imposes a **10% window on a
quantity measured to ~0.5%** using a 1996 COBE normalisation. That soft cut on a redundant
direction biases the r-histogram their headline claim rests on.

**Search in w-space; do SR in V-space.** w-space is for conditioning and sampling; the deliverable
the SR consumes is V(φ), which is what you asked for and what a reader wants.

*Caveat to handle explicitly:* the paper shifts V so V(φ_m) = 0, which makes ln V singular at the
minimum. Irrelevant on the inflating plateau (w → ∞ near φ_m just means ε is huge, i.e. not
inflating). Design: parametrise w over the inflating window, reconstruct V there, then attach the
descent to a minimum as a separate step. This decouples cleanly and is worth a paragraph, since
V(φ_m) = 0 is an assumption that only enters via N_* = 58 + ln V(φ_e)/6 and the eternal-inflation cut.

---

## Phase 0 — Environment and repo foundation

The box is a DGX Spark: **aarch64**, NVIDIA GB10, 20 cores, 121 GB unified memory.
`/home/zotar/.venv` currently has **only jupyter** — numpy is not even installed. Everything is greenfield.

- Fresh venv. Core: `numpy scipy sympy matplotlib pandas joblib pyarrow pytest`.
- SR tooling — **validate both install on aarch64 early, this is the main environment risk**:
  - **ESR** (Bartlett/Desmond/Ferreira, arXiv:2211.11461) — pure Python, exhaustive enumeration +
    MDL. This is the tool your thesis used; needed for rigorous complexity/description-length claims.
  - **PySR** — needs Julia via `juliaup` (aarch64 builds exist). Needed to go past the complexity
    ceiling where exhaustive enumeration dies.
- Repo: promote `inflation-miner` to a real package (`pyproject.toml`, `src/`, `tests/`).
  It currently has no `requirements.txt`, no tests, no driver — just a fitness function.
- **Literature check, do this first:** confirm whether Abel/Constantin/Harvey/Lukas already
  published the GP follow-up they promised. An independent preprint that duplicates your former
  supervisor's in-progress work is the worst outcome here.
- **Verify current CMB constraints** rather than assuming. The paper uses Planck 2018
  (n_s = 0.9649 ± 0.0042, r < 0.061). ACT DR6 (2025) moved n_s up substantially (I recall
  P-ACT-LB ≈ 0.974 ± 0.003, and r < 0.038 with BK18) — **check these numbers, don't take mine**.
  If that shift holds it is a >2σ move and the entire viable region relocates, which is a
  publishable result on its own.

## Phase 1 — Trustworthy fitness + reproduction baseline

`inflation_fitness.py` (665 lines, cursor-generated, committed Aug 2026) is genuinely good — I read
it in full. Vectorised, frozen dataclasses, every magic number named in `FitnessConfig`, scoped
`np.errstate` instead of blanket suppression, structured failure tags, and an honest docstring
listing its three deliberate deviations from the paper (lines 40–49). It is a faithful and careful
implementation of §II C. **Build on it.**

Defects to fix (all verified against the source):

| Location | Problem |
|---|---|
| [inflation_fitness.py:380](jupyterlab/inflation-miner/inflation_fitness.py#L380) | `if not np.all(np.isfinite(v)): return None` — one overflow anywhere in [−30,30] rejects the whole potential and mislabels it `"no-minimum"`. `exp(1.2φ)` overflows at φ=30, so this **systematically culls the entire poly+exp class** that `polynomial_plus_exp_polynomial` exists to build. Highest-priority fix. |
| [inflation_fitness.py:486](jupyterlab/inflation-miner/inflation_fitness.py#L486) | `np.where(np.isfinite(integrand), integrand, 0.0)` silently zeroes the divergent V/V′ at a maximum while the neighbouring finite point is still enormous. N_total near a maximum is grid-dependent, and f_N compares N_total to N_*. |
| [inflation_fitness.py:390](jupyterlab/inflation-miner/inflation_fitness.py#L390) | `[:16]` magic cap on minima refinement — a narrow deep minimum poorly sampled by the 0.005 grid can rank out and be missed, so the "global" minimum may be the wrong one. |
| [inflation_fitness.py:341](jupyterlab/inflation-miner/inflation_fitness.py#L341) | `_bisect_root` returns the bracket midpoint on failure and the caller treats it as a genuine maximum, unflagged. |
| [inflation_fitness.py:567](jupyterlab/inflation-miner/inflation_fitness.py#L567), [:582](jupyterlab/inflation-miner/inflation_fitness.py#L582) | `except (ArithmeticError, ValueError, TypeError)` — an `AttributeError`/`IndexError` from a user-supplied `Potential` aborts a whole generation. |

Then:

- **Regression suite.** [inflation_fitness.py:602-665](jupyterlab/inflation-miner/inflation_fitness.py#L602-L665)
  already runs the paper's five worked examples against quoted values and prints a comparison
  table — it has the fixtures but asserts nothing. Convert to `pytest` with tolerances; fill or
  skip the `float("nan")` placeholders at lines 626 and 647.
- **Check the paper's Example 1.** Eq. (III.2) is presented as viable but quotes r ≈ 0.063 > r_max = 0.061,
  and `f_r` has no `f_term` saturation, so it should score non-terminal. Either their figure rounds
  differently from their code or the r cut is applied inconsistently. Resolve before building on it.
- **Reproduce the paper.** Implement the binary-encoded GA (N_pop = 100, 50 generations, 1% mutation,
  single-point crossover, rank-linear selection with α = 3, 6 bits/coefficient log-spaced) over the
  sextic class. Confirm O(10) viable per run and the r-histogram peaked at 1×10⁻⁴–4×10⁻⁴
  (`rStats3.pdf`). This is the anchor.
- **Retire `fitness.py`.** It is untracked, thesis-era, and its vacuum search at
  [fitness.py:17](jupyterlab/inflation-miner/fitness.py#L17) is genuinely wrong — it minimises
  *over* [−5,5] rather than finding the global minimum in [−30,30] and *checking* it lands there,
  so for monotonic V it returns the boundary point ±5 and treats a non-stationary point as φ_m.
  Combined with the module-level `warnings.filterwarnings("ignore")` at line 10 hiding every
  `log(0)`, division-by-zero and `quad` convergence failure. Move to `legacy/` with a README note.

## Phase 2 — The free-form point-cloud engine (the new core)

- `w(φ)` cloud over K knots spanning the inflating window, interpolated by a C¹ spline (only one
  derivative is needed now).
- **Prior over w is a Gaussian process with an explicit length scale ℓ.** This is deliberate and is
  the methodological centrepiece: the smoothness prior is exactly what determines which functional
  forms can appear, so it must be a stated, tunable, reported quantity rather than hidden inside a
  regularisation constant.
- Impose viability by construction: `|w| < 1` over a stretch long enough that ∫dφ/w ≥ N_*.
- Solve V₀ analytically from A_s. Reconstruct V = V₀exp(∫w dφ).
- **Sample density-correctly.** Do *not* use a GA to generate the ensemble — a GA is a biased
  sampler and reusing one would reproduce the paper's central flaw. With acceptance near 1 in
  w-space, straightforward constrained sampling from the GP prior gives a well-defined object:
  *uniform on the viable set under a stated prior*. That is the fix, and it costs nothing extra.
- **Cross-validate:** push every reconstructed V back through the Phase-1 paper fitness and confirm
  f = 0. Two independent codepaths agreeing is what makes this defensible.
- **Canonicalise before anything downstream:** φ-shift (φ_m → 0 or φ_* → 0), V-scale
  (V(φ_*) → 1), and φ → −φ reflection. The paper explicitly did not do this and admits a "seven-fold
  degeneracy" with "the number of non-degenerate solutions may be an order of magnitude fewer than
  the raw number" — which quietly undercuts their own headline 300k/1M counts.

## Phase 3 — Symbolic regression as discovery

Two distinct modes, and the second is the actual novelty:

**(a) Single-curve SR.** ESR at complexity ≤7–9 on canonicalised curves; PySR beyond.
This is what your thesis did, and it failed there: the winner `|arctan(a₀φ)|^a₁` had
L(D) = 48945.66 and fitness −1.0745 (non-viable), only becoming viable after you hand-added a
prefactor. The orphaned `esrfit.png` — present in `ISR_latex/` but never referenced in `main.tex` —
shows the fit visibly diverging for φ > 7 with the posterior collapsed onto essentially one
expression. **The diagnosis is a normalisation bug, not a method failure:** fitting V ~ 10⁻¹¹ with
an unstated σ makes MDL meaningless. Fix by fitting the *shape* only (canonicalised, unit-scaled)
with an explicitly stated noise model. Free result: the φ-shift and V-scale symmetries stop
consuming SR complexity budget on offset and prefactor constants.

**(b) Ensemble SR — shared skeleton, per-curve constants.** Fit one expression *skeleton* to many
viable curves simultaneously, each curve getting its own constants, and rank skeletons by **total
description length over the ensemble**. This falls out naturally from ESR because it already
enumerates skeletons exhaustively — you fit each skeleton to every curve and sum. It answers a
different and better question than (a): not "what expression fits this one potential" but "what
functional family *is* the viable region". It is also, without needing their machinery, the
rigorous version of the class-level fitness your thesis proposed in §3.2 ("normalized n-volume of
viable potentials within the parameter space") and never computed, **and** it dissolves the nested
GP+GA speciation scheme the parent paper proposed in §IV.B to handle constant-fitting.

**Validation the thesis completely lacked:** seed the ensemble with known potentials
(Starobinsky, α-attractor, natural, hilltop) and confirm SR recovers them. Ground truth first,
discovery second.

Deliverable: a ranked catalogue of discovered forms with description lengths, plus the
**complexity floor** for viable inflation — which is the parent paper's own stated open question
("what is the *simplest* form of potential that can produce inflation?").

## Phase 4 — Physics results

1. **The honest r prediction.** Show how the r-distribution moves as the GP length scale ℓ varies.
   The paper's `0.0001 ≤ r ≤ 0.0004` is a density in their discretised log-coefficient prior
   sampled by a biased optimiser, not a posterior. Making the prior-dependence explicit and
   measurable is the correct scientific move — partly deflationary toward the paper, but it lets you
   state the part that *is* robust.
2. **The updated-data result.** Re-run under current constraints. If the ACT DR6 n_s shift is as
   large as I recall, the viable landscape relocates and much of the paper's catalogue moves out of
   bounds — a clean standalone result.
3. **Deferred likelihood** (your call to postpone): compressed (n_s, r, ln A_s) Gaussian with the
   real covariance, N_* marginalised over reheating rather than fixed at 58 + ln V(φ_e)/6.
   This is what would eventually upgrade "catalogue" into "posterior over functional forms".

## Phase 5 — Preprint

---

## Other ideas worth researching (not in the critical path)

You asked for a brainstorm. These are the ones I think are actually worth something, ordered by
value-per-risk. None are required for the above to succeed.

1. **SR in ε₁(N) Hubble-flow space.** You chose V(φ) and the w = V′/V device recovers most of the
   numerical benefit, so this is no longer urgent — but the *scientific* argument stands
   independently: Starobinsky and α-attractors are simple curves in ε(N) and ugly in V(φ), so
   "minimum description length" only tracks physics in that space. Running SR in both and comparing
   which known models come out simple is a sharp, cheap test of whether your complexity ranking
   means anything. Nobody has done SR in this space.
2. **Beyond slow-roll: features and PBHs.** Everything here is first-order slow-roll. Free-form
   clouds will happily produce potentials with features that first-order slow-roll cannot see.
   Solving the exact background + Mukhanov–Sasaki equations for survivors and searching directly on
   the P(k) *shape* is where a genuinely major discovery could live — e.g. a symbolic potential
   fitting CMB while producing PBH-scale enhancement. Highest upside, highest cost. Natural second paper.
3. **Observational equivalence classes.** Observables constrain V only through a narrow window, so
   infinitely many potentials are observationally identical. Rather than counting "viable models"
   (a meaningless number, as the paper's own degeneracy admission shows), characterise the
   degeneracy manifold and report the *simplest representative per class*. Turns a counting exercise
   into a model-selection statement.
4. **Experimental design.** With a posterior over functional forms in hand, ask which future
   observable (r from LiteBIRD/CMB-S4, running α_s, non-Gaussianity, features) most reduces entropy
   over form-space. High-impact and hard to attack, but needs Phase 4.3 done first.
5. **The GNN from your uncompiled section.** `ISR_latex/main.tex` lines 491–496 sit *after*
   `\end{document}` — an entire proposed Graph Neural Network section that maps a function-class tree
   to a tensor of viable-region densities over log-coefficient boxes, invisible in the compiled PDF.
   The idea is sound (learn where in coefficient space a given tree is viable, to pre-localise
   search) and it is now much more tractable, because the Phase-2 ensemble is exactly the training
   set it needs. Worth resurrecting once the ensemble exists.
6. **LLM/transformer-guided SR** (LaSR-style priors over expressions) as an alternative to ESR's
   exhaustive enumeration past complexity ~9. Currently fashionable, genuinely useful for breaking
   the enumeration ceiling, but adds a hard-to-defend prior — keep it out of any headline claim.

---

## Things in the existing work that must not be repeated

Recorded here because the preprint will be judged against them:

- **No baseline.** The thesis's state space was 64³ = 262,144 points — brute-force scannable on a
  laptop. It was never scanned, so the "arm" and the "gap" in `3dscatter1/2.png` could be GA
  sampling bias rather than physics. Any structural claim needs a brute-force or density-correct
  control.
- **Selection artefacts reported as findings.** The thesis's `nstar.png` and `mstar.png` span
  exactly the windows the fitness function hard-cut on, and "both classes favor n_s ≈ 0.9685" was
  reported as a result. Never report a distribution over a quantity you cut on without
  deconvolving the cut.
- **A figure contradicting its own text.** In `rstar.png` the legend and caption both say
  grey = large field / orange = small field, but grey is the spike at r ≈ 0 and orange the hump at
  ≈ 0.002 — the exact opposite of §4.4's claim, on which the headline "10 folds larger or smaller
  than Ref [inflatga]" comparison rests. Check every figure against its caption mechanically.
- **No error bars anywhere.** 1768, 3565, 599, 0.002, 0.00005, 0.9685 are all quoted bare.
- **Five orphaned figures.** `2dplot1/2/3.png`, `esr1.png`, `esrfit.png` exist in `ISR_latex/` and are
  never referenced — including the only visual evidence of ESR fit quality, which shows it failing.
- **Unreproducible.** No seeds, no commit hash, no environment, no runtime.

---

## Verification

- `pytest` green on the paper's five worked examples with stated tolerances (Phase 1).
- GA reproduction recovers the paper's qualitative results: O(10) viable models per 50-generation
  run, r-histogram peaked at 1–4×10⁻⁴, bimodal field excursion (`FieldRangeStats2.pdf`), flat n_s.
- **Round-trip test:** every V reconstructed from a w-cloud scores f = 0 under the independent
  Phase-1 paper fitness. Two codepaths, one answer.
- **SR ground-truth test:** injected Starobinsky / α-attractor / natural / hilltop potentials are
  recovered by the Phase-3 pipeline at their known complexity before any novel form is believed.
- **Prior-sensitivity sweep:** all headline distributions reported as a function of the GP length
  scale ℓ, not at a single hidden value.
- Seeds, environment lockfile, and commit hash recorded for every figure in the preprint.
