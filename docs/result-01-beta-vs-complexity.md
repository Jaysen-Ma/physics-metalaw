# Result 1 — beta is not a corpus constant; it falls steeply with expression complexity

Date: 2026-08-27. Code: `metalaw.py`, `analyse_corpora.py`. Data: the authors'
own released corpus, `github.com/andrei-const/physics-formulae-patterns`
(cloned to `.scratch/`, gitignored).

This was obtained before generating a single potential, from the released data
of Constantin, Bartlett, Desmond & Ferreira 2024 (arXiv:2408.11065).

## The claim being tested

The paper reports one `beta` per corpus and treats the Encyclopaedia
Inflationaris value as anomalous:

| corpus | published beta | reproduced here |
|---|---|---|
| Feynman Lectures | 0.342 +/- 0.010 | 0.352 +/- 0.011 |
| Wikipedia named equations | 0.35 +/- 0.03 | 0.338 +/- 0.016 |
| Cambridge Handbook | 0.33 +/- 0.04 | (data not released) |
| Encyclopaedia Inflationaris | 0.266 +/- 0.014 | 0.282 +/- 0.011 |

Reproduction is within roughly one sigma throughout, and the reproduced gap
between Inflationaris and Feynman (0.070) matches the published gap (0.076).
The residual offset is most likely a fitting-method difference: least squares
here, versus the notebook's Bayesian model comparison.

## Finding

**`beta` depends strongly on the complexity of the expressions it is measured
over, in every corpus.** Splitting each corpus by complexity and refitting:

| complexity window | Inflationaris | Feynman | Wikipedia |
|---|---|---|---|
| 3-7   | (n=2)                  | 0.545 [0.494, 0.591] | 0.525 [0.471, 0.575] |
| 8-11  | 0.329 [0.311, 0.344]   | 0.413 [0.374, 0.458] | 0.353 [0.331, 0.387] |
| 12-17 | 0.321 [0.296, 0.345]   | 0.303 [0.283, 0.320] | 0.343 [0.306, 0.386] |
| 18-25 | 0.296 [0.278, 0.313]   | (n=7)                | (n=3)                |
| 26+   | 0.259 [0.248, 0.278]   | (n=4)                | (n=4)                |

Brackets are 68% intervals from 3000 bootstrap resamples over formulas.

`beta` falls from about 0.54 at complexity 3-7 to about 0.26 at complexity 26+.
That range is roughly four times the Inflationaris-versus-Feynman gap the paper
builds its discussion on.

**In the complexity windows where the corpora actually overlap, they agree.**
At complexity 12-17 all three intervals overlap. At 8-11 Inflationaris sits
below Feynman but overlaps Wikipedia at the margin. Only Inflationaris
populates complexity 18 and above, and that is where the low values live.

The median complexities explain the aggregate result directly: Inflationaris
17, Feynman 11, Wikipedia 10. Inflationaris does not decay more slowly because
it is cosmology. It decays more slowly because its expressions are longer, and
an aggregate `beta` averages over whatever complexity distribution a corpus
happens to have.

## Why this is not an artifact of the fit

Two checks, both in `analyse_corpora.py`:

* **Rank count.** Short expressions contain fewer distinct operator types, and
  the normalised exponential's denominator runs over the ranks present.
  Refitting every bin using only its top 10, or only its top 13, ranks leaves
  the trend intact (Feynman 0.545 -> 0.374 -> 0.297 on the top 10).
* **Arity mix.** In Inflationaris the nullary share is flat across bins
  (0.436, 0.427, 0.441) while `beta` drops 0.332 -> 0.320 -> 0.260, so the
  decline is not simply the nullary fraction falling as trees grow.

## What is still needed before interpreting it

The within-window comparison across corpora stands on its own. The *shape* of
`beta(k)` does not yet, because tree combinatorics alone constrains the arity
mix at fixed node count, so some decline is expected with no physics and no
human involvement. Quantifying that requires the matched control the proposal
already specifies: random expression trees at the same complexity, using the
paper's own generator (a one-dimensional random walk over +1 binary / 0 unary
/ -1 nullary steps, `RandomExpressions.nb`). Until `beta_null(k)` is measured,
report `beta(k)` as description, not as mechanism.

## Counting convention, recovered and validated

The paper does not state its convention exhaustively. It was recovered from the
released formulas and checked against the per-formula complexity column the
authors publish:

* Complexity is total node count of the expression as written.
* Nullary splits into three ranked types: variables, numbers, constants. The
  variable/constant split is *not* recoverable from formula text (`p` in LFI is
  filed as a variable, `A` in RCHI as a constant), so the authors' own
  per-formula columns are used and only operators are recounted.
* An integer exponent of 2, 3 or 4 collapses to a single unary `Pow2`/`Pow3`/
  `Pow4` node and the literal is not counted.
* Unary minus is its own operator.
* **The amplitude prefactor is counted.** Every Inflationaris entry carries
  `(M**4)*(...)`, contributing `M` + `Pow4` + `Mult` = 3 nodes. This settles an
  open question in PROPOSAL.md Section 7.2 for the primary convention.

Validation: 59 of 70 Inflationaris formulas, 97 of 100 Feynman and 38 of 41
Wikipedia reproduce the published complexity exactly. Aggregated, the recovered
counts reproduce the published rank list of Inflationaris on the nose for the
six largest ranks: 362, 220, 168, 167, 81, 70.

## Data-quality notes on the released corpus

* Six Inflationaris entries have unbalanced brackets and do not parse (CNAI,
  CNBI, TMI, CNCI, HBI, SATI). Minimal repairs are recorded in
  `metalaw.EI_BRACKET_REPAIRS`.
* The residual complexity mismatches are hand-counting slips in the source, of
  size 1 to 3.
* The released `Wiki.csv` holds 41 formulas totalling 550 nodes, and carries its
  own totals row agreeing at 551. Whether this matches the corpus size behind
  the published Wikipedia `beta` has not been confirmed against the paper text.
* No Cambridge Handbook data is released, so the largest corpus cannot be
  checked.

## Consequence for this project

The comparison target changes. Comparing a machine-generated catalogue against
Inflationaris's aggregate `beta = 0.266` would be comparing against a number
that is partly an artifact of Inflationaris's complexity distribution. The
comparison has to be made complexity window by complexity window, which is what
PROPOSAL.md Section 5 already argues for on theoretical grounds and this result
now demands on empirical ones.
