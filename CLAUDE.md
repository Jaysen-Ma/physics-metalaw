# inflation-miner

Testing whether the operator-frequency "meta-law of nature" (`beta`) is a
property of physics or an artifact of human convention. Originally scoped as
bottom-up generation of viable inflationary potentials compared against
Encyclopaedia Inflationaris; **the current direction (decided 2026-08-27) is
narrower and cheaper — "stream A"**: `beta` is a function of expression
complexity, not a corpus constant (see `docs/result-01-beta-vs-complexity.md`),
and the immediate work is the random-tree null model that determines how much
of that curve is forced by tree combinatorics alone. No inflation physics is
needed for stream A. The generation programme (PROPOSAL.md Steps 0-5) is
retained as a possible follow-up, not current work.

## Read these first, in order

1. **`.notes/STATE.md`** — session handoff: current status, next actions, open
   decisions. Gitignored and machine-local; **on a fresh clone it does not
   exist** — bootstrap from this file, `docs/result-01-beta-vs-complexity.md`
   and PROPOSAL.md, then create it. Read it at the start of every session and
   verify its claims against ground truth before acting on them. Update it as
   milestones land, not only at session end.
2. **`PROPOSAL.md`** — the authoritative research spec: question, hypothesis,
   methodology, outcomes, caveats.

## Fresh-machine setup

```bash
python3 -m venv .venv && .venv/bin/pip install -r requirements.txt
git clone https://github.com/andrei-const/physics-formulae-patterns \
    .scratch/physics-formulae-patterns
.venv/bin/python analyse_corpora.py   # reproduces Result 1 in ~1 min
```

The corpus clone in `.scratch/` is gitignored on purpose — it is the meta-law
authors' data and code, re-cloned rather than vendored. `analyse_corpora.py`
ending with the matched-window table matching
`docs/result-01-beta-vs-complexity.md` confirms the environment is right.

CPU-only; any desktop is sufficient for stream A.

## Standing warnings

- **"300+ viable potentials at complexity 5, ~80% unary" is a recollection, not
  data.** It is not in the thesis; it is from unpublished 2023 work whose code and
  data no longer exist (verified — see PROPOSAL.md Section 3.1). Never cite the
  number; never treat the hypothesis as already evidenced.
- **Viability is probably rare** and the per-expression cost is an inner
  parameter optimisation, not a single evaluation. PROPOSAL.md Section 3.2.
- **`beta` is not a corpus constant.** It falls steeply with expression
  complexity in all three released corpora, and the corpora agree once
  complexity is matched (`docs/result-01-beta-vs-complexity.md`, 2026-08-27).
  Never compare anything against Encyclopaedia Inflationaris's aggregate
  `beta = 0.266`; compare complexity window by complexity window.
- **`beta` also depends on the counting convention, which the paper never
  states.** Folding the power family into one binary `Pow` with its exponent
  counted (`sqrt(x)` = `Pow(x, 0.5)`) leaves Encyclopaedia Inflationaris at
  0.289 but drops Feynman 0.352 -> 0.324, halving the gap the meta-law paper is
  built on (`docs/result-02-counting-convention.md`, 2026-08-31). Never report a
  `beta` without saying which convention produced it, and run the null model
  under the same convention as the measurement it is compared against.
- **The variable/constant split is a per-formula judgement, not a property of
  the symbol**, and the released CSVs give only counts. It is not recoverable;
  take the authors' three columns as given (`metalaw.load_corpus` does).
- **When fitting `beta`, normalise the exponential over the ranks present.** The
  authors' released notebook does not, and the unnormalised fit gives ~1.07
  instead of ~0.27. `metalaw.fit_beta` gets this right.
- Use `.venv/bin/python` (see Fresh-machine setup). Never assume system
  python3 has the scientific stack.

## Layout

- `inflation_fitness.py` — viability criterion (Abel, Constantin, Harvey & Lukas,
  arXiv:2208.13804). `evaluate()` / `fitness()` on a `Potential`; `from_sympy()`
  for symbolic input; all cuts in `FitnessConfig`.
- `fitness.py` — thesis-era original, untracked and superseded. Reference only.
- `metalaw.py` — operator counting under the meta-law paper's convention
  (reverse-engineered and validated against their published data) and the
  normalised-exponential `beta` fit.
- `analyse_corpora.py` — complexity-stratified `beta` for the three released
  corpora, with bootstrap intervals and artifact checks.
- `conventions.py` — the same corpora recounted under alternative counting
  conventions, so the share of `beta` that is a fact about the counting rule can
  be measured. `metalaw.py` stays untouched: it reproduces the authors'
  published numbers and must not drift. Run it directly for the comparison.
- `export_corpus.py` — parses each formula once and writes one JSON record per
  formula to `runs/corpora.json`: tree included, every node tagged with the rank
  bucket it contributes to, alongside the authors' declared complexity, our
  recount, and a flag wherever the two disagree. Nothing is silently corrected.
- `build_viewer.py` — renders that to a standalone page in `runs/`, every
  formula as the nodes it is charged for. Run `export_corpus.py` first.
- `docs/result-01-beta-vs-complexity.md` — Result 1. Read before designing any
  comparison against Encyclopaedia Inflationaris.
- `docs/result-02-counting-convention.md` — Result 2: the convention's share of
  `beta`, the data-quality audit of all three corpora, and the open
  25-vs-23 rank discrepancy in Inflationaris.
- `docs/point-cloud-sr-plan.md` — prior plan for an adjacent project; source of
  the `w = V'/V` device and the "must not be repeated" methodology list.
- `../ciga_tex/InflationGA.tex` — LaTeX source of the ACHL paper, the origin of
  our fitness. Authoritative; prefer it over the arXiv PDF. **DGX box only**;
  on other machines fall back to arXiv:2208.13804.
- `../ISR_latex/` — the MSc thesis (`main.tex`). Sec 4.5 "Function Relatives" is
  the symbolic-degeneracy taxonomy needed for deduplication. **DGX box only.**
- Scratch output goes in `./runs/`, `./logs/` or `./.scratch/` — never `/tmp`.
