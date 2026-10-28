# Exploring variation in the exponential metalaw of physics

<div align="center">

<a href="writeup/main.pdf"><img src="https://img.shields.io/badge/Write--up-PDF-red" alt="Write-up PDF"></a>
<a href="https://jaysen-ma.github.io/physics-metalaw/data/catalog/index.html"><img src="https://img.shields.io/badge/Formula-ledger-blue" alt="Formula ledger"></a>

</div>

An ongoing personal research project asking whether the mathematical structures used in physics equations vary across domains, and whether those differences could eventually help guide symbolic regression.

> **Work in progress.** The analysis, [write-up](writeup/main.pdf), and [formula ledger](https://jaysen-ma.github.io/physics-metalaw/data/catalog/index.html) are still being updated. Reported numbers and wording may change.

## Motivation

Symbolic regression searches over mathematical expressions that fit observed data. The operators and functional forms worth exploring are often domain-dependent: an expression familiar in one branch of science may be an unusual candidate in another. Quantifying these preferences could provide a statistical prior that directs more of the search towards structures typical of the relevant domain.

Previous work found that, after mathematical operator types are ranked from most to least frequent, their frequencies in several physics collections decline approximately exponentially:

$$p_r(\beta)=\frac{e^{-\beta r}}{\sum_s e^{-\beta s}}.$$

The decay rate $\beta$ describes how concentrated the distribution is among the most frequent operator types. That work reported a lower rate for Encyclopaedia Inflationaris than for three broader physics collections. This project uses that difference as a starting point for asking whether operator statistics vary across physical domains.

**This project has not yet produced or validated a symbolic-regression prior.** Its current contribution is an exploratory analysis of corpus-level and per-expression patterns, together with the limitations that would need to be addressed before those patterns could be used as search guidance.

## What the project does

1. Re-estimates decay rates using one fitting procedure and whole-equation resampling.
2. Reconsiders proposed explanations for the earlier inflationary result by changing representation conventions and constructing a comparison collection with similar scope.
3. Expands the Feynman benchmark with expressions from all three volumes of the Feynman Lectures and compares seven physics branches.
4. Fits each expression separately and examines its decay rate against expression-tree complexity, including a random-tree baseline.

The expanded Feynman collection contains 1,016 expressions after deduplication, of which 968 are assigned to the seven branches. The project also maintains sourced catalogues for the original comparison collections, 35 additional Wikipedia expressions, 61 Leemis probability distributions and 37 proposed $f(R)$ curvature functions.

## Estimation

| Collection | Expressions | Fitted β | Nominal individual 95% interval |
|---|---:|---:|---|
| Feynman benchmark | 100 | 0.342 | [0.319, 0.371] |
| Feynman benchmark + Lectures | 1,016 | 0.354 | [0.333, 0.366] |
| Wikipedia | 76 | 0.315 | [0.284, 0.346] |
| Encyclopaedia Inflationaris | 70 | 0.269 | [0.253, 0.293] |
| Leemis probability distributions | 61 | 0.269 | [0.250, 0.304] |
| f(R) gravity | 37 | 0.270 | [0.264, 0.347] |

For each collection or branch, positive operator counts are ranked and a normalised exponential is fitted by **log least squares**. Every observed rank receives equal weight. Uncertainty is estimated from **10,000 whole-equation bootstrap resamples** using seed `20261013`; every resample is recounted, reranked and refitted. A separate chapter bootstrap examines dependence among expressions from the same source chapter.

## Current findings

- The selected collections and branches have different fitted rates, but the comparisons do not isolate one physical explanation for that variation.
- The lower inflationary estimate is reproduced. The tested representation conventions account for little of its difference from the Feynman benchmark, while a comparison collection with similar scope has a close point estimate.
- Within the expanded Feynman collection, quantum mechanics, electromagnetism, waves and optics, and classical mechanics form a lower-rate group in the main specification; thermodynamics and fluid mechanics have higher estimates.
- This grouping is suggestive rather than definitive. Some estimates depend strongly on branch classification, and the reported intervals are less reliable than their nominal 95% label implies.

| Branch | Expressions | Fitted β | Nominal individual 95% interval |
|---|---:|---:|---|
| Special relativity | 89 | 0.382 | [0.322, 0.432] |
| Quantum mechanics | 221 | 0.327 | [0.286, 0.349] |
| Electromagnetism | 312 | 0.325 | [0.313, 0.353] |
| Waves and optics | 126 | 0.346 | [0.303, 0.367] |
| Classical mechanics | 105 | 0.360 | [0.297, 0.377] |
| Thermodynamics and kinetic theory | 73 | 0.403 | [0.371, 0.452] |
| Fluid and continuum mechanics | 42 | 0.422 | [0.377, 0.483] |

![Feynman branch estimates and five comparison collections, with individual bootstrap intervals](writeup/figures/beta_by_branch.png)

The figure shows seven branches as circles and five comparison collections as diamonds: the Feynman benchmark, expanded Feynman collection, f(R) gravity, Encyclopaedia Inflationaris and Wikipedia. These comparisons overlap in their source expressions; they are not all independent datasets.

## Important limitations

- **Restricted equation forms.** The additional Feynman collection contains explicit algebraic expressions. Differential and integral forms, sums, limits, operator hats and complex conjugates are excluded, so the results do not describe the full mathematical language of each branch.
- **Imperfect exponential model.** Rare operator types produce flat steps in the observed tail, whereas a positive-rate exponential declines continuously. Equal-weight log least squares allows these rare ranks to influence the fit despite their few occurrences.
- **Dependence on corpus construction.** Algebraic representation, branch labels and related model variants affect the counts. Removing or reassigning ambiguous chapters moves the special-relativity estimate from 0.382 to 0.231.
- **Approximate intervals.** The bootstrap captures variation from resampling the selected expressions, not uncertainty from omitted equation formats, systematic selection or transcription.

The interval calibration check makes the last limitation concrete. In artificial collections where the true rate is known, an interval labelled 95% should contain that rate about 190 times out of 200. It did so 170–189 times when equations were resampled and only 149–169 times when dependent equations were grouped by chapter. The intervals are therefore often too short. These checks use invented operator counts rather than real expression trees and do not provide corrected intervals for the observed collections.

The write-up presents descriptive comparisons and does not report a formal test that all branch rates are equal.

## Per-expression patterns and symbolic regression

The [formula ledger](https://jaysen-ma.github.io/physics-metalaw/data/catalog/index.html?page=beta) also fits each expression on its own and plots its decay rate $\hat\beta_g$ against its expression-tree complexity $k_g$. The apparent concentration of expressions in parts of this plane suggests a possible expression-level statistic for symbolic regression:

$$\mathcal{L}(g)=\mathcal{L}_{\mathrm{data}}(g)+\lambda\,\mathcal{R}_{d}(\hat\beta_g,k_g).$$

Here, $\mathcal{R}_{d}$ would penalise combinations of decay rate and complexity that are unusual in domain $d$. Conditioning on complexity avoids requiring a short candidate and a long candidate to have the same preferred rate.

The visible regions may partly be combinatorial artefacts: short expressions allow only a limited set of integer operator counts. The ledger therefore includes random expression trees matched to the Feynman complexity distribution as a baseline. More analysis is needed to determine whether the observed regions contain domain information beyond complexity, operator vocabulary and tree structure.

**The proposed penalty has not been implemented or evaluated.** A proper test would learn it from one set of expressions and compare equation-recovery rates and computational effort on separate symbolic-regression problems, with and without the penalty.

## Formula ledger

The public ledger shows every catalogued expression with its source, counted tree and operator totals. Its *β of each expression* page contains corpus, branch and random-tree complexity plots.

To serve it locally:

```bash
.venv/bin/python -m physics_metalaw.catalog --serve --port 8766
```

Then open `http://127.0.0.1:8766/`.

## Reproduce the analysis

Run commands from the repository root. To prepare a Python environment:

```bash
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt matplotlib
```

Run the fitting pipeline and the separate calibration analysis:

```bash
.venv/bin/python scripts/exponential_bootstrap.py
.venv/bin/python scripts/exponential_coverage_check.py
```

The primary pipeline uses 10,000 equation resamples and joint stratified chapter resampling. It writes `results/exponential_fit/summary.json`, bootstrap draws, generated LaTeX tables and manuscript figures. The calibration defaults are 200 datasets per scenario, 400 bootstrap samples and three workers. Smaller diagnostic runs write to separate subdirectories.

To refresh Figures 1 and 2 from the saved results **without refitting or resampling**:

```bash
.venv/bin/python scripts/exponential_bootstrap.py --corpus-plot-only
.venv/bin/python scripts/exponential_bootstrap.py --branch-plot-only
```

These commands update both PDF and PNG figure files; they do not compile the manuscript. Figure 2 includes the seven branches and all five comparison collections listed above.

The existing test suite can be run separately:

```bash
.venv/bin/python -m unittest discover -s tests -t .
```

## Data sources and attribution

- **The Feynman Lectures on Physics**, Vol. I–III, R. P. Feynman, R. B. Leighton and M. Sands, New Millennium Edition, © California Institute of Technology, [feynmanlectures.caltech.edu](https://www.feynmanlectures.caltech.edu/). Each extracted expression links to its equation (e.g. `https://www.feynmanlectures.caltech.edu/I_11.html#Eq:I:11:6`) and records volume, chapter, section and equation number.
- **AI Feynman benchmark**: S.-M. Udrescu and M. Tegmark, *AI Feynman: a physics-inspired method for symbolic regression*, Science Advances 6, eaay2631 (2020), as distributed with [arXiv:2408.11065](https://arxiv.org/abs/2408.11065).
- **Wikipedia**, [List of scientific equations named after people](https://en.wikipedia.org/wiki/List_of_scientific_equations_named_after_people) (CC BY-SA 4.0); each expression links to its article.
- **Encyclopaedia Inflationaris**: J. Martin, C. Ringeval and V. Vennin, Phys. Dark Univ. 5–6, 75 (2014), [arXiv:1303.3787](https://arxiv.org/abs/1303.3787).
- **Leemis**: L. M. Leemis and J. T. McQueston, *Univariate distribution relationships*, The American Statistician 62, 45 (2008); [chart](https://www.math.wm.edu/~leemis/chart/UDR/UDR.html).
- **f(R)**: T. P. Sotiriou and V. Faraoni, *f(R) theories of gravity*, Rev. Mod. Phys. 82, 451 (2010), [arXiv:0805.1726](https://arxiv.org/abs/0805.1726), and papers cited per expression in the ledger.

