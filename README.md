# Exploring variation in the exponential metalaw of physics

<div align="center">

<a href="writeup/main.pdf"><img src="https://img.shields.io/badge/Write--up-PDF-red" alt="Write-up PDF"></a>
<a href="https://jaysen-ma.github.io/physics-metalaw/data/catalog/index.html"><img src="https://img.shields.io/badge/Formula-ledger-blue" alt="Formula ledger"></a>

</div>

An ongoing personal research project exploring how the ranked frequencies of mathematical operators vary across collections of physics equations. It revisits the lower exponential decay rate reported for Encyclopaedia Inflationaris, examines proposed explanations, and compares physics branches within an expanded Feynman collection.

> **Work in progress.** The analysis, [write-up](writeup/main.pdf), and [formula ledger](https://jaysen-ma.github.io/physics-metalaw/data/catalog/index.html) are still being updated. Reported numbers and wording may change.

## Methods and collections

The reported analysis fits a normalised exponential distribution by **log least squares**, with equal weight for each observed operator rank. All positive-count operator types are included, and each collection or branch is ranked separately.

Uncertainty is estimated using **10,000 whole-equation bootstrap resamples**, with seed `20261013`. Each resample is recounted, reranked and refitted. The writeup reports fitted rates and **individual percentile 95% bootstrap intervals**. A joint chapter bootstrap examines dependence among equations sharing source chapters, including chapters represented in more than one branch.

| Collection | Expressions | Fitted β | Nominal individual 95% interval |
|---|---:|---:|---|
| Feynman benchmark | 100 | 0.342 | [0.319, 0.371] |
| Feynman benchmark + Lectures | 1,016 | 0.354 | [0.333, 0.366] |
| Wikipedia | 76 | 0.315 | [0.284, 0.346] |
| Encyclopaedia Inflationaris | 70 | 0.269 | [0.253, 0.293] |
| Leemis probability distributions | 61 | 0.269 | [0.250, 0.304] |
| f(R) gravity | 37 | 0.270 | [0.264, 0.347] |

The expanded Feynman collection adds 916 expressions to the original benchmark; 968 expressions are assigned to seven physics branches. The benchmark is contained in the expanded collection. The f(R) catalogue includes related model variants, so its 37 expressions do not represent 37 independent model families. Source references and transcription conventions are retained in the catalogue and manuscript appendix.

## Current findings

- **The lower inflationary estimate is reproduced.** Its point estimate differs from the Feynman benchmark by approximately 0.073.
- **The tested constant conventions leave most of this gap intact.** Setting `c = ħ = 1` in textbook expressions and removing the common inflationary amplitude prefactor changes the gap to approximately 0.071.
- **The modified-gravity estimate is close to the inflationary value.** Both collections contain hypothetical models focused on one physical concept. The comparison does not establish which features explain their similar rates.
- **Branch estimates suggest a broad grouping.** Quantum mechanics, electromagnetism, waves and optics, and classical mechanics have lower estimates than thermodynamics and fluid mechanics. Their intervals and the sensitivity checks do not establish a sharp division between branches or an energy-scale explanation.

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

## Sensitivity and uncertainty

Removing or reassigning ambiguously classified chapters moves the special-relativity estimate from 0.382 to 0.231. Fluid and continuum mechanics and quantum mechanics have separated individual intervals under both equation and chapter resampling, although the fluid collection contains only four chapters. Appendix B reports the classification, notation, volume and chapter-resampling checks.

Across 200 synthetic datasets in each of three scenarios, individual interval coverage ranged from 74.5% to 94.5%, below the nominal 95% level. These checks use synthetic count vectors rather than real expression trees. The intervals are approximate and do not account for systematic selection or transcription errors.

The writeup presents descriptive comparisons and does not report a formal test that all branch rates are equal. Ranking removes operator identities, and each branch has its own observed operator set. A fitted decay rate alone cannot identify which operators or recurring expression structures explain the differences.

## Motivation: guiding symbolic regression

A longer-term aim is to quantify mathematical structures typical of different physical domains and use them to guide symbolic regression. Examples motivating this idea include products of variables in thermodynamic and fluid equations, and recurring Lorentz-factor structures in special relativity.

One proposal is to add a domain-dependent structural penalty to a data-fitting objective:

$$
\mathcal{L}(g)=\mathcal{L}_{\mathrm{data}}(g)+\lambda\mathcal{R}_{d}(g).
$$

Here, `g` is a candidate expression and the penalty measures how unusual its structure is for domain `d`. Such a penalty could be explored through a custom objective in a framework such as PySR. The aim would be to recover generating equations with fewer candidate evaluations or less computation.

**This application has not been implemented or evaluated in this project.** It would require statistics of operator identities, combinations and recurring subexpressions, beyond the ranked decay rate. A single candidate equation should not be required to reproduce frequencies pooled across a whole collection. Evaluation would compare recovery rates and computational effort on equations separate from those used to learn the penalty.

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

