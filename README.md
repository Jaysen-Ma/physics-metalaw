<div align="center">
<h1>The grammar of physics is branch-dependent: evidence from the Feynman Lectures</h1>

<!-- TODO(author): name and link, e.g. Chun Hin Ma? -->

<a href="writeup/main.pdf"><img src="https://img.shields.io/badge/Write--up-PDF-red" alt="Write-up PDF"></a>
<a href="data/catalog/index.html"><img src="https://img.shields.io/badge/Formula-ledger-blue" alt="Formula ledger"></a>
</div>

> **Work in progress.** This README, the [write-up](writeup/main.pdf), and the project page are under active revision; reported numbers and wording may change.

## Overview

Prior work ([arXiv:2408.11065](https://arxiv.org/abs/2408.11065)) counts operators in equations of physics and finds that, ranked by frequency, operator-type probabilities follow an exponential law, $p_r \propto e^{-\beta r}$, with decay rate $\beta$ near 0.34 for the Feynman Lectures and Wikipedia and lower for inflaton potentials in Encyclopaedia Inflationaris.

Our likelihood estimate confirms the gap: $\beta = 0.260 \pm 0.006$ for Encyclopaedia Inflationaris against $0.338 \pm 0.010$ for the Feynman benchmark (likelihood ratio 38.9, permutation $p < 0.00025$). We start from that puzzle and test the explanations suggested for it in turn:

1. **Constant conventions.** Rewriting the Feynman benchmark with $c = \hbar = 1$ moves its $\beta$ from 0.338 to $0.348 \pm 0.011$; setting the amplitude prefactor $M = 1$ in all 70 potentials moves Encyclopaedia Inflationaris from 0.260 to $0.257 \pm 0.006$ ([`results/units_check.json`](results/units_check.json)). The gap widens rather than closes.
2. **Lagrangian form.** Ten $f(R)$ gravity Lagrangians, also candidates for one unproven concept, give $0.346 \pm 0.030$, in line with the textbook corpora; random subsets of ten potentials have median 0.257 and rarely reach that value ([`results/ei_sample10.json`](results/ei_sample10.json)).
3. **Physical regime, recast as branch of physics.** The remaining suggestion, a different energy scale or regime, can be tested within one source:

> **Is $\beta$ a single parameter for physics, or does it vary by branch?**

Using one source, we extend the 100 Feynman expressions of the AI Feynman benchmark with every closed-form `quantity = expression` equation in the Feynman Lectures on Physics, Vol. I–III: 1,016 expressions after deduplication, 968 assigned to seven branches. Each extraction was checked by hand against the source LaTeX; 83 corrections and rejections are documented in [`data/feynman_lectures_fixes.csv`](data/feynman_lectures_fixes.csv).

![β by branch, Feynman corpus](writeup/figures/beta_by_branch.png)

> *Figure: $\beta$ (multinomial maximum likelihood) by branch in the Feynman corpus, with bootstrap 1σ errors (Fisher σ for special relativity, marked †). The dashed line and grey band pool all 1,016 expressions.*

**A single $\beta$ is rejected.** Branch-specific $\beta$ fits substantially better than a shared value: likelihood ratio 102.1 on 6 degrees of freedom, permutation $p < 0.00025$ (4,000 permutations; none as extreme), $\Delta$AIC $= -90$, $\Delta$BIC $= -47$.

**Quantum mechanics is low; thermodynamics and fluid mechanics are high.** Against the other branches combined, and after Bonferroni correction for seven tests, quantum mechanics ($0.291 \pm 0.007$) has lower $\beta$; thermodynamics ($0.408 \pm 0.019$), fluid and continuum mechanics ($0.433 \pm 0.020$), and, by a smaller margin, electromagnetism ($0.352 \pm 0.007$) have higher $\beta$. Special relativity ($0.308 \pm 0.011$, with relativistic electrodynamics coded as relativity) is low before correction but not after (corrected $p = 0.05$). Electromagnetism, waves and optics ($0.343 \pm 0.008$) and classical mechanics ($0.359 \pm 0.012$) are mutually consistent (pairwise permutation $p$ between 0.24 and 0.64).

| Branch (Feynman corpus) | Expressions | β | σ |
|---|---:|---:|---:|
| Special relativity | 89 | 0.308 | 0.011† |
| Quantum mechanics | 221 | 0.291 | 0.007 |
| Electromagnetism | 312 | 0.352 | 0.007 |
| Waves and optics | 126 | 0.343 | 0.008 |
| Classical mechanics | 105 | 0.359 | 0.012 |
| Thermodynamics and kinetic theory | 73 | 0.408 | 0.019 |
| Fluid and continuum mechanics | 42 | 0.433 | 0.020 |
| All 1,016 expressions | 1,016 | 0.334 | 0.003 |

Special relativity includes relativistic electrodynamics (Vol. I ch. 34; Vol. II ch. 25, 26, 28). σ: bootstrap over expressions (2,000 replicates, seed 20261013); † Fisher σ for special relativity, whose bootstrap is unreliable without those chapters (32 expressions over 9 ranks); here bootstrap and Fisher agree (bootstrap 0.012, Fisher 0.011).

**Estimating $\beta$.** Operator instances (variables, constants, numbers, $\sin$, $+$, $\times$, …) are counted per group and operator types are ranked by frequency. The model assigns rank $r$ probability $p_r = e^{-\beta r} / \sum_s e^{-\beta s}$, so each rank down multiplies expected frequency by $e^{-\beta}$ (at $\beta = 0.34$, roughly 71% of the previous rank). Treating the $N$ instances as multinomial draws, $\beta$ maximises $\sum_r n_r \log p_r(\beta)$—equivalently, matches the model mean rank to the observed mean rank. Uncertainty comes from a bootstrap over expressions: the expressions of a group are resampled with replacement 2,000 times, counted, re-ranked and refitted, and σ is the standard deviation of those $\beta$. Groups are compared via likelihood-ratio tests with permutation-calibrated $p$-values.

## Contributions

1. **Direct likelihood estimation.** $\beta$ is the multinomial MLE with bootstrap and Fisher errors; group comparisons use permutation-calibrated likelihood-ratio tests. On the published corpora this confirms the lower value of Encyclopaedia Inflationaris: $0.260 \pm 0.006$ against $0.338 \pm 0.010$ for the Feynman benchmark (likelihood ratio 38.9, permutation $p < 0.00025$), compared with the simulation-based $0.266 \pm 0.014$ and $0.342 \pm 0.010$ of the prior study.
2. **Constant conventions ruled out**: natural units and the energy-density prefactor shift $\beta$ by at most 0.01 and do not close the gap.
3. **Lagrangian form ruled out**: ten $f(R)$ Lagrangians have the textbook value of $\beta$, and subsets of ten potentials rarely reach it.
4. **Branch dependence of $\beta$** within one textbook source, the main result: a catalog of 1,016 deduplicated Feynman expressions and a likelihood-ratio test that rejects one $\beta$ for all branches.
5. **Expression catalog with provenance.** The [formula ledger](data/catalog/index.html) builds on the data release accompanying [arXiv:2408.11065](https://arxiv.org/abs/2408.11065) ([physics-formulae-patterns, commit 2aa199b](https://github.com/andrei-const/physics-formulae-patterns/commit/2aa199b5a1d210dacd122f97c3e8086d3e54f752)): 100 Feynman Lectures expressions, 41 Wikipedia equations named after people, and 70 Encyclopaedia Inflationaris potentials. We add 35 post-cutoff Wikipedia entries (76 total), Leemis univariate distributions, ten $f(R)$ Lagrangians, and 1,012 expressions extracted from Feynman Lectures Vol. I–III, each hand-checked against source. Operator counts are recomputed from parsed expressions (physical constants distinguished from variables) and compared to the release. Each entry includes LaTeX rendering, an expression tree, and a source link.

## Formula ledger

```bash
.venv/bin/python -m physics_metalaw.catalog --serve --port 8766
```

Then open `http://127.0.0.1:8766/`. The *Formulas* page lists each expression with tree, counts and source; *β of each expression* plots per-expression $\beta$.

## Reproduce

```bash
python -m venv .venv && .venv/bin/pip install -r requirements.txt
.venv/bin/python scripts/branch_test.py     # β by branch and corpus, all tests (about 80 s) -> results/
.venv/bin/python scripts/units_check.py     # natural-unit and prefactor checks -> results/units_check.json
.venv/bin/python scripts/make_plots.py      # figures -> writeup/figures/, ledger data -> data/catalog/expression-beta.json
.venv/bin/python -m unittest discover -s tests -t .
cd writeup && latexmk -pdf main.tex          # write-up
```

`scripts/branch_test.py` writes all numbers quoted here and in the write-up to `results/branch_test.json` and `writeup/generated/*.tex`. Random seeds: 20261013 (bootstrap, permutations, power simulation) and 20260930 (Encyclopaedia Inflationaris subsamples).

Catalog inputs live in `data/catalog/`. Rebuilding them requires material outside this repository: released CSVs from [arXiv:2408.11065](https://arxiv.org/abs/2408.11065) (for `python -m physics_metalaw.catalog --seed`) and locally saved Feynman Lectures chapter pages with extraction scripts (for `scripts/build_feynman_lectures_catalog.py`). `scripts/build_leemis_catalog.py` and `scripts/build_fr_catalog.py` rebuild the other corpora.

## Data sources and attribution

- **The Feynman Lectures on Physics**, Vol. I–III, R. P. Feynman, R. B. Leighton and M. Sands, New Millennium Edition, © California Institute of Technology, [feynmanlectures.caltech.edu](https://www.feynmanlectures.caltech.edu/). Each extracted expression links to its equation (e.g. `https://www.feynmanlectures.caltech.edu/I_11.html#Eq:I:11:6`) and records volume, chapter, section and equation number.
- **AI Feynman benchmark**: S.-M. Udrescu and M. Tegmark, *AI Feynman: a physics-inspired method for symbolic regression*, Science Advances 6, eaay2631 (2020), as distributed with [arXiv:2408.11065](https://arxiv.org/abs/2408.11065).
- **Wikipedia**, [List of scientific equations named after people](https://en.wikipedia.org/wiki/List_of_scientific_equations_named_after_people) (CC BY-SA 4.0); each expression links to its article.
- **Encyclopaedia Inflationaris**: J. Martin, C. Ringeval and V. Vennin, Phys. Dark Univ. 5–6, 75 (2014), [arXiv:1303.3787](https://arxiv.org/abs/1303.3787).
- **Leemis**: L. M. Leemis and J. T. McQueston, *Univariate distribution relationships*, The American Statistician 62, 45 (2008); [chart](https://www.math.wm.edu/~leemis/chart/UDR/UDR.html).
- **f(R)**: T. P. Sotiriou and V. Faraoni, *f(R) theories of gravity*, Rev. Mod. Phys. 82, 451 (2010), [arXiv:0805.1726](https://arxiv.org/abs/0805.1726), and papers cited per expression in the ledger.

## Citation

This work builds on

```bibtex
@article{constantin2024metalaw,
  title   = {Statistical Patterns in the Equations of Physics and the Emergence of a Meta-Law of Nature},
  author  = {Constantin, Andrei and Bartlett, Deaglan J. and Desmond, Harry and Ferreira, Pedro G.},
  journal = {arXiv preprint arXiv:2408.11065},
  year    = {2024}
}
```
