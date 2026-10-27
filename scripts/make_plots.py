"""Draw the figures of the README and writeup, and the ledger's expression-beta.json.

Run from the repository root, after scripts/branch_test.py:

    .venv/bin/python scripts/make_plots.py

Every number drawn here is read from results/branch_test.json (written by
scripts/branch_test.py), except the Encyclopaedia Inflationaris subsample histogram,
which is recomputed with seed physics_metalaw.sampling.SAMPLE_SEED.

Outputs (PNG for the README, PDF for LaTeX), in results/legacy_figures/:
  beta_by_branch          β by branch, Feynman corpus (main figure)
  rank_frequency_branches rank-frequency of operator counts per branch, with the fitted exponential
  beta_by_corpus          β by corpus
  beta_sensitivity        β by branch under the chapter-assignment and notation checks
  ei_sample10_beta        β of 2000 random subsets of 10 Encyclopaedia Inflationaris potentials
and data/catalog/expression-beta.json (the ledger's "β of each expression" page).
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from physics_metalaw.catalog import CATALOG_DIR, expression_scatter, record_as_entry, scatter_records  # noqa: E402
from physics_metalaw.sampling import N_DRAWS, SAMPLE_SEED, SAMPLE_SIZE, inflation_sample_betas  # noqa: E402

RESULTS = ROOT / "results" / "branch_test.json"
FIG = ROOT / "results" / "legacy_figures"
LONG = {"SR": "Special relativity", "QM": "Quantum mechanics", "EM": "Electromagnetism",
        "Waves+optics": "Waves and optics", "Classical": "Classical mechanics",
        "Thermo": "Thermodynamics", "Fluid+continuum": "Fluid + continuum"}
COLORS = {"SR": "#d62728", "QM": "#9467bd", "EM": "#1f77b4", "Waves+optics": "#17becf",
          "Classical": "#7f7f7f", "Thermo": "#ff7f0e", "Fluid+continuum": "#2ca02c"}
plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})


def save(fig, name):
    FIG.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG / f"{name}.pdf", bbox_inches="tight")
    fig.savefig(FIG / f"{name}.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"{FIG}/{name}.{{pdf,png}}")


def beta_by_branch(res):
    table = res["main"]["table"]
    pooled = res["pooled_all"]
    fig, ax = plt.subplots(figsize=(6.4, 3.6))
    y = np.arange(len(table))[::-1]
    for yy, d in zip(y, table):
        ax.errorbar(d["beta"], yy, xerr=d["quoted_sigma"], fmt="o", color=COLORS[d["branch"]], capsize=3, ms=6)
        ax.annotate(f"{d['beta']:.3f} ± {d['quoted_sigma']:.3f}", (d["beta"] + d["quoted_sigma"] + 0.004, yy),
                    va="center", fontsize=8)
    ax.axvspan(pooled["beta"] - pooled["boot_sigma"], pooled["beta"] + pooled["boot_sigma"], color="0.85", zorder=0)
    ax.axvline(pooled["beta"], color="0.4", lw=1, ls="--", zorder=0)
    ax.text(pooled["beta"] + 0.003, -0.45, f"all {pooled['n']}: {pooled['beta']:.3f} ± {pooled['boot_sigma']:.3f}",
            ha="left", fontsize=7.5, color="0.3")
    ax.set_yticks(y)
    ax.set_yticklabels([f"{LONG[d['branch']]} (n={d['n']})" for d in table])
    ax.set_xlabel(r"$\beta$ (log least squares, $\pm1\sigma$)")
    right = max(d["beta"] + d["quoted_sigma"] for d in table)
    ax.set_xlim(min(d["beta"] - d["quoted_sigma"] for d in table) - 0.02, right + 0.09)
    ax.set_ylim(-0.7, len(table) - 0.4)
    g = res["main"]["tests"][0]
    ax.set_title(f"Feynman corpus by branch: shared vs separate β, LR = {g['LR']:.1f} (df {g['df']}), "
                 f"ΔAIC = {g['dAIC']:.0f}", fontsize=9)
    fig.text(0.01, -0.07, "Bars: bootstrap σ (2000 resamples of expressions).", fontsize=7, color="0.3")
    save(fig, "beta_by_branch")


def rank_frequency(res):
    table = res["main"]["table"]
    fig, axes = plt.subplots(2, 4, figsize=(9.6, 4.8), sharey=True)
    for ax, d in zip(axes.flat, table):
        rf = res["rank_frequency"][d["branch"]]
        ax.plot(rf["rank"], rf["freq"], "o", color="black", ms=3.5)
        ax.plot(rf["rank"], rf["model"], color=COLORS[d["branch"]], lw=1.5)
        ax.set_yscale("log")
        ax.set_title(f"{LONG[d['branch']]}\nn={d['n']}, β={d['beta']:.3f}", fontsize=8.5)
        ax.set_xlim(0, 22)
    ax = axes.flat[-1]
    for d in table:
        rf = res["rank_frequency"][d["branch"]]
        ax.plot(rf["rank"], rf["model"], color=COLORS[d["branch"]], lw=1.2, label=d["branch"])
    ax.set_yscale("log")
    ax.set_title("Fitted exponentials", fontsize=8.5)
    ax.legend(fontsize=6, frameon=False)
    ax.set_xlim(0, 22)
    for a in axes[1]:
        a.set_xlabel("Rank")
    for a in axes[:, 0]:
        a.set_ylabel("Frequency")
    fig.tight_layout()
    save(fig, "rank_frequency_branches")


def beta_by_corpus(res):
    corp = res["corpora"]
    fig, ax = plt.subplots(figsize=(6.4, 3.0))
    y = np.arange(len(corp))[::-1]
    for yy, d in zip(y, corp):
        c = "#d62728" if d["key"] == "inflation" else "black"
        ax.errorbar(d["beta"], yy, xerr=d["boot_sigma"], fmt="o", color=c, capsize=3, ms=6)
        ax.annotate(f"{d['beta']:.3f} ± {d['boot_sigma']:.3f}", (d["beta"] + d["boot_sigma"] + 0.004, yy),
                    va="center", fontsize=8)
    ax.set_yticks(y)
    ax.set_yticklabels([f"{d['label']} (n={d['n']})" for d in corp])
    ax.set_xlabel(r"$\beta$ (log least squares, bootstrap $\pm1\sigma$)")
    right = max(d["beta"] + d["boot_sigma"] for d in corp)
    ax.set_xlim(min(d["beta"] - d["boot_sigma"] for d in corp) - 0.02, right + 0.08)
    save(fig, "beta_by_corpus")


def sensitivity(res):
    sens = res["sensitivity"]
    order = [d["branch"] for d in res["main"]["table"]]
    fig, ax = plt.subplots(figsize=(6.4, 3.8))
    off = np.linspace(-0.3, 0.3, len(sens))
    marks = ["o", "s", "^", "D", "v"]
    for j, s in enumerate(sens):
        bm = {d["branch"]: d for d in s["table"]}
        for i, k in enumerate(order[::-1]):
            if k in bm:
                d = bm[k]
                ax.errorbar(d["beta"], i + off[j], xerr=d["quoted_sigma"], fmt=marks[j], ms=4, capsize=2,
                            color=f"C{j}", label=s["label"] if i == 0 else None)
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels([LONG[k] for k in order[::-1]])
    ax.set_xlabel(r"$\beta$")
    ax.legend(fontsize=7, frameon=False, loc="upper center", bbox_to_anchor=(0.45, -0.16), ncol=3)
    save(fig, "beta_sensitivity")


def ei_sample():
    entries = [record_as_entry(r) for r in scatter_records("inflation")]
    betas = inflation_sample_betas(entries)
    fig, ax = plt.subplots(figsize=(4.8, 3.2))
    ax.hist(betas, bins=30, color="black", alpha=0.85)
    ax.axvline(np.median(betas), color="#d62728", lw=1.5, label=f"median {np.median(betas):.3f}")
    ax.set_xlabel(r"$\beta$ of a random subset")
    ax.set_ylabel("Subsets")
    ax.set_title(f"EI, {N_DRAWS} subsets of {SAMPLE_SIZE} (seed {SAMPLE_SEED})", fontsize=9)
    ax.legend(frameon=False, fontsize=8)
    save(fig, "ei_sample10_beta")
    return dict(median=float(np.median(betas)), p16=float(np.percentile(betas, 16)),
                p84=float(np.percentile(betas, 84)), max=float(betas.max()))


def main() -> int:
    res = json.loads(RESULTS.read_text())
    if res.get("quick"):
        print("warning: results come from a --quick run")
    (CATALOG_DIR / "expression-beta.json").write_text(json.dumps(expression_scatter()) + "\n")
    print("data/catalog/expression-beta.json")
    beta_by_branch(res)
    rank_frequency(res)
    beta_by_corpus(res)
    sensitivity(res)
    ei = ei_sample()
    print("EI subsample:", {k: round(v, 3) for k, v in ei.items()})
    (ROOT / "results" / "ei_sample10.json").write_text(json.dumps(ei, indent=1) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
