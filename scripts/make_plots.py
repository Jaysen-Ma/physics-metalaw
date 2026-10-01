"""Redraw the figures embedded in data/catalog/index.html.

Run from the repository root:

    .venv/bin/python scripts/make_plots.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from physics_metalaw import metalaw as M
from physics_metalaw.catalog import CATALOG_DIR, record_as_entry, scatter_records
from physics_metalaw.sampling import N_DRAWS, SAMPLE_SEED, SAMPLE_SIZE, inflation_sample_betas

OUT = CATALOG_DIR
GREEN = "#2ca02c"


def _freq(rows):
    entries = [record_as_entry(r) for r in rows]
    counts = list(M.rank_counts(entries).values())
    obs = np.sort(np.asarray(counts, dtype=float))[::-1]
    obs = obs[obs > 0]
    freq = obs / obs.sum()
    rank = np.arange(1, obs.size + 1, dtype=float)
    return rank, freq


def fit_exp(rank, freq):
    logf = np.log(freq)

    def model(r, beta):
        logy = -beta * r
        return logy - np.logaddexp.reduce(logy)

    popt, pcov = curve_fit(model, rank, logf, p0=[0.5])
    beta = float(popt[0])
    se = float(np.sqrt(pcov[0, 0]))
    curve = np.exp(model(rank, beta))
    return beta, se, curve


def fit_zipf(rank, freq):
    logf = np.log(freq)

    def model(r, alpha, b):
        logy = -alpha * np.log(r + b)
        return logy - np.logaddexp.reduce(logy)

    popt, _ = curve_fit(model, rank, logf, p0=[2.0, 2.0])
    return np.exp(model(rank, *popt))


def panel(rows, title, path):
    rank, freq = _freq(rows)
    beta, _, curve = fit_exp(rank, freq)
    fig, ax = plt.subplots(figsize=(6.36, 6.0), dpi=100)
    ax.plot(rank, freq, "o", color="black", markersize=5, label="Observed")
    ax.plot(rank, curve, color=GREEN, label=f"Exponential\n$\\beta={beta:.3f}$")
    ax.set_yscale("log")
    ax.set_xlabel("Rank")
    ax.set_ylabel("Frequency")
    ax.set_title(f"{title}\nn={len(rows)}")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)
    print(path.name, len(rows), f"{beta:.3f}")


def zipf_figure():
    names = [
        ("feynman", "Feynman", "zipf-feynman.png"),
        ("wiki", "Wikipedia", "zipf-wiki.png"),
        ("inflation", "Encyclopaedia Inflationaris", "zipf-inflation.png"),
    ]
    for key, title, filename in names:
        rows = scatter_records(key)
        rank, freq = _freq(rows)
        beta, _, exp_curve = fit_exp(rank, freq)
        zm = fit_zipf(rank, freq)
        fig, ax = plt.subplots(figsize=(6.36, 6.0), dpi=100)
        ax.plot(rank, freq, "o", color="black", markersize=5, label="Observed")
        ax.plot(rank, exp_curve, color=GREEN, label=f"Exponential\n$\\beta={beta:.3f}$")
        ax.plot(rank, zm, color="#1f77b4", label="Zipf–Mandelbrot")
        ax.set_yscale("log")
        ax.set_xlabel("Rank")
        ax.set_ylabel("Frequency")
        ax.set_title(f"{title}\nn={len(rows)}")
        ax.legend(frameon=False)
        fig.tight_layout()
        path = OUT / filename
        fig.savefig(path)
        plt.close(fig)
        print(path.name, len(rows), f"{beta:.3f}")


def sample_histogram():
    entries_rows = scatter_records("inflation")
    entries = [record_as_entry(r) for r in entries_rows]
    betas = inflation_sample_betas(entries)
    fig, ax = plt.subplots(figsize=(6.36, 6.0), dpi=100)
    ax.hist(betas, bins=30, color="black", histtype="stepfilled", alpha=0.85)
    ax.axvline(np.median(betas), color=GREEN, lw=1.5,
               label=f"median {np.median(betas):.3f}")
    ax.axvline(0.320, color="0.45", ls="--", lw=1.2, label="f(R) 0.320")
    ax.set_xlabel(r"$\beta$")
    ax.set_ylabel("Samples")
    ax.set_title(f"EI, samples of {SAMPLE_SIZE}\n{N_DRAWS} draws, seed {SAMPLE_SEED}")
    ax.legend(frameon=False)
    fig.tight_layout()
    path = OUT / "ei-sample10-beta.png"
    fig.savefig(path)
    plt.close(fig)
    print(path.name, f"median {np.median(betas):.3f}")


def comparison():
    labels = {}
    for name in ("feynman.json", "wiki.json"):
        payload = json.loads((OUT / name).read_text())
        for rec in payload["records"]:
            if rec.get("branch"):
                labels[rec["id"]] = rec["branch"]
    fw = scatter_records("feynman") + scatter_records("wiki")
    merge = [
        ("Classical mechanics + special relativity",
         ["Classical mechanics", "Special relativity"]),
        ("Electromagnetism + waves and optics",
         ["Electromagnetism", "Waves and optics"]),
        ("Thermodynamics and kinetic theory",
         ["Thermodynamics and kinetic theory"]),
        ("Quantum mechanics", ["Quantum mechanics"]),
        ("Fluid + continuum mechanics",
         ["Fluid mechanics", "Continuum mechanics"]),
    ]
    items = []
    for key, title in [
        ("feynman", "Feynman"),
        ("wiki", "Wikipedia"),
        ("inflation", "Encyclopaedia Inflationaris"),
        ("leemis", "Leemis"),
        ("fr", "f(R)"),
    ]:
        rows = scatter_records(key)
        beta, se, _ = fit_exp(*_freq(rows))
        items.append((title, len(rows), beta, se, "corpus"))
    branch_ids = set()
    for title, parts in merge:
        rows = [r for r in fw if labels.get(r.id) in parts]
        beta, se, _ = fit_exp(*_freq(rows))
        items.append((title, len(rows), beta, se, "branch"))
        branch_ids.update(parts)
        panel(rows, title, OUT / {
            "Classical mechanics + special relativity": "branch-classical-mechanics.png",
            "Electromagnetism + waves and optics": "branch-electromagnetism.png",
            "Thermodynamics and kinetic theory": "branch-thermodynamics.png",
            "Quantum mechanics": "branch-quantum.png",
            "Fluid + continuum mechanics": "branch-fluids.png",
        }[title])
    rows = [r for r in fw if labels.get(r.id) in branch_ids]
    beta, se, _ = fit_exp(*_freq(rows))
    items.append(("These branches together", len(rows), beta, se, "pool"))

    ylabels = [f"{name}  (n={n})" for name, n, _, _, _ in items][::-1]
    betas = np.array([b for _, _, b, _, _ in items][::-1])
    ses = np.array([s for _, _, _, s, _ in items][::-1])
    kinds = [k for *_, k in items][::-1]
    colors = {"corpus": "black", "branch": GREEN, "pool": "#1f77b4"}
    fig, ax = plt.subplots(figsize=(8.2, 6.4), dpi=110)
    y = np.arange(len(ylabels))
    for yy, b, s, k in zip(y, betas, ses, kinds):
        ax.errorbar(b, yy, xerr=s, fmt="o", color=colors[k], capsize=3, markersize=5)
    ax.set_yticks(y)
    ax.set_yticklabels(ylabels)
    ax.set_xlabel(r"$\beta$")
    ax.set_xlim(0.22, 0.46)
    fig.tight_layout()
    path = OUT / "beta-comparison.png"
    fig.savefig(path)
    plt.close(fig)
    print(path.name)


def main() -> int:
    import json
    from physics_metalaw.catalog import expression_scatter
    (OUT / "expression-beta.json").write_text(json.dumps(expression_scatter()) + "\n")
    print("expression-beta.json")
    zipf_figure()
    panel(scatter_records("fr"), "f(R)", OUT / "fr-exponential.png")
    panel(scatter_records("leemis"), "Leemis", OUT / "leemis-exponential.png")
    sample_histogram()
    comparison()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
