#!/usr/bin/env python3
"""Exponential fitting with nonparametric bootstrap inference.

Primary fit: unweighted least squares of log frequency against a normalised
exponential, retaining every positive-count operator rank.
Uncertainty: 10,000 equation-level bootstrap replicates, seed 20261013.
Branch differences: 21 pairwise contrasts with centred simultaneous nominal 95% intervals.

The previous log-least-squares permutation results are left in
results/previous_unweighted_logls/ and in results/branch_test.json.

Run from the repository root:

    .venv/bin/python scripts/exponential_bootstrap.py

Redraw only the corpus figure from saved results:

    .venv/bin/python scripts/exponential_bootstrap.py --corpus-plot-only

Redraw the branch figure, including comparison collections, from saved results:

    .venv/bin/python scripts/exponential_bootstrap.py --branch-plot-only
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import shutil
import sys
import time
import subprocess
import platform
import scipy
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from physics_metalaw import metalaw as M  # noqa: E402
from physics_metalaw.catalog import record_as_entry, scatter_records  # noqa: E402
from physics_metalaw.exponential_fit import (  # noqa: E402
    BOOTSTRAP_SEED,
    PRIMARY_THRESHOLD,
    QUANTILE_METHOD,
    bootstrap_type_matrix,
    fit_from_type_counts,
    simultaneous_intervals,
    summarise_bootstrap,
    stratified_cluster_weights,
)

SPEC = importlib.util.spec_from_file_location("branch_test", ROOT / "scripts" / "branch_test.py")
BT = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(BT)
USPEC = importlib.util.spec_from_file_location("units_check", ROOT / "scripts" / "units_check.py")
UNITS = importlib.util.module_from_spec(USPEC)
USPEC.loader.exec_module(UNITS)

OUT = ROOT / "results" / "exponential_fit"
PREV = ROOT / "results" / "previous_unweighted_logls"
FIG = ROOT / "writeup" / "figures"
GEN = ROOT / "writeup" / "generated"
ORDER = list(BT.ORDER)
LONG = dict(BT.LONG)
COLORS = {
    "SR": "#d62728", "QM": "#9467bd", "EM": "#1f77b4", "Waves+optics": "#17becf",
    "Classical": "#7f7f7f", "Thermo": "#ff7f0e", "Fluid+continuum": "#2ca02c",
}
THRESHOLDS = [("none", None)]
CORPORA = (
    ("feynman", "Feynman benchmark"),
    ("feynman_all", "Feynman benchmark + Lectures"),
    ("wiki", "Wikipedia"),
    ("inflation", "Encyclopaedia Inflationaris"),
    ("leemis", "Leemis"),
    ("fr", "f(R)"),
)
plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    digest.update(path.read_bytes())
    return digest.hexdigest()


def archive_previous() -> None:
    if PREV.exists():
        return
    PREV.mkdir(parents=True)
    for name in ("branch_test.json", "branch_beta.csv", "branch_tests.csv",
                 "units_check.json", "ei_sample10.json"):
        src = ROOT / "results" / name
        if src.exists():
            shutil.copy2(src, PREV / name)
    (PREV / "generated").mkdir()
    if GEN.exists():
        for path in GEN.glob("*.tex"):
            shutil.copy2(path, PREV / "generated" / path.name)
    (PREV / "figures").mkdir()
    if FIG.exists():
        for path in FIG.glob("*"):
            if path.is_file():
                shutil.copy2(path, PREV / "figures" / path.name)


def matrix_from_entries(entries) -> np.ndarray:
    per = [M.rank_counts([entry]) for entry in entries]
    keys = sorted({key for counts in per for key in counts})
    if not keys:
        return np.zeros((len(entries), 0))
    return np.array([[counts.get(key, 0) for key in keys] for counts in per], dtype=float)


def chapter_key(record, meta) -> tuple[str, int] | None:
    if meta.get("volume") and meta.get("chapter"):
        return (str(meta["volume"]), int(meta["chapter"]))
    return BT._vc(record.source_url)


def clusters_for(rows, indices: np.ndarray) -> tuple[list[np.ndarray], dict]:
    """Cluster equations by source chapter.

    Equations with no volume/chapter on the catalog record and no chapter in
    the source URL are placed together in one cluster. They are not treated
    as separate independent chapters.
    """
    grouped: dict[tuple[str, int], list[int]] = {}
    missing: list[int] = []
    for local, global_index in enumerate(indices):
        record, _branch, meta = rows[int(global_index)]
        key = chapter_key(record, meta)
        if key is None:
            missing.append(local)
        else:
            grouped.setdefault(key, []).append(local)
    clusters = [np.asarray(members, dtype=int) for members in grouped.values()]
    if missing:
        clusters.append(np.asarray(missing, dtype=int))
    info = {
        "n_chapters": len(grouped),
        "n_missing_chapter": len(missing),
        "missing_policy": (
            "Equations without a chapter identifier share one cluster. "
            "They are not treated as independent chapters."
            if missing else
            "Every equation in this branch has a chapter identifier."
        ),
        "chapter_sizes": sorted((len(members) for members in grouped.values()), reverse=True),
    }
    return clusters, info


def fit_group(matrix, threshold, n_boot, rng, clusters=None):
    fit, betas, statuses = bootstrap_type_matrix(
        matrix, n_boot=n_boot, threshold=threshold, clusters=clusters, rng=rng)
    summary = summarise_bootstrap(fit, betas, statuses)
    summary["n_equations"] = int(matrix.shape[0])
    return summary, betas, np.asarray(statuses, dtype=str)


def difference_block(order, summaries: dict[str, dict], draws: dict[str, np.ndarray]) -> dict:
    names = []
    hats = []
    columns = []
    for i, left in enumerate(order):
        for right in order[i + 1:]:
            names.append({"a": left, "b": right, "label": f"{left} - {right}"})
            beta_a = summaries[left]["beta"]
            beta_b = summaries[right]["beta"]
            if beta_a is None or beta_b is None:
                hats.append(np.nan)
                columns.append(np.full(next(iter(draws.values())).shape, np.nan))
            else:
                hats.append(beta_a - beta_b)
                columns.append(draws[left] - draws[right])
    result = simultaneous_intervals(np.asarray(hats, dtype=float), np.column_stack(columns))
    rows = []
    for meta, interval, individual in zip(
            names, result["intervals"], result["individual_percentile_ci95"]):
        rows.append({**meta, **interval, "individual_percentile_ci95": individual})
    n_exclude = sum(1 for row in rows if row["excludes_zero"])
    validity = result["interval_validity"]
    return {
        "critical_value": result["critical_value"],
        "quantile_method": result["quantile_method"],
        "centring": result["centring"],
        "level": result["level"],
        "interval_validity": validity,
        "n_incomplete": result["n_incomplete"],
        "n_complete": result["n_complete"],
        "n_excluding_zero": n_exclude,
        "any_excludes_zero": None if validity == "unresolved" else bool(n_exclude),
        "pairs": rows,
    }


def json_ready(value):
    if isinstance(value, dict):
        return {str(key): json_ready(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_ready(item) for item in value]
    if isinstance(value, (np.floating, float)):
        number = float(value)
        return None if not np.isfinite(number) else number
    if isinstance(value, (np.integer, int)):
        return int(value)
    if isinstance(value, (np.bool_, bool)):
        return bool(value)
    if value is None:
        return None
    return value


def save_fig(fig, name: str) -> None:
    FIG.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG / f"{name}.pdf", bbox_inches="tight")
    fig.savefig(FIG / f"{name}.png", dpi=200, bbox_inches="tight")
    plt.close(fig)


def fmt(value, digits=3):
    if value is None:
        return "---"
    return f"{value:.{digits}f}"


def fmt_ci(interval) -> str:
    if not interval or interval[0] is None:
        return "---"
    return f"[{fmt(interval[0])}, {fmt(interval[1])}]"


def write_tex(payload: dict) -> None:
    GEN.mkdir(parents=True, exist_ok=True)
    primary = "none"
    lines = ["% Generated by scripts/exponential_bootstrap.py. Do not edit."]
    def macro(name, value):
        lines.append(f"\\newcommand{{\\{name}}}{{{value}}}")
    for key, summary in payload["branches"][primary].items():
        tag = {"SR": "SR", "QM": "QM", "EM": "EM", "Waves+optics": "Waves",
               "Classical": "Classical", "Thermo": "Thermo", "Fluid+continuum": "Fluid"}[key]
        macro(f"beta{tag}", fmt(summary["beta"]))
        macro(f"sig{tag}", fmt(summary["boot_se"]))
        macro(f"n{tag}", summary["n_equations"])
        macro(f"lo{tag}", fmt(summary["ci95"][0]))
        macro(f"hi{tag}", fmt(summary["ci95"][1]))
    diff = payload["differences"][primary]
    macro("simCrit", fmt(diff["critical_value"], 2))
    macro("nExclude", diff["n_excluding_zero"])
    macro("bootN", payload["n_boot"])
    macro("bootSeed", payload["seed"])
    for key, summary in payload["corpora"][primary].items():
        tag = {"feynman": "FeynBench", "feynman_all": "FeynAll", "wiki": "Wiki",
               "inflation": "EI", "leemis": "Leemis", "fr": "fR"}[key]
        macro(f"beta{tag}", fmt(summary["beta"]))
        macro(f"sig{tag}", fmt(summary["boot_se"]))
        macro(f"n{tag}", summary["n_equations"])
        macro(f"lo{tag}", fmt(summary["ci95"][0]))
        macro(f"hi{tag}", fmt(summary["ci95"][1]))
    pooled = payload["pooled"]
    macro("betaPooled", fmt(pooled["beta"]))
    macro("sigPooled", fmt(pooled["boot_se"]))
    macro("nPooled", pooled["n_equations"])
    macro("nLabelled", payload["n_labelled"])
    macro("nFeynAllExpr", payload["n_feynman"])
    macro("nLectAdded", payload["n_lectures_added"])
    (GEN / "numbers.tex").write_text("\n".join(lines) + "\n")

    rows = []
    for key in ORDER:
        summary = payload["branches"][primary][key]
        rows.append(
            f"{LONG[key]} & {summary['n_equations']} & {summary['n_instances']} & {summary['n_ranks']} & "
            f"{fmt(summary['beta'])} & "
            f"{fmt(summary['boot_se'])} & {fmt_ci(summary['ci95'])} \\\\")
    (GEN / "branch_table.tex").write_text("\n".join(rows) + "\n")

    rows = []
    for row in diff["pairs"]:
        rows.append(
            f"{LONG[row['a']]} & {LONG[row['b']]} & {fmt(row['delta'])} & "
            f"{fmt_ci(row['individual_percentile_ci95'])} & {fmt_ci(row['simultaneous_ci95'])} & "
            f"{'yes' if row['excludes_zero'] else 'no'} \\\\")
    (GEN / "pairwise_table.tex").write_text("\n".join(rows) + "\n")

    rows = []
    for key, label in CORPORA:
        summary = payload["corpora"][primary][key]
        rows.append(
            f"{label} & {summary['n_equations']} & {summary['n_instances']} & {summary['n_ranks']} & "
            f"{fmt(summary['beta'])} & "
            f"{fmt(summary['boot_se'])} & {fmt_ci(summary['ci95'])} \\\\")
    (GEN / "corpus_table.tex").write_text("\n".join(rows) + "\n")

    rows = []
    for threshold, _value in THRESHOLDS:
        cells = [fmt(payload["branches"][threshold][key]["beta"]) for key in ORDER]
        n_ex = payload["differences"][threshold]["n_excluding_zero"]
        rows.append(f"{threshold} & " + " & ".join(cells) + f" & {n_ex} \\\\")
    (GEN / "sensitivity_table.tex").write_text("\n".join(rows) + "\n")


def interval_marker(ax, x, y, interval, color, marker="o"):
    lo, hi = interval
    if lo is None or hi is None or x is None:
        return
    ax.plot([lo, hi], [y, y], color=color, lw=1.6, solid_capstyle="butt")
    ax.plot([lo, hi], [y, y], "|", color=color, ms=8)
    ax.plot(x, y, marker, color=color, ms=5)


def plot_branches(payload) -> None:
    """Compare seven branches with five reference collections using saved results."""
    primary = payload["branches"]["none"]
    corpora = payload["corpora"]["none"]
    rows = [(LONG[key], primary[key], COLORS[key], "o") for key in ORDER]
    rows.extend((label, corpora[key], "0.25", "D") for key, label in (
        ("feynman", "Feynman benchmark"),
        ("feynman_all", "Feynman benchmark + Lectures"),
        ("fr", r"$f(R)$"),
        ("inflation", "Encyclopaedia Inflationaris"),
        ("wiki", "Wikipedia"),
    ))
    fig, ax = plt.subplots(figsize=(6.6, 6.5))
    y = np.arange(len(rows), dtype=float)[::-1]
    y[len(ORDER):] -= 0.6
    for yy, (_label, summary, color, marker) in zip(y, rows):
        lo, hi = summary["ci95"]
        interval_marker(ax, summary["beta"], yy, (lo, hi), color, marker)
        if lo is not None and hi is not None:
            ax.text(1.02, yy, f"{summary['beta']:.3f}  [{lo:.3f}, {hi:.3f}]",
                    transform=ax.get_yaxis_transform(), va="center", fontsize=8,
                    clip_on=False)
    ax.axhline((y[len(ORDER)-1] + y[len(ORDER)]) / 2, color="0.8", lw=0.8)
    ax.set_yticks(y)
    ax.set_yticklabels([f"{label} (n={summary['n_equations']})"
                        for label, summary, _color, _marker in rows])
    ax.set_xlabel(r"$\beta$  (nominal individual 95% bootstrap interval)")
    ax.set_title("Feynman branches and comparison collections", fontsize=10)
    save_fig(fig, "beta_by_branch")


def plot_all(payload, rank_curves) -> None:
    primary = payload["branches"]["none"]
    plot_branches(payload)

    diff = payload["differences"]["none"]
    fig, ax = plt.subplots(figsize=(6.6, 6.2))
    pairs = diff["pairs"]
    y = np.arange(len(pairs))[::-1]
    for yy, row in zip(y, pairs):
        lo, hi = row["simultaneous_ci95"]
        color = "#d62728" if row["excludes_zero"] else "0.35"
        interval_marker(ax, row["delta"], yy, (lo, hi), color)
    ax.axvline(0, color="0.5", lw=1)
    ax.set_yticks(y)
    ax.set_yticklabels([f"{LONG[row['a']]} $-$ {LONG[row['b']]}" for row in pairs], fontsize=8)
    ax.set_xlabel(r"Difference of $\beta$  (simultaneous nominal 95% interval)")
    save_fig(fig, "beta_differences")

    fig, axes = plt.subplots(2, 4, figsize=(9.6, 5.0), sharey=True)
    for ax, key in zip(axes.flat, ORDER):
        curve = rank_curves[key]
        rank = np.asarray(curve["rank"])
        freq = np.asarray(curve["frequency"])
        kept = np.asarray(curve["retained"], dtype=bool)
        ax.plot(rank[kept], freq[kept], "o", color=COLORS[key], ms=3.5, label="retained")
        if (~kept).any():
            ax.plot(rank[~kept], freq[~kept], "o", mfc="none", color=COLORS[key], ms=3.5, label="excluded")
        if curve["model"] is not None:
            ax.plot(rank, curve["model"], color="black", lw=1)
        ax.set_yscale("log")
        ax.set_title(f"{LONG[key]}\n" + rf"$\beta$={primary[key]['beta']:.3f}", fontsize=8)
    axes.flat[-1].axis("off")
    axes.flat[0].legend(fontsize=7, frameon=False)
    for ax in axes[1, :3]:
        ax.set_xlabel("Original rank")
    for ax in axes[:, 0]:
        ax.set_ylabel("Frequency")
    fig.tight_layout()
    save_fig(fig, "rank_frequency_branches")

    fig, ax = plt.subplots(figsize=(6.6, 3.8))
    y = np.arange(len(ORDER))[::-1]
    for yy, key in zip(y, ORDER):
        equation = primary[key]
        chapter = payload["chapter_bootstrap"][key]
        interval_marker(ax, equation["beta"], yy + 0.12, equation["ci95"], COLORS[key], "o")
        interval_marker(ax, chapter["beta"], yy - 0.12, chapter["ci95"], COLORS[key], "s")
        if key == ORDER[0]:
            ax.plot([], [], "o", color="0.2", label="equations")
            ax.plot([], [], "s", color="0.2", label="chapters")
    ax.set_yticks(y)
    ax.set_yticklabels([f"{LONG[key]} ({payload['chapter_info'][key]['n_chapters']} ch.)" for key in ORDER])
    ax.set_xlabel(r"$\beta$  (individual nominal 95% interval)")
    ax.legend(fontsize=8, frameon=False)
    save_fig(fig, "chapter_bootstrap")

    plot_corpora(payload)


def plot_corpora(payload) -> None:
    """Render the corpus figure from saved fit and bootstrap summaries."""
    fig, ax = plt.subplots(figsize=(6.6, 3.2))
    corp = payload["corpora"]["none"]
    labels = [label for _key, label in CORPORA]
    keys = [key for key, _label in CORPORA]
    y = np.arange(len(keys))[::-1]
    for yy, key in zip(y, keys):
        summary = corp[key]
        lo, hi = summary["ci95"]
        color = "#d62728" if key == "inflation" else "black"
        interval_marker(ax, summary["beta"], yy, (lo, hi), color)
        if hi is not None:
            ax.annotate(f"{summary['beta']:.3f}", (hi + 0.006, yy), va="center", fontsize=8)
    ax.set_yticks(y)
    ax.set_yticklabels([f"{label} (n={corp[key]['n_equations']})" for label, key in zip(labels, keys)])
    ax.set_xlabel(r"$\beta$  (individual percentile nominal 95% interval)")
    save_fig(fig, "beta_by_corpus")


def plot_ei_subsets(betas, fr_beta):
    fig, ax = plt.subplots(figsize=(5.4, 3.5))
    finite = np.asarray(betas)[np.isfinite(betas)]
    ax.hist(finite, bins=35, color="0.5", edgecolor="white")
    ax.axvline(fr_beta, color="C3", label=r"$f(R)$ fit")
    ax.set_xlabel(r"$\beta$ (all observed ranks)")
    ax.set_ylabel("Number of subsets")
    ax.legend(frameon=False)
    save_fig(fig, "ei_sample10_beta")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--boot", type=int, default=10000)
    parser.add_argument("--corpus-plot-only", action="store_true",
                        help="Redraw the corpus figure from saved summary.json without refitting")
    parser.add_argument("--branch-plot-only", action="store_true",
                        help="Redraw the branch and comparison figure from saved summary.json without refitting")
    args = parser.parse_args()
    if args.boot < 2:
        parser.error("--boot must be at least 2")
    global OUT
    if args.corpus_plot_only or args.branch_plot_only:
        if args.boot != 10000:
            parser.error("Plot-only options use the saved primary results; omit --boot")
        payload = json.loads((OUT / "summary.json").read_text())
        if args.corpus_plot_only:
            plot_corpora(payload)
        if args.branch_plot_only:
            plot_branches(payload)
        print(f"Updated requested figures from {OUT / 'summary.json'}")
        return 0
    if args.boot != 10000:
        OUT = OUT / "diagnostics" / f"boot_{args.boot}"
    t0 = time.time()
    archive_previous()
    OUT.mkdir(parents=True, exist_ok=True)
    rows, dedupe, added = BT.feynman_corpus()
    records = [record for record, _branch, _meta in rows]
    full_matrix = matrix_from_entries([record_as_entry(record) for record in records])
    groups = BT.groups_of(rows)
    labelled = np.concatenate([groups[key] for key in ORDER])
    n_boot = args.boot
    seed_seq = np.random.SeedSequence(BOOTSTRAP_SEED)
    spawn_index = 0

    def generator():
        nonlocal spawn_index
        child = seed_seq.spawn(1)[0]
        spawn_index += 1
        return np.random.default_rng(child), spawn_index

    draws = {}
    statuses = {}

    def store(name, summary, beta_draws, status_draws):
        summary = dict(summary)
        summary["spawn_index"] = summary.pop("_spawn", None)
        draws[name] = np.asarray(beta_draws, dtype=float)
        statuses[name] = np.asarray([{"ok": 0, "boundary": 1, "degenerate": 2, "failed": 3}[s] for s in status_draws], dtype=np.uint8)
        return summary

    def one(name, matrix, threshold, clusters=None):
        rng, index = generator()
        print(f"fit {name} threshold={threshold} n={matrix.shape[0]}", flush=True)
        summary, beta_draws, status_draws = fit_group(matrix, threshold, n_boot, rng, clusters)
        summary["_spawn"] = index
        return store(name, summary, beta_draws, status_draws)

    corpora = {label: {} for label, _value in THRESHOLDS}
    corpus_matrices = {}
    for key, _label in CORPORA:
        if key == "feynman_all":
            corpus_matrices[key] = full_matrix
        else:
            corpus_matrices[key] = matrix_from_entries(
                [record_as_entry(record) for record in scatter_records(key)])
    for label, threshold in THRESHOLDS:
        for key, _name in CORPORA:
            corpora[label][key] = one(f"corpus__{key}__{label}", corpus_matrices[key], threshold)

    branches = {label: {} for label, _value in THRESHOLDS}
    branch_draws = {label: {} for label, _value in THRESHOLDS}
    for label, threshold in THRESHOLDS:
        for key in ORDER:
            summary = one(f"branch__{key}__{label}", full_matrix[groups[key]], threshold)
            branches[label][key] = summary
            branch_draws[label][key] = draws[f"branch__{key}__{label}"]
    differences = {
        label: difference_block(ORDER, branches[label], branch_draws[label])
        for label, _value in THRESHOLDS
    }

    corpus_order = ["feynman", "inflation", "fr"]
    corpus_differences = difference_block(corpus_order, corpora["none"], {key: draws[f"corpus__{key}__none"] for key in corpus_order})
    pooled = one("pooled_all", full_matrix, PRIMARY_THRESHOLD)
    # A source chapter can occur in more than one branch. Resample it once
    # and apply its multiplicity to every branch containing its equations.
    chapter_info = {}
    chapter_boot, chapter_draws = {}, {}
    chapter_rows = {}
    for idx in labelled:
        record, branch, meta = rows[int(idx)]
        key = chapter_key(record, meta)
        if key is None:
            raise ValueError("Chapter resampling requires source chapter identifiers")
        chapter_rows.setdefault(key, []).append(int(idx))
    cluster_keys = list(chapter_rows)
    clusters = [np.asarray(chapter_rows[key], dtype=int) for key in cluster_keys]
    memberships = {key: set(groups[key].tolist()) for key in ORDER}
    strata_by_branch = {key: [] for key in ORDER}
    shared = []
    for c, indices in enumerate(clusters):
        counts = {key: sum(int(i) in memberships[key] for i in indices) for key in ORDER}
        dominant = max(ORDER, key=lambda key: counts[key])
        strata_by_branch[dominant].append(c)
        if sum(n > 0 for n in counts.values()) > 1:
            shared.append({"chapter": cluster_keys[c], "branch_counts": counts})
    strata = [np.asarray(ids) for ids in strata_by_branch.values() if ids]
    chapter_rng, chapter_spawn = generator()
    raw = {key: np.full(n_boot, np.nan) for key in ORDER}
    raw_status = {key: np.full(n_boot, "failed", dtype=object) for key in ORDER}
    for b in range(n_boot):
        weights = stratified_cluster_weights(clusters, strata, len(rows), chapter_rng)
        for key in ORDER:
            fit = fit_from_type_counts(weights[groups[key]] @ full_matrix[groups[key]], PRIMARY_THRESHOLD)
            raw[key][b] = np.nan if fit.beta is None else fit.beta
            raw_status[key][b] = fit.status
    for key in ORDER:
        _, info = clusters_for(rows, groups[key])
        info["n_stratum_chapters"] = len(strata_by_branch[key])
        chapter_info[key] = info
        original = fit_from_type_counts(full_matrix[groups[key]].sum(axis=0), PRIMARY_THRESHOLD)
        summary = summarise_bootstrap(original, raw[key], raw_status[key])
        summary.update(n_equations=len(groups[key]), chapters=info, _spawn=chapter_spawn)
        chapter_boot[key] = store(f"chapter__{key}", summary, raw[key], raw_status[key])
        chapter_draws[key] = raw[key]
    chapter_differences = difference_block(ORDER, chapter_boot, chapter_draws)
    chapter_design = {
        "method": "joint stratified source-chapter bootstrap",
        "strata": "Each chapter belongs to the branch containing most of its retained equations; ties use ORDER. Resample the observed number of chapters within each stratum.",
        "shared_chapters": shared,
        "stratum_chapter_counts": {key: len(ids) for key, ids in strata_by_branch.items()},
        "dependence": "All equations in a selected chapter receive the same multiplicity across branches."
    }

    scenarios = [
        ("ambiguous_dropped", None, lambda meta, branch: meta["ambiguous"]),
        ("ambiguous_reassigned",
         lambda meta, branch: (BT.ALT.get((meta["volume"], meta["chapter"]), branch) if meta["ambiguous"] else branch),
         lambda meta, branch: meta["ambiguous"] and (meta["volume"], meta["chapter"]) in BT.ALT and BT.ALT[(meta["volume"], meta["chapter"])] is None),
        ("notation_dropped", None, lambda meta, branch: meta["notation"]),
        ("volume_iii_dropped", None, lambda meta, branch: meta["volume"] == "III"),
    ]
    ablations = {}
    for name, relabel, drop in scenarios:
        grouped = BT.groups_of(rows, relabel, drop)
        summaries = {}
        local_draws = {}
        for key in ORDER:
            if key not in grouped:
                continue
            summary = one(f"ablation__{name}__{key}", full_matrix[grouped[key]], PRIMARY_THRESHOLD)
            summaries[key] = summary
            local_draws[key] = draws[f"ablation__{name}__{key}"]
        present = [key for key in ORDER if key in summaries]
        ablations[name] = {
            "n": {key: int(len(grouped[key])) for key in present},
            "branches": summaries,
            "differences": difference_block(present, summaries, local_draws) if len(present) >= 2 else None,
        }

    units = {}
    for key in ("feynman", "wiki"):
        records_c = scatter_records(key)
        rewritten = [UNITS.C.rewrite_units(record, "c1_hbar1") for record in records_c]
        units[key] = {
            "as_written": corpora["none"][key],
            "c1_hbar1": one(
                f"units__{key}__c1",
                matrix_from_entries([UNITS.entry(record, formula) for record, formula in zip(records_c, rewritten)]),
                PRIMARY_THRESHOLD),
            "n_changed": int(sum(formula != record.formula for record, formula in zip(records_c, rewritten))),
        }
    inflation_records = scatter_records("inflation")
    rewritten = [M.set_names_to_one(record.formula, ("M",)) for record in inflation_records]
    units["inflation"] = {
        "as_written": corpora["none"]["inflation"],
        "M1": one(
            "units__inflation__M1",
            matrix_from_entries([UNITS.entry(record, formula) for record, formula in zip(inflation_records, rewritten)]),
            PRIMARY_THRESHOLD),
        "n_changed": int(sum(formula != record.formula for record, formula in zip(inflation_records, rewritten))),
    }

    inflation_matrix = corpus_matrices["inflation"]
    subset_rng = np.random.default_rng(20260930)
    subset_betas = []
    subset_status = []
    for _draw in range(2000):
        take = subset_rng.choice(inflation_matrix.shape[0], size=10, replace=False)
        fit = fit_from_type_counts(inflation_matrix[take].sum(axis=0), PRIMARY_THRESHOLD)
        subset_betas.append(np.nan if fit.beta is None else fit.beta)
        subset_status.append(fit.status)
    subset_betas = np.asarray(subset_betas, dtype=float)
    draws["ei_sample10"] = subset_betas
    statuses["ei_sample10"] = np.asarray([{"ok": 0, "boundary": 1, "degenerate": 2, "failed": 3}[v] for v in subset_status], dtype=np.uint8)
    finite_subsets = subset_betas[np.isfinite(subset_betas)]
    ei_sample = {
        "seed": 20260930,
        "n_draws": 2000,
        "sample_size": 10,
        "n_finite": int(finite_subsets.size),
        "n_degenerate_or_failed": int(subset_betas.size - finite_subsets.size),
        "median": float(np.median(finite_subsets)) if finite_subsets.size else None,
        "p16": float(np.quantile(finite_subsets, 0.16)) if finite_subsets.size else None,
        "p84": float(np.quantile(finite_subsets, 0.84)) if finite_subsets.size else None,
        "max": float(finite_subsets.max()) if finite_subsets.size else None,
    }

    rank_curves = {}
    for key in ORDER:
        fit = fit_from_type_counts(full_matrix[groups[key]].sum(axis=0), PRIMARY_THRESHOLD)
        counts = fit.counts
        freq = counts / counts.sum()
        retained = np.zeros(counts.size, dtype=bool)
        retained[np.asarray(fit.retained_ranks) - 1] = True
        if fit.beta is None:
            model = None
        else:
            logp = -fit.beta * np.arange(1, counts.size + 1)
            logp -= np.logaddexp.reduce(logp)
            model = np.exp(logp).tolist()
        rank_curves[key] = {
            "rank": list(range(1, counts.size + 1)),
            "frequency": freq.tolist(),
            "retained": retained.tolist(),
            "model": model,
            "beta": fit.beta,
        }

    code_path = ROOT / "physics_metalaw" / "exponential_fit.py"
    data_files = [
        "data/catalog/feynman.json",
        "data/catalog/feynman_lectures_v1.json",
        "data/catalog/feynman_lectures_v2.json",
        "data/catalog/feynman_lectures_v3.json",
        "data/catalog/wiki.json",
        "data/catalog/inflation.json",
        "data/catalog/leemis.json",
        "data/catalog/fr.json",
    ]
    payload = {
        "method": "exponential fitting with nonparametric bootstrap inference",
        "seed": BOOTSTRAP_SEED,
        "seed_sequence": "numpy.random.SeedSequence spawn order recorded as spawn_index",
        "n_boot": n_boot,
        "primary_threshold": PRIMARY_THRESHOLD,
        "thresholds": {label: value for label, value in THRESHOLDS},
        "bounds": [0.0, 5.0],
        "bound_tol": 1e-4,
        "quantile_method": QUANTILE_METHOD,
        "centring": "original_estimate",
        "resampling_unit": "equation",
        "chapter_resampling_unit": "source chapter",
        "overlap": {
            "feynman_benchmark_inside_feynman_all": True,
            "formal_difference": (
                "No difference interval is reported between the Feynman benchmark and the "
                "expanded Feynman corpus, because the benchmark equations are part of the "
                "expanded corpus. Equation resampling assumes independence between disjoint "
                "branches. Joint chapter resampling preserves shared-chapter dependence."
            ),
        },
        "provenance": {
            "fit_module_sha256": sha256(code_path),
            "code_sha256": {name: sha256(ROOT / name) for name in ["physics_metalaw/exponential_fit.py", "physics_metalaw/metalaw.py", "physics_metalaw/catalog.py", "scripts/exponential_bootstrap.py", "scripts/branch_test.py", "scripts/build_feynman_lectures_catalog.py", "scripts/units_check.py"]},
            "git_head": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
            "versions": {"python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__},
            "data_sha256": {name: sha256(ROOT / name) for name in data_files},
            "dedupe": dedupe,
            "lectures_added": added,
        },
        "n_feynman": len(rows),
        "n_lectures_added": int(sum(1 for _r, _b, meta in rows if meta["source"] != "benchmark")),
        "n_labelled": int(len(labelled)),
        "corpora": corpora,
        "corpus_differences": corpus_differences,
        "branches": branches,
        "differences": differences,
        "pooled": pooled,
        "chapter_info": chapter_info,
        "chapter_design": chapter_design,
        "chapter_bootstrap": chapter_boot,
        "chapter_differences": chapter_differences,
        "ablations": ablations,
        "units": units,
        "ei_sample10": ei_sample,
        "rank_curves": rank_curves,
        "runtime_s": time.time() - t0,
    }
    safe = json_ready(payload)
    (OUT / "summary.json").write_text(json.dumps(safe, indent=1) + "\n")
    np.savez_compressed(OUT / "draws.npz", **{key: value for key, value in draws.items()},
                        **{f"status__{key}": value for key, value in statuses.items()})
    if n_boot == 10000:
        write_tex(safe)
        plot_all(safe, rank_curves)
        plot_ei_subsets(subset_betas, corpora["none"]["fr"]["beta"])
    print(f"done in {time.time() - t0:.0f}s  primary exclusions "
          f"{safe['differences']['none']['n_excluding_zero']}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
