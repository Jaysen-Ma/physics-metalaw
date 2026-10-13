#!/usr/bin/env python3
"""Build ``data/catalog/feynman_lectures_v{1,2,3}.json`` from the Feynman Lectures HTML.

Input: ``.scratch/feynman_v1/rows_{I,II,III}.json`` (gitignored), produced by
``.scratch/feynman_v1/extract.py VOL`` (raw MathJax TeX from the user's saved
pages in ``data/feynman/``, per-page macros expanded) and
``.scratch/feynman_v1/build_rows.py VOL`` (the
``quantity = closed-form RHS`` selection rule, sympy parse, operator counts,
within-volume dedupe).

The text of the Feynman Lectures is (c) California Institute of Technology. The
user has decided that the extracted expressions may be shown publicly with a
source link to feynmanlectures.caltech.edu for every expression, so records and
payloads carry ``public: true`` (``PUBLIC`` below). ``--counts-only DIR`` still
writes copies without formula text.
Records that are one of the 100 AI-Feynman benchmark expressions
(``feynman.json``) carry ``benchmark_overlap`` so they are not counted twice.

Run: .venv/bin/python scripts/build_feynman_lectures_catalog.py [--counts-only DIR]
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from physics_metalaw import metalaw as M
from physics_metalaw.catalog import CATALOG_DIR, ROOT, FormulaRecord, apply_constants, attach_counts
from symbol_tex import symbol_tex

PUBLIC = True  # decided by the user (2026-10-13): show expressions with per-expression source links
PRIVATE_FIELDS = ("formula", "original_formula", "tree", "latex", "notes", "name",
                  "constant_names", "variable_names", "c1", "hbar1", "c1_hbar1")


def strip_private(payload: dict) -> dict:
    out = {k: v for k, v in payload.items() if k != "records"}
    out["release_formulas"] = False
    out["records"] = [{k: v for k, v in r.items() if k not in PRIVATE_FIELDS} for r in payload["records"]]
    return out

SCRATCH = ROOT / ".scratch" / "feynman_v1"
FIXES = ROOT / "data" / "feynman_lectures_fixes.csv"
# Main chapter -> branch mapping overrides (user decision 2026-10-13): relativistic electrodynamics counts as
# special relativity. The chapters stay flagged ambiguous; the sensitivity analysis reassigns them to EM.
BRANCH_OVERRIDE = {("I", 34): "Special relativity", ("II", 25): "Special relativity",
                   ("II", 26): "Special relativity", ("II", 28): "Special relativity"}


def load_fixes() -> dict:
    """Manual post-extraction fixes from the full recheck against the source LaTeX (tracked, one row per record)."""
    import csv
    with FIXES.open(newline="", encoding="utf-8") as f:
        return {r["id"]: r for r in csv.DictReader(f)}


def apply_fixes(recs: list, rows: list, vol: str, fixes: dict) -> list:
    out = []
    for rec, row in zip(recs, rows):
        fx = fixes.get(rec["id"])
        if fx is None:
            out.append(rec)
            continue
        if fx["action"] == "reject":
            continue
        row = dict(row)
        if fx["action"] == "replace":
            row["py"] = fx["new_formula"]
        elif fx["action"] == "constants":
            row["constant_names"] = sorted(set(row.get("constant_names", [])) | set(fx["add_constants"].split(";")))
            row["variable_names"] = [v for v in row.get("variable_names", []) if v not in row["constant_names"]]
        new = row_to_json(row, vol)
        new["id"] = rec["id"]
        new["fix"] = {k: fx[k] for k in ("action", "error_type", "reason")} | {"previous_formula": rec["formula"]}
        out.append(new)
    return out
VOLS = {"I": 1, "II": 2, "III": 3}
CITATION = {
    "I": "Feynman, R. P., Leighton, R. B., Sands, M. The Feynman Lectures on Physics, Vol. I (New Millennium Edition, Caltech HTML).",
    "II": "Feynman, R. P., Leighton, R. B., Sands, M. The Feynman Lectures on Physics, Vol. II (New Millennium Edition, Caltech HTML).",
    "III": "Feynman, R. P., Leighton, R. B., Sands, M. The Feynman Lectures on Physics, Vol. III (New Millennium Edition, Caltech HTML).",
}


def benchmark_index() -> dict:
    """(vol, chapter) -> list of benchmark dicts with eq numbers and operator signature."""
    d = json.loads((CATALOG_DIR / "feynman.json").read_text())
    out = {}
    for r in d["records"]:
        if not r.get("in_released"):
            continue
        m = re.search(r"/(I{1,3})_(\d+)\.html", r["source_url"])
        if not m:
            continue
        vol, ch = m.group(1), int(m.group(2))
        eqs = {f"{vol}.{a}.{b}" for a, b in re.findall(r"\((\d+)\.(\d+)[a-z]?", r["source_locator"])}
        sig = (signature(r["formula"]), consts(r["formula"]))
        out.setdefault((vol, ch), []).append(dict(id=r["id"], eqs=eqs, sig=sig, branch=r.get("branch")))
    return out


def signature(py: str):
    """Shape of the expression with every symbol renamed: equal only for the same formula up to notation."""
    import sympy
    try:
        return re.sub(r"Symbol\('[^']*'(, [^)]*)?\)", "S", sympy.srepr(sympy.sympify(py)))
    except Exception:  # noqa: BLE001
        pass
    c = M.count_operators(py)
    return (frozenset((k, v) for k, v in c.items() if not k.startswith("name:") and k != M.NULLARY_NUMBER),
            sum(v for k, v in c.items() if k.startswith("name:") or k == M.NULLARY_NUMBER))


CONST_ALIASES = {"kb": "k", "k_B": "k", "epsilon": "epsilon_0", "epsO": "epsilon_0", "pi": "pi", "c": "c", "k": "k",
                 "hbar": "hbar", "h": "h", "G": "G", "epsilon_0": "epsilon_0", "q_e": "q_e", "mu_0": "mu_0"}


def consts(py: str) -> frozenset:
    """Physical constants named in a formula (benchmark and lecture spellings unified)."""
    names = {k[5:] for k in M.count_operators(py) if k.startswith("name:")}
    return frozenset(CONST_ALIASES[n] for n in names if n in CONST_ALIASES)


def row_to_json(row: dict, vol: str) -> dict:
    eq = row["eq"]
    num = eq.split(".", 1)[1] if eq else None
    rid = f"flp{VOLS[vol]}.{num}" if eq else f"flp{VOLS[vol]}.{row['key'].split('.', 1)[1].replace('#', '.u')}"
    where = f"Eq. ({eq})" if eq else (f"unnumbered display after Eq. ({row['after_eq']})" if row.get("after_eq") else "unnumbered display")
    rec = FormulaRecord(
        corpus=f"feynman_lectures_v{VOLS[vol]}", original_id=eq or row["key"], id=rid,
        formula=row["py"], original_formula=row["py"], output=None,
        n_variables=0, n_numbers=0, n_constants=0, complexity=0,
        original_n_variables=0, original_n_numbers=0, original_n_constants=0, original_complexity=0,
        status="extension", family=None, family_role=None, representative=True, in_released=False,
        source_kind="html", source_url=row["url"],
        source_locator=f"Vol. {vol}, ch. {row['chapter']} {row['chapter_title']}, sec. {row['section']}, {where}",
        notes="; ".join(p for p in (f"LaTeX: {row['tex']}", f"RHS used: {row['rhs']}") if p),
        name="", constant_names=row.get("constant_names", []), variable_names=row.get("variable_names", []),
        consideration="in", branch=BRANCH_OVERRIDE.get((vol, row["chapter"]), row["branch"]), public=PUBLIC,
    )
    rec = attach_counts(apply_constants(rec))
    rec.original_n_variables, rec.original_n_numbers = rec.n_variables, rec.n_numbers
    rec.original_n_constants, rec.original_complexity = rec.n_constants, rec.complexity
    d = rec.to_json()
    d["latex"] = M.formula_latex(rec.formula, symbol_tex=symbol_tex)
    d.update(volume=vol, chapter=row["chapter"], chapter_name=row["chapter_title"], section=row["section"],
             equation=eq, after_equation=row.get("after_eq"), feynman_latex=row["tex"], source_file=row["file"],
             ambiguous_chapter=row.get("ambiguous_chapter", False),
             benchmark_overlap=None, benchmark_overlap_rule=None,
             operator_counts={**{k: int(v) for k, v in rec.operator_counts.items()},
                              M.NULLARY_VARIABLE: rec.n_variables, M.NULLARY_NUMBER: rec.n_numbers,
                              M.NULLARY_CONSTANT: rec.n_constants})
    return d


def assign_overlap(recs: list, vol: str, bench: dict) -> list:
    """Each benchmark expression is matched to at most one kept record: its equation number first,
    else the record in the same chapter with the same operator signature (closest equation number).
    Returns benchmark ids of this volume with no kept record (rejected or never displayed)."""
    unmatched = []
    for (v, ch), bs in sorted(bench.items()):
        if v != vol:
            continue
        cand = [r for r in recs if r["chapter"] == ch and not r["benchmark_overlap"]]
        for b in bs:
            hit = next((r for r in cand if r["equation"] and r["equation"] in b["eqs"] and not r["benchmark_overlap"]), None)
            rule = "same equation number"
            if hit is None:
                sigm = [r for r in cand if not r["benchmark_overlap"] and (signature(r["formula"]), consts(r["formula"])) == b["sig"]]
                if sigm:
                    target = min((int(re.match(r"\d+", e.split(".")[2]).group()) for e in b["eqs"]), default=0)
                    def dist(r):
                        anchor = r["equation"] or r["after_equation"] or ""
                        m = re.match(r"\d+", anchor.split(".")[2]) if anchor.count(".") >= 2 else None
                        return abs(int(m.group()) - target) if m else 999
                    hit = min(sigm, key=dist)
                    rule = "same chapter and operator signature"
            if hit is None:
                unmatched.append(b["id"])
            else:
                hit["benchmark_overlap"], hit["benchmark_overlap_rule"] = b["id"], rule
    return unmatched


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--counts-only", type=Path, help="also write counts+refs-only copies (no formula text) into this directory")
    args = ap.parse_args()
    bench = benchmark_index()
    fixes = load_fixes()
    for vol, n in VOLS.items():
        path = SCRATCH / f"rows_{vol}.json"
        if not path.exists():
            continue
        rows = json.loads(path.read_text())["rows"]
        keep = [r for r in rows if r["status"] == "kept"]
        recs = [row_to_json(r, vol) for r in keep]
        seen = Counter()
        for r in recs:  # several kept rows under one equation number (aligned chains): suffix a, b, ...
            seen[r["id"]] += 1
        dup = {k for k, v in seen.items() if v > 1}
        k2 = Counter()
        for r in recs:
            if r["id"] in dup:
                k2[r["id"]] += 1
                r["id"] = f"{r['id']}.r{k2[r['id']]}"
        n_extracted = len(recs)
        recs = apply_fixes(recs, keep, vol, fixes)
        unmatched = assign_overlap(recs, vol, bench)
        bench_eqs = {e for (v, _), bs in bench.items() if v == vol for b in bs for e in b["eqs"]}
        payload = {
            "corpus": f"feynman_lectures_v{n}", "citation": CITATION[vol], "source_kind": "html",
            "source_url": f"https://www.feynmanlectures.caltech.edu/{vol}_toc.html",
            "board_url": f"https://www.feynmanlectures.caltech.edu/{vol}_toc.html", "board_label": f"Feynman Lectures Vol. {vol}",
            "public": PUBLIC, "release_formulas": PUBLIC,
            "copyright_note": "(c) California Institute of Technology. Expressions shown with a link to the source equation at feynmanlectures.caltech.edu.",
            "extraction": "user-saved chapter HTML (data/feynman/) -> MathJax TeX + page macros -> display equations "
                          "-> selection rule quantity = closed-form RHS (+ Feynman-specific hazards) -> sympy parse_latex -> metalaw.count_operators; "
                          "see .scratch/feynman_v1/",
            "chapters_present": sorted({r["chapter"] for r in rows}),
            "n_equations_seen": len(rows), "n_extracted": n_extracted, "n_records": len(recs),
            "n_fixes": Counter(fx["action"] for k, fx in fixes.items() if fx["volume"] == vol),
            "fixes_file": "data/feynman_lectures_fixes.csv",
            "n_benchmark_overlap": sum(1 for r in recs if r["benchmark_overlap"]),
            "benchmark_ids_in_volume": sorted({b["id"] for (v, _), bs in bench.items() if v == vol for b in bs}),
            "benchmark_ids_unmatched": unmatched,
            "policies": ["released"],
            "records": recs,
        }
        out = CATALOG_DIR / f"feynman_lectures_v{n}.json"
        out.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
        print(f"wrote {out} ({len(recs)} records, {payload['n_benchmark_overlap']} overlap the benchmark)")
        if args.counts_only:
            args.counts_only.mkdir(parents=True, exist_ok=True)
            pub = strip_private(payload)
            for r in pub["records"]:
                for k in ("feynman_latex",):
                    r.pop(k, None)
            (args.counts_only / out.name).write_text(json.dumps(pub, indent=2) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
