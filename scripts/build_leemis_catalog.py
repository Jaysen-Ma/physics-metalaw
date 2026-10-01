#!/usr/bin/env python3
"""Build ``data/catalog/leemis.json`` from ``leemis_appendix.py``.

Run: .venv/bin/python scripts/build_leemis_catalog.py
"""
from __future__ import annotations

import json
import sys
from dataclasses import asdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from physics_metalaw import metalaw as M
from physics_metalaw.catalog import CATALOG_DIR, FormulaRecord, attach_counts, apply_constants
from physics_metalaw.leemis_appendix import LAWS, LeemisLaw

SOURCE_PDF = (
    "https://www.stat.rice.edu/~dobelman/courses/texts/"
    "leemis.distributions.2008amstat.pdf"
)
CHART_URL = "https://www.math.wm.edu/~leemis/chart/UDR/UDR.html"
CITATION = (
    "Leemis, L. M., & McQueston, J. T. (2008). Univariate distribution "
    "relationships. The American Statistician, 62(1), 39–44."
)

# Random variable plus free parameters are nullary variables. pi and e stay
# constants (apply_constants also forces those two).
RV = "x"
UNIVERSAL = frozenset({"pi", "e"})


def param_names(formula: str) -> list[str]:
    """Symbolic leaves other than the random variable and pi, e."""
    if not formula:
        return []
    names = {
        (key[5:-4] if key.endswith("_sym") else key[5:])
        for key in M.count_operators(formula)
        if key.startswith("name:")
    }
    names -= {RV, *UNIVERSAL}
    return sorted(names)


def validate_parse(formula: str) -> None:
    M.count_operators(formula, variables={RV})


def law_to_record(law: LeemisLaw) -> FormulaRecord:
    rid = f"leemis.{law.slug}"
    formula = law.formula.strip()
    if law.parseable and formula:
        validate_parse(formula)
    elif law.parseable and not formula:
        raise ValueError(f"{law.slug}: parseable but empty formula")

    params = [name for name in param_names(formula) if name not in UNIVERSAL]
    rec = FormulaRecord(
        corpus="leemis",
        original_id=rid,
        id=rid,
        formula=formula,
        original_formula=formula,
        output="f",
        n_variables=0,
        n_numbers=0,
        n_constants=0,
        complexity=0,
        original_n_variables=0,
        original_n_numbers=0,
        original_n_constants=0,
        original_complexity=0,
        status="extension",
        family=None,
        family_role=None,
        representative=law.consideration == "in",
        in_released=law.consideration == "in",
        source_kind="pdf",
        source_url=SOURCE_PDF,
        source_locator=f"Appendix A ({law.kind}); chart {CHART_URL}",
        notes="; ".join(
            p
            for p in (
                law.notes,
                f"support: {law.support}",
                f"aliases: {', '.join(law.aliases)}" if law.aliases else "",
            )
            if p
        ),
        name=law.name,
        variable_names=[RV, *params],
        constant_names=[],
        operator_counts=None,
        consideration=law.consideration,
    )
    if formula:
        rec = apply_constants(rec)
        rec = attach_counts(rec)
        rec.original_n_variables = rec.n_variables
        rec.original_n_numbers = rec.n_numbers
        rec.original_n_constants = rec.n_constants
        rec.original_complexity = rec.complexity
    return rec


def record_to_json(rec: FormulaRecord) -> dict:
    if not rec.formula.strip():
        d = asdict(rec)
        d.pop("operator_counts", None)
        d["tree"] = None
        d["latex"] = ""
        d["csv_disagreements"] = {}
        return d
    return rec.to_json()


def main() -> int:
    errors: list[str] = []
    records: list[FormulaRecord] = []
    for law in LAWS:
        try:
            records.append(law_to_record(law))
        except Exception as exc:  # noqa: BLE001 — collect all failures
            errors.append(f"{law.slug}: {exc}")

    n_in = sum(1 for r in records if r.consideration == "in" and r.formula)
    out = CATALOG_DIR / "leemis.json"
    payload = {
        "corpus": "leemis",
        "citation": CITATION,
        "source_kind": "pdf",
        "source_url": SOURCE_PDF,
        "chart_url": CHART_URL,
        "board_label": "Leemis UDR chart",
        "n_records": len(records),
        "n_parseable": n_in,
        "n_out": sum(1 for r in records if r.consideration == "out"),
        "policies": ["released"],
        "records": [record_to_json(r) for r in records],
    }
    out.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {out} ({len(records)} laws, {n_in} in consideration)")
    if errors:
        print("errors:", file=sys.stderr)
        for line in errors:
            print(f"  {line}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
