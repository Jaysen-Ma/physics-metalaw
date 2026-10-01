"""Build data/catalog/fr.json from the closed f(R) forms in arXiv:0805.1726.

These are examples the review writes down, not a released catalog.
Parameters are variables, including R, so the curvature is not painted
as the gas constant.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from physics_metalaw.catalog import (
    CATALOG_DIR,
    FormulaRecord,
    _name_occurrences,
    apply_constants,
    attach_counts,
)

ABS = "https://arxiv.org/abs/0805.1726"
HTML = "https://arxiv.org/html/0805.1726"

# slug, name, formula, locator, HTML anchor on the formula, notes.
FORMS = [
    ("power", "R^n",
     "R**n",
     "f(R) ∝ R^n",
     "S4.SS1.SSS1.p7.m2",
     "Sec. IV.1.1, the first demonstrative example. "
     "Cited there to Capozziello, Carloni and Troisi (2003)."),
    ("inverse_power", "R − μ^{2(n+1)}/R^n",
     "R - mu**(2*(n+1))/R**n",
     "f(R) = R − μ^{2(n+1)}/R^n",
     "S4.SS1.SSS1.p8.m1",
     "Sec. IV.1.1, the second demonstrative example. "
     "Cited there to Carroll, Duvvuri, Trodden and Turner (2004). "
     "The case n = 1 is R − μ^4/R."),
    ("quadratic", "R + α R^2",
     "R + alpha*R**2",
     "f(R) = R + α R^2",
     "S5.SS1.SSS2.p17.m2",
     "Sec. V.1.2, short-range scalar field. "
     "Cited there to Starobinsky (1980)."),
    ("inverse_plus_quadratic", "R − μ^4/R + a R^2",
     "R - mu**4/R + a*R**2",
     "f(R) = R − μ^4/R + a R^2",
     "S5.E175.m1",
     "Sec. V.2.2, the quadratic correction added to the prototype "
     "R − μ^4/R when discussing stability of de Sitter space."),
    ("chameleon", "R − (1−n) μ^2 (R/μ^2)^n",
     "R - (1-n)*mu**2*(R/mu**2)**n",
     "f(R) = R − (1−n) μ^2 (R/μ^2)^n",
     "S5.E134.m1",
     "Sec. V.1.2, the displayed chameleon example, eq. (134). "
     "Cited there to Carroll et al. (2004) and Amendola et al. (2007)."),
    ("logarithm", "R + α′ ln(R/μ²) + β R^m",
     "R + alpha*log(R/mu**2) + beta*R**m",
     "F(R) = R + α′ ln(R/μ²) + β R^m",
     "https://arxiv.org/pdf/hep-th/0308176.pdf#page=4",
     ""),
    ("exponential", "exp(λ R)",
     "exp(lambd*R)",
     "f(R) = exp(λ R)",
     "https://arxiv.org/pdf/gr-qc/0510130.pdf#page=6",
     ""),
    ("sine", "R − μ^2 sin(μ^2/(R−λ))",
     "R - mu**2*sin(mu**2/(R-lambd))",
     "f(R) = R − μ^2 sin(μ^2/(R−λ))",
     "S5.SS1.SSS2.p16.m11",
     "Sec. V.1.2, the essential-singularity example. "
     "Cited there to Jin, Liu and Li (2006)."),
    ("two_scales", "R − a(R−Λ_1)^{−m} + b(R−Λ_2)^n",
     "R - a*(R-Lambda1)**(-m) + b*(R-Lambda2)**n",
     "f(R) = R − a(R−Λ_1)^{−m} + b(R−Λ_2)^n",
     "S7.SS1.p2.m5",
     "Sec. VII.1, listed among typical examples. "
     "Cited there to Nojiri and Odintsov (2003)."),
    ("starobinsky_2007", "R + λ R_0 [(1+R^2/R_0^2)^{−n} − 1]",
     "R + lambd*R0*((1 + R**2/R0**2)**(-n) - 1)",
     "f(R) = R + λ R_0 [(1+R^2/R_0^2)^{−n} − 1]",
     "S7.SS1.p2.m7",
     "Sec. VII.1, the model tailored to the known constraints. "
     "Cited there to Starobinsky (2007). R_0 is of order H_0^2."),
]


def one(slug: str, name: str, formula: str, locator: str, anchor: str, notes: str) -> FormulaRecord:
    rid = f"fr.{slug}"
    rec = FormulaRecord(
        corpus="fr",
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
        representative=True,
        in_released=False,
        source_kind="compiled",
        source_url=anchor if anchor.startswith("http") else f"{HTML}#{anchor}",
        source_locator=locator,
        notes="",
        name=name,
        variable_names=sorted(_name_occurrences(formula)),
        constant_names=[],
        operator_counts=None,
        consideration="in",
    )
    rec = apply_constants(rec)
    rec = attach_counts(rec)
    rec.original_n_variables = rec.n_variables
    rec.original_n_numbers = rec.n_numbers
    rec.original_n_constants = rec.n_constants
    rec.original_complexity = rec.complexity
    return rec


def main() -> int:
    records = [one(*row) for row in FORMS]
    payload = {
        "corpus": "fr",
        "citation": (
            "Ten closed forms of f(R) written as examples in "
            "Sotiriou and Faraoni, arXiv:0805.1726. Not a released catalog."
        ),
        "source_kind": "compiled",
        "source_url": ABS,
        "board_url": ABS,
        "board_label": "arXiv:0805.1726",
        "n_records": len(records),
        "policies": ["extension"],
        "records": [r.to_json() for r in records],
    }
    out = CATALOG_DIR / "fr.json"
    out.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n",
                   encoding="utf-8")
    print(f"wrote {out} ({len(records)} forms)")
    for rec in records:
        consts = rec.constant_names
        print(f"  {rec.id:28} k={rec.complexity:3d}  "
              f"var={rec.n_variables:2d} num={rec.n_numbers:2d} "
              f"const={consts or '-'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
