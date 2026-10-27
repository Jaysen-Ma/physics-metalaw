"""Build data/catalog/fr.json from closed f(R) forms.

The first ten are examples written in Sotiriou and Faraoni, arXiv:0805.1726.
The rest are later closed forms, each taken from the paper cited on that
record. Parameters are variables, including R, so the curvature is not
painted as the gas constant.
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
    ("hu_sawicki", "Hu–Sawicki",
     "R - c1*R_HS*(R/R_HS)**p/(c2*(R/R_HS)**p + 1)",
     "f(R) = R − c1 R_HS (R/R_HS)^p / [c2 (R/R_HS)^p + 1]",
     "https://arxiv.org/abs/0705.1158",
     "Hu and Sawicki, Phys. Rev. D 76, 064004 (2007), eq. (4) of the arXiv version. "
     "Also written as eq. (2.18) of arXiv:1610.07518."),
    ("tsujikawa", "Tsujikawa",
     "R - mu*R_c*tanh(R/R_c)",
     "f(R) = R − μ R_c tanh(R/R_c)",
     "https://arxiv.org/abs/1002.4928",
     "De Felice and Tsujikawa, Living Rev. Rel. 13, 3 (2010), eq. (4.89). "
     "Also eq. (2.22) of arXiv:1610.07518."),
    ("exponential_gravity", "exponential gravity",
     "R - beta*R_E*(1 - exp(-R/R_E))",
     "f(R) = R − β R_E (1 − exp(−R/R_E))",
     "https://arxiv.org/abs/1610.07518",
     "Eq. (2.24) of arXiv:1610.07518. The same family as "
     "f(R) = R + c1 (1 − exp(−c2 R)) in arXiv:1005.4574."),
    ("appleby_battye", "Appleby–Battye",
     "R/2 + (1/(2*a))*log(cosh(a*R) - tanh(b)*sinh(a*R))",
     "F(R) = R/2 + (1/(2a)) log[cosh(aR) − tanh(b) sinh(aR)]",
     "https://arxiv.org/abs/0705.3199",
     "Appleby and Battye, Phys. Lett. B 654, 7 (2007), eq. (12)."),
    ("miranda_log", "Miranda–Jorás–Waga–Quartin",
     "R - alpha*R_ast*log(1 + R/R_ast)",
     "f(R) = R − α R_* ln(1 + R/R_*)",
     "https://arxiv.org/abs/0905.1941",
     "Miranda, Jorás, Waga and Quartin, Phys. Rev. Lett. 102, 221101 (2009), eq. (3)."),
    ("rational_nojiri", "rational Nojiri–Odintsov",
     "R + (alpha*R**(2*n) - beta*R**n)/(1 + gamma*R**n)",
     "F(R) = R + (α R^{2n} − β R^n)/(1 + γ R^n)",
     "https://arxiv.org/abs/1008.4275",
     "Nojiri and Odintsov, arXiv:1008.4275, eq. (50) "
     "(Prog. Theor. Phys. Suppl. 190, 155 (2011), eq. (6.5)). "
     "The paper writes the fraction as the correction f(R) in F(R) = R + f(R). "
     "The catalog stores F."),
    ("two_tanh", "two hyperbolic steps",
     "R - alpha0*(tanh(b0*(R - R0)/2) + tanh(b0*R0/2)) - alphaI*(tanh(bI*(R - RI)/2) + tanh(bI*RI/2))",
     "F(R) = R − α0 [tanh(b0(R−R0)/2) + tanh(b0 R0/2)] − αI [tanh(bI(R−RI)/2) + tanh(bI RI/2)]",
     "https://arxiv.org/abs/1008.4275",
     "Nojiri and Odintsov, arXiv:1008.4275, the display after eq. (53) "
     "(journal eq. (6.9)). The paper writes the tanh terms as the correction "
     "f(R) in F(R) = R + f(R). The catalog stores F."),
    ("odd_rational", "odd rational",
     "R - ((R - R0)**(2*k + 1) + R0**(2*k + 1))/(f0 + f1*((R - R0)**(2*k + 1) + R0**(2*k + 1)))",
     "F(R) = R − [(R−R0)^{2k+1} + R0^{2k+1}] / {f0 + f1[(R−R0)^{2k+1} + R0^{2k+1}]}",
     "https://arxiv.org/abs/1008.4275",
     "Nojiri and Odintsov, arXiv:1008.4275, eq. (66) "
     "(journal eq. (6.22)). The paper writes the fraction as the correction "
     "f(R) in F(R) = R + f(R). The catalog stores F."),
    ("kruglov_arctan", "Kruglov arctan",
     "R - (b/beta)*atan(beta*R)",
     "F(R) = R − (b/β) arctan(β R)",
     "https://arxiv.org/abs/1310.6915",
     "Kruglov, Phys. Rev. D 89, 064004 (2014)."),
    ("kruglov_arctan_quadratic", "Kruglov arctan of a quadratic",
     "(1/beta)*atan(beta*R - beta**2*R**2)",
     "F(R) = (1/β) arctan(β R − β^2 R^2)",
     "https://doi.org/10.3390/universe1010082",
     "Kruglov, Universe 1, 82 (2015), eq. (1). No separate Einstein–Hilbert term."),
    ("pn_sqrt", "Pérez-Romero y = (Λ/R)^{1/2}",
     "R - 2*Lambda/(1 + b*sqrt(Lambda/R))",
     "f(R) = R − 2Λ / [1 + b (Λ/R)^{1/2}]",
     "https://arxiv.org/abs/1710.05634",
     "Pérez-Romero and Nesseris, Phys. Rev. D 97, 023525 (2018). "
     "y(Λ/R) = (Λ/R)^{1/2} in their eq. (24)."),
    ("pn_inverse", "Pérez-Romero y = Λ/R",
     "R - 2*Lambda/(1 + b*(Lambda/R))",
     "f(R) = R − 2Λ / [1 + b (Λ/R)]",
     "https://arxiv.org/abs/1710.05634",
     "Pérez-Romero and Nesseris, eq. (24) with y = Λ/R. "
     "They identify the n = 1 power with Hu–Sawicki."),
    ("pn_three_halves", "Pérez-Romero y = (Λ/R)^{3/2}",
     "R - 2*Lambda/(1 + b*(Lambda/R)**(3/2))",
     "f(R) = R − 2Λ / [1 + b (Λ/R)^{3/2}]",
     "https://arxiv.org/abs/1710.05634",
     "Pérez-Romero and Nesseris, eq. (24) with y = (Λ/R)^{3/2}."),
    ("pn_arctanh", "Pérez-Romero y = arctanh(Λ/R)",
     "R - 2*Lambda/(1 + b*atanh(Lambda/R))",
     "f(R) = R − 2Λ / [1 + b arctanh(Λ/R)]",
     "https://arxiv.org/abs/1710.05634",
     "Pérez-Romero and Nesseris, eq. (24) with y = arctanh(Λ/R)."),
    ("pn_sin", "Pérez-Romero y = sin(Λ/R)",
     "R - 2*Lambda/(1 + b*sin(Lambda/R))",
     "f(R) = R − 2Λ / [1 + b sin(Λ/R)]",
     "https://arxiv.org/abs/1710.05634",
     "Pérez-Romero and Nesseris, eq. (24) with y = sin(Λ/R)."),
    ("pn_sinh", "Pérez-Romero y = sinh(Λ/R)",
     "R - 2*Lambda/(1 + b*sinh(Lambda/R))",
     "f(R) = R − 2Λ / [1 + b sinh(Λ/R)]",
     "https://arxiv.org/abs/1710.05634",
     "Pérez-Romero and Nesseris, eq. (24) with y = sinh(Λ/R)."),
    ("pn_exp", "Pérez-Romero y = exp(Λ/R)",
     "R - 2*Lambda/(1 + b*exp(Lambda/R))",
     "f(R) = R − 2Λ / [1 + b exp(Λ/R)]",
     "https://arxiv.org/abs/1710.05634",
     "Pérez-Romero and Nesseris, eq. (24) with y = exp(Λ/R)."),
    ("pn_tanh", "Pérez-Romero y = tanh(Λ/R)",
     "R - 2*Lambda/(1 + b*tanh(Lambda/R))",
     "f(R) = R − 2Λ / [1 + b tanh(Λ/R)]",
     "https://arxiv.org/abs/1710.05634",
     "Pérez-Romero and Nesseris, eq. (24) with y = tanh(Λ/R)."),
    ("pn_asin", "Pérez-Romero y = arcsin(Λ/R)",
     "R - 2*Lambda/(1 + b*asin(Lambda/R))",
     "f(R) = R − 2Λ / [1 + b arcsin(Λ/R)]",
     "https://arxiv.org/abs/1710.05634",
     "Pérez-Romero and Nesseris, eq. (24) with y = arcsin(Λ/R)."),
    ("pn_log", "Pérez-Romero y = ln(Λ/R)",
     "R - 2*Lambda/(1 + b*log(Lambda/R))",
     "f(R) = R − 2Λ / [1 + b ln(Λ/R)]",
     "https://arxiv.org/abs/1710.05634",
     "Pérez-Romero and Nesseris, eq. (24) with y = ln(Λ/R)."),
    ("pn_tan", "Pérez-Romero y = tan(Λ/R)",
     "R - 2*Lambda/(1 + b*tan(Lambda/R))",
     "f(R) = R − 2Λ / [1 + b tan(Λ/R)]",
     "https://arxiv.org/abs/1710.05634",
     "Pérez-Romero and Nesseris, eq. (24) with y = tan(Λ/R)."),
    ("pn_arctan", "Pérez-Romero y = arctan(Λ/R)",
     "R - 2*Lambda/(1 + b*atan(Lambda/R))",
     "f(R) = R − 2Λ / [1 + b arctan(Λ/R)]",
     "https://arxiv.org/abs/1710.05634",
     "Pérez-Romero and Nesseris, eq. (24) with y = arctan(Λ/R)."),
    ("pn_erf", "Pérez-Romero y = erf(Λ/R)",
     "R - 2*Lambda/(1 + b*erf(Lambda/R))",
     "f(R) = R − 2Λ / [1 + b erf(Λ/R)]",
     "https://arxiv.org/abs/1710.05634",
     "Pérez-Romero and Nesseris, eq. (24) with y = erf(Λ/R)."),
    ("pn_asinh", "Pérez-Romero y = arcsinh(Λ/R)",
     "R - 2*Lambda/(1 + b*asinh(Lambda/R))",
     "f(R) = R − 2Λ / [1 + b arcsinh(Λ/R)]",
     "https://arxiv.org/abs/1710.05634",
     "Pérez-Romero and Nesseris, eq. (24) with y = arcsinh(Λ/R)."),
    ("pn_polynomial", "Pérez-Romero polynomial",
     "R - 2*Lambda/(1 + b*(b*(Lambda/R) + c*(Lambda/R)**2))",
     "f(R) = R − 2Λ / {1 + b [b (Λ/R) + c (Λ/R)^2]}",
     "https://arxiv.org/abs/1710.05634",
     "Pérez-Romero and Nesseris, polynomial y(Λ/R) below eq. (24)."),
    ("pn_pade", "Pérez-Romero Padé",
     "R - 2*Lambda/(1 + b*((b*(Lambda/R) + c*(Lambda/R)**2)/(1 + p*(Lambda/R) + n*(Lambda/R)**2)))",
     "f(R) = R − 2Λ / {1 + b [b Λ/R + c (Λ/R)^2] / [1 + p Λ/R + n (Λ/R)^2]}",
     "https://arxiv.org/abs/1710.05634",
     "Pérez-Romero and Nesseris, Padé y(Λ/R) below eq. (24)."),
    ("pn_hypergeometric", "Pérez-Romero hypergeometric",
     "R - 2*Lambda + b*H0**2*(Lambda/(R - 3*Lambda))**((1/12)*(-7 + sqrt(73)))*hyp2f1((1/12)*(-7 + sqrt(73)), 3/2 + (1/12)*(-7 + sqrt(73)), 13/6 + 2*(1/12)*(-7 + sqrt(73)), Lambda/(R - 3*Lambda))",
     "f(R) = R − 2Λ + b H0^2 (Λ/(R−3Λ))^α 2F1(α, 3/2+α, 13/6+2α, Λ/(R−3Λ)), α=(−7+√73)/12",
     "https://arxiv.org/abs/1710.05634",
     "Pérez-Romero and Nesseris, eq. (23). Built so the expansion history matches ΛCDM."),
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
        notes=notes,
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
            "Closed forms of f(R). The first ten are examples in "
            "Sotiriou and Faraoni, arXiv:0805.1726. Later records cite "
            "Hu–Sawicki, Tsujikawa, exponential gravity, Appleby–Battye, "
            "Miranda et al., Nojiri–Odintsov, Kruglov, and "
            "Pérez-Romero and Nesseris. Not a released catalog."
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
