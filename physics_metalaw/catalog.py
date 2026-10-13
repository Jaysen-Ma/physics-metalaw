"""Formula catalog for the HTML ledger.

The Constantin release lives at
`.scratch/physics-formulae-patterns/data/` and is never edited.
`python -m physics_metalaw.catalog --seed` rebuilds
`data/catalog/{feynman,wiki,inflation}.json` from those CSVs.

Source URLs, by corpus:
  feynman   — one HTML chapter per id, `https://www.feynmanlectures.caltech.edu/{I,II,III}_{cc}.html`
  wiki      — per-equation Wikipedia page once identified; until then the named-equations list
  inflation — one PDF for every model, arXiv:1303.3787; locator is the model acronym
  feynman_lectures_v1/v2/v3 — Feynman Lectures chapter HTML, one link per equation
              ((c) Caltech; shown with source links; scripts/build_feynman_lectures_catalog.py)
"""
from __future__ import annotations

import csv
import json
import re
from collections import Counter
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Callable, Iterable

from . import metalaw as M

ROOT = Path(__file__).resolve().parents[1]
CATALOG_DIR = ROOT / "data" / "catalog"
RELEASED_CLONE = ROOT / ".scratch" / "physics-formulae-patterns"
RELEASED_DIR = RELEASED_CLONE / "data"
RELEASED_REPO_URL = "https://github.com/andrei-const/physics-formulae-patterns"


def released_clone_commit(clone: Path | None = None) -> str | None:
    """HEAD of the authors' CSV clone. None if `.scratch/` was not cloned."""
    root = clone or RELEASED_CLONE
    git_dir = root / ".git"
    if git_dir.is_file():
        spec = git_dir.read_text(encoding="utf-8").strip()
        if spec.lower().startswith("gitdir:"):
            git_dir = Path(spec.split(":", 1)[1].strip())
            if not git_dir.is_absolute():
                git_dir = root / git_dir
    head_file = git_dir / "HEAD"
    if not git_dir.is_dir() or not head_file.is_file():
        return None
    raw = head_file.read_text(encoding="utf-8").strip()
    if raw.startswith("ref:"):
        ref = git_dir / raw.split(":", 1)[1].strip()
        if not ref.is_file():
            return None
        raw = ref.read_text(encoding="utf-8").strip()
    return raw or None


def released_csv_meta(csv_name: str) -> dict:
    """Provenance stamped into each catalog JSON at seed time."""
    meta = {
        "released_csv": f"physics-formulae-patterns/data/{csv_name}",
        "released_repo": RELEASED_REPO_URL,
    }
    sha = released_clone_commit()
    if sha:
        meta["released_commit"] = sha
    return meta

WIKI_LIST_URL = (
    "https://en.wikipedia.org/wiki/List_of_scientific_equations_named_after_people"
)


def wiki_article_url(title: str) -> str:
    from urllib.parse import quote
    return "https://en.wikipedia.org/wiki/" + quote(
        title.replace(" ", "_"), safe=":_()'-*"
    )


INFLATION_ABS_URL = "https://arxiv.org/abs/1303.3787"
INFLATION_PDF_URL = "https://arxiv.org/pdf/1303.3787"
FEYNMAN_HOME_URL = "https://www.feynmanlectures.caltech.edu/"
FEYNMAN_CHAPTER_URL = (
    "https://www.feynmanlectures.caltech.edu/{volume}_{chapter:02d}.html"
)

# Independent-study slices. `released` is the paper's n=100, IDs notwithstanding.
Policy = Callable[["FormulaRecord"], bool]
POLICIES: dict[str, Policy] = {
    "released": lambda r: r.in_released,
    "primary": lambda r: r.representative,
    "lecture": lambda r: r.status == "original",
}


@dataclass
class FormulaRecord:
    corpus: str
    original_id: str
    id: str
    formula: str
    original_formula: str
    output: str | None
    n_variables: int
    n_numbers: int
    n_constants: int
    complexity: int
    original_n_variables: int
    original_n_numbers: int
    original_n_constants: int
    original_complexity: int
    status: str
    family: str | None
    family_role: str | None
    representative: bool
    in_released: bool
    source_kind: str
    source_url: str | None
    source_locator: str
    notes: str
    name: str = ""
    constant_names: list[str] = field(default_factory=list)
    variable_names: list[str] = field(default_factory=list)
    operator_counts: Counter = field(default_factory=Counter, repr=False, compare=False)
    consideration: str = ""
    branch: str = ""
    # False: formula text may not be released; publish counts + refs only.
    public: bool = True

    def to_json(self) -> dict:
        d = asdict(self)
        d.pop("operator_counts")
        constants = {"pi", "e", *self.constant_names}
        d["tree"] = M.expression_tree(self.formula, constants=constants)
        d["latex"] = M.formula_latex(self.formula)
        d["csv_disagreements"] = csv_disagreements(self)
        d.update(_unit_views(self))
        return d


# Universal / lecture constants. `pi` and `e` are already constants in
# `metalaw.expression_tree`; these extra names are otherwise painted as
# variables. `g_` is the Landé g-factor, not gravity, so it stays a variable.
PHYSICS_CONSTANTS = frozenset({
    "G", "c", "epsilon", "epsilon0", "h", "kb", "kB", "g", "R", "M_Pl",
})


def _name_occurrences(formula: str) -> dict[str, int]:
    return {k[5:]: v for k, v in M.count_operators(formula).items()
            if k.startswith("name:")}


def _view_payload(record: FormulaRecord, formula: str) -> dict:
    names = set(_name_occurrences(formula))
    constant_names = sorted(
        ((set(record.constant_names) | PHYSICS_CONSTANTS)
         - set(record.variable_names)) & names
    )
    constants = {"pi", "e", *constant_names}
    tree = M.expression_tree(formula, constants=constants)
    n_var = n_num = n_const = 0

    def walk(node):
        nonlocal n_var, n_num, n_const
        kind = node.get("kind")
        if kind == "variable":
            n_var += 1
        elif kind == "number":
            n_num += 1
        elif kind == "constant":
            n_const += 1
        for child in node.get("children") or []:
            walk(child)

    walk(tree)
    return {
        "formula": formula,
        "latex": M.formula_latex(formula),
        "tree": tree,
        "n_variables": n_var,
        "n_numbers": n_num,
        "n_constants": n_const,
        "complexity": M.tree_size(tree),
        "constant_names": constant_names,
    }


def rewrite_units(record: FormulaRecord, mode: str) -> str:
    """Drop painted physics constants. Names in `variable_names` stay.

    Modes:
      as_written — identity
      c1         — speed of light → 1
      hbar1      — unreduced h → 2πħ, then ħ = 1
      c1_hbar1   — both
      natural    — the paper's EI convention applied here: every painted
                   physics constant except M_Pl is set to 1 (ħ cancelled
                   via hbar1 when h is a constant)
    """
    if mode == "as_written":
        return record.formula
    consts = set(record.constant_names)
    formula = record.formula
    drop: set[str] = set()
    if mode in {"hbar1", "c1_hbar1", "natural"} and "h" in consts:
        formula = M.set_hbar_to_one(formula)
        consts.discard("h")
    if mode in {"c1", "c1_hbar1"}:
        drop.add("c")
    if mode == "natural":
        drop |= PHYSICS_CONSTANTS - {"M_Pl"}
    drop &= consts
    if drop:
        formula = M.set_names_to_one(formula, drop)
    return formula


def record_as_entry(record: FormulaRecord, mode: str = "as_written") -> M.EIEntry:
    """Catalog row as a metalaw entry, under a unit-rewrite mode."""
    formula = rewrite_units(record, mode)
    if formula == record.formula:
        ops = record.operator_counts
        return M.EIEntry(
            name=record.id, formula=formula,
            n_variables=record.n_variables, n_numbers=record.n_numbers,
            n_constants=record.n_constants,
            complexity_declared=record.complexity,
            operator_counts=ops, repaired=False)
    payload = _view_payload(record, formula)
    counts = M.count_operators(formula)
    ops = Counter(
        {k: v for k, v in counts.items()
         if not k.startswith("name:") and k != M.NULLARY_NUMBER}
    )
    return M.EIEntry(
        name=record.id, formula=formula,
        n_variables=payload["n_variables"],
        n_numbers=payload["n_numbers"],
        n_constants=payload["n_constants"],
        complexity_declared=payload["complexity"],
        operator_counts=ops, repaired=False)


def _unit_views(record: FormulaRecord) -> dict:
    """Ledger-only c=1 / ħ=1 rewrites. Keys omitted when that name is absent."""
    views = {}
    consts = set(record.constant_names)
    if "c" in consts:
        views["c1"] = _view_payload(record, rewrite_units(record, "c1"))
    if "h" in consts:
        views["hbar1"] = _view_payload(record, rewrite_units(record, "hbar1"))
    if "c" in consts and "h" in consts:
        views["c1_hbar1"] = _view_payload(
            record, rewrite_units(record, "c1_hbar1"))
    return views


def apply_constants(record: FormulaRecord) -> FormulaRecord:
    """Colour physics constants and recount n_variables / n_constants.

    Complexity is unchanged: a leaf that was a variable is still a leaf.
    `original_n_*` keep the CSV split. `g_` is not `g`.
    """
    if not record.formula.strip():
        return record
    names = set(_name_occurrences(record.formula))
    record.constant_names = sorted(
        ((set(record.constant_names) | PHYSICS_CONSTANTS)
         - set(record.variable_names)) & names
    )
    constants = {"pi", "e", *record.constant_names}
    n_var = n_const = 0
    for name, n in _name_occurrences(record.formula).items():
        if name in constants:
            n_const += n
        else:
            n_var += n
    record.n_variables = n_var
    record.n_constants = n_const
    return record


def csv_disagreements(record: FormulaRecord) -> dict[str, dict[str, int]]:
    """Catalog tally vs the released CSV. Empty iff they agree."""
    pairs = (
        ("variables", record.n_variables, record.original_n_variables),
        ("numbers", record.n_numbers, record.original_n_numbers),
        ("constants", record.n_constants, record.original_n_constants),
        ("complexity", record.complexity, record.original_complexity),
    )
    return {name: {"csv": orig, "catalog": cat}
            for name, cat, orig in pairs if cat != orig}


def attach_counts(record: FormulaRecord) -> FormulaRecord:
    if not record.formula.strip():
        record.operator_counts = Counter()
        return record
    counts = M.count_operators(record.formula)
    # CSV sometimes counts a collapsed Pow2/3/4 exponent as a number.
    record.n_numbers = counts.get(M.NULLARY_NUMBER, 0)
    record.operator_counts = Counter(
        {k: v for k, v in counts.items()
         if not k.startswith("name:") and k != M.NULLARY_NUMBER}
    )
    record.complexity = (sum(record.operator_counts.values())
                         + record.n_variables
                         + record.n_numbers
                         + record.n_constants)
    return record


def select(records: Iterable[FormulaRecord], policy: str) -> list[FormulaRecord]:
    try:
        keep = POLICIES[policy]
    except KeyError as e:
        raise ValueError(f"unknown policy {policy!r}; "
                         f"choose from {sorted(POLICIES)}") from e
    return [r for r in records if keep(r)]


def feynman_chapter_url(formula_id: str) -> str | None:
    """Chapter page for an id like `I.6.20`, `I.15.3x`, `II.6.15a`."""
    m = re.match(r"^(I{1,3})\.(\d+)\.", formula_id)
    if not m:
        return None
    return FEYNMAN_CHAPTER_URL.format(volume=m.group(1), chapter=int(m.group(2)))


def feynman_locator(formula_id: str) -> str:
    m = re.match(r"^I{1,3}\.(\d+\.\d+[a-z]?)", formula_id, re.I)
    return f"({m.group(1)})" if m else formula_id


# Overlay on the released CSV. Keys are original Filename values.
_FEYNMAN_PATCHES: dict[str, dict] = {
    "I.6.2a": {
        "id": "I.6.20a",
        "status": "original",
        "family": "I.6.20",
        "family_role": "special_case",
        "representative": False,
        "source_url": feynman_chapter_url("I.6.20a"),
        "source_locator": "(6.20) at σ=1; not a numbered equation",
        "notes": ("Unit-width Gaussian. σ=1 case of lecture (6.20); not a "
                  "numbered equation of its own. AI Feynman Table 4 already "
                  "used I.6.20a; the CSV dropped the zero."),
    },
    "I.6.2": {
        "id": "I.6.20",
        "status": "original",
        "family": "I.6.20",
        "family_role": "centred",
        "representative": False,
        "source_url": feynman_chapter_url("I.6.20"),
        "source_locator": "(6.20)",
        "notes": ("CSV id I.6.2 is lecture (6.2), the nuclear cross-section "
                  "N_C/N = nσ/A. This formula is lecture (6.20) in §6–4 — "
                  "the centred Gaussian as introduced. Family I.6.20 also "
                  "has the σ=1 case (a) and the mean-shifted form (b)."),
    },
    "I.6.2b": {
        "id": "I.6.20b",
        "status": "original",
        "family": "I.6.20",
        "family_role": "generalization",
        "representative": True,
        "source_url": feynman_chapter_url("I.6.20b"),
        "source_locator": "mean-shifted (6.20)",
        "notes": ("θ → θ-θ₁ shift of lecture (6.20). Working default for "
                  "this family (k=17)."),
    },
    "I.9.18": {
        "constant_names": ["G"],
        "status": "original",
        "family": "I.9.18",
        "family_role": "expanded",
        "representative": False,
        "notes": ("Newton gravity with r² written out as coordinate "
                  "differences. G is Newton's constant. Default family "
                  "member is I.9.18b, G m1 m2 / r²."),
    },
    "I.10.7": {
        "constant_names": ["c"],
        "notes": ("Relativistic mass factor. c is the speed of light, not a "
                  "variable. CSV already has 1 constant / 2 variables "
                  "(m_0, v); the tree now colours c as constant."),
    },
    "I.12.2": {
        "constant_names": ["epsilon"],
        "status": "original",
        "family": "I.12.2",
        "family_role": "vector",
        "representative": False,
        "source_locator": "(12.3); r/r³ vector encoding",
        "notes": ("Coulomb force. Lecture prose is 1/r²; the numbered "
                  "vector form is (12.3), F = q1 q2 r / (4 π ε0 r³). The "
                  "paper kept r/r³. Same choice as Newton: no vector, so "
                  "the family default is I.12.2b, q1 q2 / (4 π ε r²). "
                  "ε is a constant; CSV already has 2 constants (π, ε)."),
    },
    "I.12.4": {
        "constant_names": ["epsilon"],
        "status": "original",
        "family": "I.12.4",
        "family_role": "vector",
        "representative": False,
        "source_locator": "(12.4); r/r³ vector encoding",
        "notes": ("Electric field. Paper/CSV uses r/r³. Same as I.12.2: no "
                  "vector, so the family default is I.12.4b, "
                  "q1 / (4 π ε r²). ε is a constant."),
    },
    "I.13.12": {
        "notes": ("Gravitational potential energy. G is Newton's constant; "
                  "CSV already has 1 constant / 4 variables, the tree now "
                  "colours G as constant."),
    },
    "I.14.3": {
        "notes": ("Gravitational potential energy m g z. g is the "
                  "gravitational acceleration, not a variable. CSV already "
                  "has 1 constant / 2 variables (m, z)."),
    },
    "I.24.6": {
        "formula": "1/2*m*(omega**2+omega_0**2)*(1/2*x**2)",
        "notes": ("CSV writes *1/2*x**2, which Python parses as (...*1)/2 "
                  "times x², so LaTeX nested the second half over the "
                  "whole product. Grouped as (1/2 x²); complexity unchanged."),
    },
    "I.30.3": {
        "formula": "I_0*sin(n*theta/2)**2/sin(theta/2)**2",
        "output": "I",
        "notes": ("CSV uses Int / Int_0 for intensity. Catalog uses I / I_0."),
    },
    "I.29.4": {
        "notes": ("Wave number k = ω/c. CSV counted c as a variable "
                  "(2 variables, 0 constants); c is the speed of light."),
    },
    "I.39.1": {
        "id": "I.39.10",
        "status": "original",
        "formula": "3/2*P_F*V",
        "source_locator": "(39.10)",
        "notes": ("CSV id I.39.1 is lecture (39.1), P = F/A. This formula "
                  "is lecture (39.10), PV = (2/3) U, written as "
                  "U = (3/2) p V. Same dropped-zero as I.6.2 → I.6.20. "
                  "CSV `pr` is pressure P_F."),
    },
    "I.39.11": {
        "formula": "1/(gamma-1)*P_F*V",
        "notes": "CSV `pr` is pressure P_F.",
    },
    "I.39.22": {
        "output": "P_F",
        "notes": "CSV output `pr` is pressure P_F.",
    },
    "I.43.16": {
        "notes": ("CSV `mu_drift` / `Volt` typeset as μ_drift / V_e."),
    },
    "I.43.31": {
        "notes": "CSV `mob` (mobility) typeset as μ.",
    },
    "I.47.23": {
        "notes": "CSV `pr` is pressure P. Output c is the speed of sound.",
    },
    "II.6.15a": {
        "status": "original",
        "family": "II.6.15",
        "family_role": "cartesian",
        "representative": False,
        "notes": ("Dipole field in cartesian form. Family II.6.15; working "
                  "default is II.6.15b (polar)."),
    },
    "II.6.15b": {
        "status": "original",
        "family": "II.6.15",
        "family_role": "generalization",
        "representative": True,
        "notes": ("Dipole field in polar form; also in FLP II.6. Working "
                  "default for this family."),
    },
}


def seed_feynman(csv_path: Path | None = None) -> list[FormulaRecord]:
    path = csv_path or RELEASED_DIR / "Feynman.csv"
    with path.open(newline="") as fh:
        rows = [r for r in csv.reader(fh, delimiter=";")
                if any(c.strip() for c in r)][1:]
    records = []
    for row in rows:
        original_id = row[0].strip()
        formula = row[3]
        if not formula.strip():
            continue
        n_var, n_num, n_const, k = (int(row[4]), int(row[5]),
                                    int(row[6]), int(row[7]))
        rec = FormulaRecord(
            corpus="feynman",
            original_id=original_id,
            id=original_id,
            formula=formula,
            original_formula=formula,
            output=row[2].strip() or None,
            n_variables=n_var, n_numbers=n_num, n_constants=n_const,
            complexity=k,
            original_n_variables=n_var, original_n_numbers=n_num,
            original_n_constants=n_const, original_complexity=k,
            status="original",
            family=None, family_role=None,
            representative=True, in_released=True,
            source_kind="html",
            source_url=feynman_chapter_url(original_id),
            source_locator=feynman_locator(original_id),
            notes="",
        )
        for key, value in _FEYNMAN_PATCHES.get(original_id, {}).items():
            setattr(rec, key, value)
        records.append(attach_counts(apply_constants(rec)))
    records.extend(attach_counts(apply_constants(r)) for r in _extra_feynman())
    return records


# CSV glued a letter onto its subscript. Catalog Python identifiers use
# an underscore so the tree and LaTeX both show a subscript. `Re`/`Pr`
# stay as Reynolds/Prandtl; `dV`/`dQ`/`dG`/`DE` are differentials / Δ.
_WIKI_NAME_RENAMES: tuple[tuple[str, str], ...] = (
    ("fD", "f_D"),
    ("Dv", "D_v"),
    ("vp", "v_p"),
    ("vs", "v_s"),
    ("Ea", "E_a"),
    ("Zn", "Z_n"),
    ("Tm", "T_m"),
)


def _rename_wiki_names(formula: str) -> str:
    for old, new in _WIKI_NAME_RENAMES:
        formula = re.sub(rf"\b{old}\b", new, formula)
    return formula


_WIKI_PATCHES: dict[str, dict] = {
    # R is the molar gas constant on these rows, not a radius / resistance.
    # Wiki_31's R0 is Callendar resistance.
    "Wiki_7": {
        "constant_names": ["R"],
        "notes": ("Arrhenius k = A exp(-E_a/(R T)). R is the gas constant; "
                  "CSV already 3 variables / 1 constant."),
    },
    "Wiki_8": {
        "constant_names": ["R"],
        "notes": "Modified Arrhenius. R is the gas constant.",
    },
    "Wiki_9": {
        "constant_names": ["R"],
        "notes": "Stretched Arrhenius. R is the gas constant.",
    },
    "Wiki_10": {
        "formula": "(kb*T/h)*exp(dG/(R*T))",
        "constant_names": ["R"],
        "notes": ("Eyring. CSV `kB` is Boltzmann's constant; catalog uses "
                  "Feynman `kb` (k_b). R is the gas constant."),
    },
    "Wiki_11": {
        "formula": "-kb*T*log(Zn)",
        "notes": ("Potts Helmholtz A_n = -kT log Z_n. CSV `k` is Boltzmann's "
                  "constant; catalog uses Feynman `kb` (k_b). Avrami Wiki_12 "
                  "keeps `k` as a rate."),
    },
    "Wiki_13": {
        "constant_names": ["R"],
        "variable_names": ["c"],
        "notes": ("Benedict–Webb–Rubin. R is the gas constant (three "
                  "occurrences). CSV `c` is a BWR parameter, not the "
                  "speed of light."),
    },
    "Wiki_16": {
        "variable_names": ["h"],
        "notes": ("Bernoulli / hydrostatic q + ρ g h. `h` is a height, not "
                  "Planck's constant; `g` is gravity."),
    },
    "Wiki_26": {
        "formula": "A*j**(-n)*exp(Q/(kb*T))",
        "notes": ("Black's equation. CSV `k` in Q/(k T) is Boltzmann's "
                  "constant; catalog uses Feynman `kb`."),
    },
    "Wiki_40": {
        "formula": "exp(-(z**2*q**2*kappa)/(8*pi*epsilonr*epsilon0*kb*T))",
        "notes": "Debye–Hückel. CSV `kB` renamed to Feynman `kb`.",
    },
}


# Per-row Wikipedia article. Locator is the article title (link text);
# `name` is the specific formula on that page.
_WIKI_PAGES: dict[str, tuple[str, str]] = {
    "Wiki_1": ("Adams–Williamson equation", "compressional-wave speed v_p"),
    "Wiki_2": ("Adams–Williamson equation", "shear-wave speed v_s"),
    "Wiki_3": ("Adams–Williamson equation", "seismic parameter Φ"),
    "Wiki_4": ("Archard equation", "wear volume Q"),
    "Wiki_5": ("Archard equation", "δW = P π a²"),
    "Wiki_6": ("Archard equation", "δV = (2/3) π a³"),
    "Wiki_7": ("Arrhenius equation", "Arrhenius rate"),
    "Wiki_8": ("Arrhenius equation", "modified Arrhenius"),
    "Wiki_9": ("Arrhenius equation", "stretched Arrhenius"),
    "Wiki_10": ("Eyring equation", "Eyring rate"),
    "Wiki_11": ("Potts model", "Helmholtz A = −kT log Z"),
    "Wiki_12": ("Avrami equation", "Avrami / KJMA"),
    "Wiki_13": ("Benedict–Webb–Rubin equation", "BWR equation of state"),
    "Wiki_14": ("Bernoulli's principle", "Bernoulli energy"),
    "Wiki_15": ("Bernoulli's principle", "steady-flow Bernoulli"),
    "Wiki_16": ("Bernoulli's principle", "q + ρgh"),
    "Wiki_17": ("Bernoulli's principle", "hydrostatic p₀ + ρgz"),
    "Wiki_18": ("Bernoulli's principle", "dynamic pressure ρv²/2"),
    "Wiki_19": ("Bernoulli's principle", "piezometric head"),
    "Wiki_20": ("Bernoulli's principle", "Bernoulli head"),
    "Wiki_21": ("Bernoulli's principle", "compressible Bernoulli"),
    "Wiki_22": ("Bernoulli's principle", "Bernoulli + enthalpy"),
    "Wiki_23": ("Birch–Murnaghan equation of state", "third-order P(V)"),
    "Wiki_24": ("Birch–Murnaghan equation of state", "Eulerian strain f"),
    "Wiki_25": ("Birch–Murnaghan equation of state", "internal energy E(V)"),
    "Wiki_26": ("Black's equation", "MTTF electromigration"),
    "Wiki_27": ("Blaney–Criddle equation", "reference evapotranspiration"),
    "Wiki_28": ("Borda–Carnot equation", "head loss ΔE"),
    "Wiki_29": ("Borda–Carnot equation", "ΔH = ΔE/(ρg)"),
    "Wiki_30": ("Schrödinger equation", "H = p²/2m + V"),
    "Wiki_31": ("Callendar–Van Dusen equation", "platinum resistance"),
    "Wiki_32": ("Carothers equation", "degree of polymerization"),
    "Wiki_33": ("Cauchy's equation", "n(λ) dispersion"),
    "Wiki_34": ("Cesàro equation", "κ = C/s"),
    "Wiki_35": ("Cesàro equation", "κ = C/√s"),
    "Wiki_36": ("Cesàro equation", "κ = Cs"),
    "Wiki_37": ("Cesàro equation", "κ = a/(s²+a²)"),
    "Wiki_38": ("Churchill–Bernstein equation", "Nusselt number"),
    "Wiki_39": ("Darcy–Weisbach equation", "pipe pressure drop"),
    "Wiki_40": ("Debye–Hückel theory", "activity coefficient"),
    "Wiki_41": ("Doppler effect", "moving source, stationary observer"),
}


def seed_wiki(csv_path: Path | None = None) -> list[FormulaRecord]:
    path = csv_path or RELEASED_DIR / "Wiki.csv"
    with path.open(newline="") as fh:
        rows = [r for r in csv.reader(fh, delimiter=";")
                if any(c.strip() for c in r)][1:]
    records = []
    for row in rows:
        original_id = row[0].strip()
        formula = row[2]
        if not original_id or not formula.strip():
            continue
        n_var, n_num, n_const, k = (int(row[3]), int(row[4]),
                                    int(row[5]), int(row[6]))
        rec = FormulaRecord(
            corpus="wiki",
            original_id=original_id,
            id=original_id,
            formula=formula,
            original_formula=formula,
            output=None,
            n_variables=n_var, n_numbers=n_num, n_constants=n_const,
            complexity=k,
            original_n_variables=n_var, original_n_numbers=n_num,
            original_n_constants=n_const, original_complexity=k,
            status="original",
            family=None, family_role=None,
            representative=True, in_released=True,
            source_kind="html",
            source_url=WIKI_LIST_URL,
            source_locator=row[1].strip(),
            notes="",
        )
        page = _WIKI_PAGES.get(original_id)
        if page:
            article, name = page
            rec.name = name
            rec.source_url = wiki_article_url(article)
            rec.source_locator = article
        for key, value in _WIKI_PATCHES.get(original_id, {}).items():
            setattr(rec, key, value)
        rec.formula = _rename_wiki_names(rec.formula)
        records.append(attach_counts(apply_constants(rec)))
    records.extend(_wiki_extension_records(start=len(records) + 1))
    records.extend(_wiki_qm_records(start=len(records) + 1))
    return records


# Named-equations list after Doppler, same cut as the released 41:
# an explicit algebraic defining expression, no derivative or integral,
# no implicit equation, no special-function ODE. Not part of Wiki.csv.
# `variable_names` stops a painted constant (R, G, …) stealing a symbol
# that the article uses for something else.
_WIKI_EXTENSION = (
    ("Drake equation", "number of communicative civilisations",
     "R_star*f_p*n_e*f_l*f_i*f_c*L", (), ()),
    ("Fenske equation", "minimum plates at total reflux",
     "log((X_d/(1-X_d))*((1-X_b)/X_b))/log(alpha)", (), ()),
    ("Fresnel equations", "s-polarised reflectance",
     "((n_1*cos(theta_i)-n_2*cos(theta_t))/(n_1*cos(theta_i)+n_2*cos(theta_t)))**2",
     (), ()),
    ("Fresnel equations", "p-polarised reflectance",
     "((n_1*cos(theta_t)-n_2*cos(theta_i))/(n_1*cos(theta_t)+n_2*cos(theta_i)))**2",
     (), ()),
    ("Goldman equation", "resting potential, Na K Cl",
     "(R*T/F)*log((P_Na*Na_out+P_K*K_out+P_Cl*Cl_in)/(P_Na*Na_in+P_K*K_in+P_Cl*Cl_out))",
     (), ("F",)),
    ("Hadamard–Rybczynski equation", "bubble terminal speed",
     "(2/3)*R**2*g*(rho_b-rho_0)/mu_0*(mu_0+mu_b)/(2*mu_0+3*mu_b)",
     ("R",), ()),
    ("Hammett equation", "log(K/K0) = σρ",
     "sigma*rho", (), ()),
    ("Hazen–Williams equation", "pipe flow speed",
     "k*C*R**0.63*S**0.54", ("R",), ()),
    ("Henderson–Hasselbalch equation", "pH from pKa",
     "pK_a+log(Base/Acid)", (), ()),
    ("Hill equation (biochemistry)", "fractional occupancy",
     "L**n/(K_d+L**n)", (), ()),
    ("Kapustinskii equation", "lattice energy",
     "K*nu*abs(z_p)*abs(z_m)/(r_p+r_m)*(1-d/(r_p+r_m))", (), ()),
    ("Karplus equation", "vicinal coupling",
     "A*cos(2*phi)+B*cos(phi)+C", (), ()),
    ("Kozeny–Carman equation", "pressure drop through a packed bed",
     "150*mu/(Phi**2*d_p**2)*((1-eps)**2)/(eps**3)*V_0", (), ()),
    ("Leonard–Merritt mass estimator", "mass inside radius R",
     "(16/(3*pi*G))*R*(2*V_R**2+V_T**2)", ("R",), ()),
    ("Lorentz–Lorenz equation", "molar refractivity",
     "N*alpha/(3*epsilon0)", (), ()),
    ("Mark–Houwink equation", "intrinsic viscosity",
     "K*M**a", (), ()),
    ("Michaelis–Menten kinetics", "enzyme rate",
     "V_max*a/(K_m+a)", (), ()),
    ("Nernst equation", "reduction potential",
     "E_0-(R*T)/(z*F)*log(a_Red/a_Ox)", (), ("F",)),
    ("Penman–Monteith equation", "reference evapotranspiration",
     "(0.408*Delta*(R_n-G_soil)+(900/T)*gamma*u_2*delta_e)/(Delta+gamma*(1+0.34*u_2))",
     (), ()),
    ("Prony equation", "hydraulic head loss",
     "(L/D)*(a*V+b*V**2)", (), ()),
    ("Rankine–Hugoniot conditions", "momentum flux",
     "rho*u**2+p", (), ()),
    ("Rankine–Hugoniot conditions", "stagnation enthalpy",
     "h+u**2/2", ("h",), ()),
    ("Saha ionization equation", "ionisation ratio",
     "(2/lambda_th**3)*(g_ip1/g_i)*exp(-(eps_ip1-eps_i)/(kb*T))", (), ()),
    ("Sakuma–Hattori equation", "radiation-thermometer signal",
     "C/(exp(c_2/(lambda_x*T))-1)", (), ()),
    ("Sellmeier equation", "refractive index squared",
     "1+B_1*lambda**2/(lambda**2-C_1)+B_2*lambda**2/(lambda**2-C_2)+B_3*lambda**2/(lambda**2-C_3)",
     (), ()),
    ("Souders–Brown equation", "flooding velocity",
     "k*sqrt((rho_L-rho_V)/rho_V)", (), ()),
    ("Starling equation", "transcapillary volume flux",
     "L_p*S*((P_c-P_i)-sigma*(Pi_p-Pi_g))", (), ()),
    ("Stokes–Einstein equation", "diffusion coefficient",
     "kb*T/(6*pi*eta*r)", (), ()),
    ("Taft equation", "log(k/k_CH3) = ρ*σ* + δ E_s",
     "rho_star*sigma_star+delta*E_s", (), ()),
    ("Torricelli's law", "efflux speed",
     "sqrt(2*g*h)", ("h",), ()),
    ("Tsiolkovsky rocket equation", "delta-v",
     "v_e*log(m_0/m_f)", (), ()),
    ("Van Deemter equation", "plate height",
     "A+B/u+C*u", (), ()),
    ("Van der Waals equation", "pressure",
     "R*T/(v-b)-a/v**2", (), ()),
    ("Young–Dupré equation", "Young: γ_SL + γ_LG cos θ",
     "gamma_SL+gamma_LG*cos(theta)", (), ()),
    ("Young–Laplace equation", "pressure jump",
     "-gamma*(1/R_1+1/R_2)", (), ()),
)


def _wiki_extension_records(start: int) -> list[FormulaRecord]:
    note = ("Algebraic defining expression from the named-equations list "
            "after Doppler. Not in the released 41.")
    records = []
    for i, (article, name, formula, variables, constants) in enumerate(_WIKI_EXTENSION):
        wid = f"Wiki_{start + i}"
        rec = FormulaRecord(
            corpus="wiki",
            original_id=wid,
            id=wid,
            formula=formula,
            original_formula=formula,
            output=None,
            n_variables=0, n_numbers=0, n_constants=0, complexity=0,
            original_n_variables=0, original_n_numbers=0,
            original_n_constants=0, original_complexity=0,
            status="extension",
            family=None, family_role=None,
            representative=True, in_released=False,
            source_kind="html",
            source_url=wiki_article_url(article),
            source_locator=article,
            notes=note,
            name=name,
            constant_names=list(constants),
            variable_names=list(variables),
        )
        rec = attach_counts(apply_constants(rec))
        rec.original_n_variables = rec.n_variables
        rec.original_n_numbers = rec.n_numbers
        rec.original_n_constants = rec.n_constants
        rec.original_complexity = rec.complexity
        records.append(rec)
    return records


# Defining equations of quantum mechanics, one per law, from the Wikipedia
# article. Equality is the root and is not an operator. A differentiation
# variable, an integral measure, and a domain label are not nodes.
# `I` is the imaginary unit. Not part of the released 41 or the algebraic
# extension, and not drawn on the Wikipedia scatter.
_WIKI_QM = (
    ("Schrödinger equation", "time-dependent Schrödinger equation",
     "I*hbar*Partial(Psi) = -hbar**2/(2*m)*Partial2(Psi) + V*Psi",
     (), ("I", "hbar")),
    ("Schrödinger equation", "time-independent Schrödinger equation",
     "E*psi = -hbar**2/(2*m)*Partial2(psi) + V*psi",
     (), ("hbar",)),
    ("Dirac equation", "Dirac equation",
     "I*hbar*Dot(alpha, Grad(psi)) + beta*m*c*psi = I*hbar/c*Partial(psi)",
     (), ("I", "hbar")),
    ("Klein–Gordon equation", "Klein–Gordon equation",
     "Partial2(phi)/c**2 - Laplacian(phi) + m**2*c**2/hbar**2*phi = 0",
     (), ("hbar",)),
    ("Pauli equation", "Pauli equation",
     "Dot(sigma, p - q*A)**2/(2*m)*psi + q*phi*psi = I*hbar*Partial(psi)",
     (), ("I", "hbar")),
    ("Heisenberg picture", "Heisenberg equation of motion",
     "D(A) = I/hbar*Commutator(H, A) + Partial(A)",
     (), ("I", "hbar")),
    ("Canonical commutation relation", "canonical commutation relation",
     "Commutator(x, p) = I*hbar*id",
     (), ("I", "hbar", "id")),
    ("Probability current", "probability current",
     "j = -I*hbar/(2*m)*(conj(Psi)*Grad(Psi) - Psi*Grad(conj(Psi)))",
     (), ("I", "hbar")),
    ("Probability current", "probability continuity equation",
     "Divergence(j) + Partial(rho) = 0",
     (), ()),
    ("Ehrenfest theorem", "Ehrenfest theorem",
     "D(Expectation(A)) = Expectation(Commutator(A, H))/(I*hbar) + Expectation(Partial(A))",
     (), ("I", "hbar")),
    ("Path integral formulation", "path integral for the wavefunction",
     "psi = 1/Z*Integral(exp(I*S)*psi_0)",
     (), ("I",)),
    ("Fermi's golden rule", "Fermi's golden rule",
     "Gamma = 2*pi/hbar*abs(H_fi)**2*rho",
     (), ("hbar",)),
)


def _wiki_qm_records(start: int) -> list[FormulaRecord]:
    note = ("Defining equation from the Wikipedia article. Equality is the "
            "root and is not an operator. Differentiation variables and the "
            "path measure are binders, not nodes. Not in the released 41 "
            "or the algebraic extension.")
    records = []
    for i, (article, name, formula, variables, constants) in enumerate(_WIKI_QM):
        wid = f"Wiki_{start + i}"
        rec = FormulaRecord(
            corpus="wiki",
            original_id=wid,
            id=wid,
            formula=formula,
            original_formula=formula,
            output=None,
            n_variables=0, n_numbers=0, n_constants=0, complexity=0,
            original_n_variables=0, original_n_numbers=0,
            original_n_constants=0, original_complexity=0,
            status="qm",
            family=None, family_role=None,
            representative=True, in_released=False,
            source_kind="html",
            source_url=wiki_article_url(article),
            source_locator=article,
            notes=note,
            name=name,
            constant_names=list(constants),
            variable_names=list(variables),
        )
        rec = attach_counts(apply_constants(rec))
        rec.original_n_variables = rec.n_variables
        rec.original_n_numbers = rec.n_numbers
        rec.original_n_constants = rec.n_constants
        rec.original_complexity = rec.complexity
        records.append(rec)
    return records


# Paper names (Martin, Ringeval, Vennin, Encyclopædia Inflationaris).
# CSV Filename is the acronym; the blank Filename between SAII and FI is VFMI.
_EI_NAMES = {
    "SI": "Starobinsky Inflation",
    "RCHI": "Radiatively Corrected Higgs Inflation",
    "LFI": "Large Field Inflation",
    "MLFI": "Mixed Large Field Inflation",
    "RCMI": "Radiatively Corrected Massive Inflation",
    "RCQI": "Radiatively Corrected Quartic Inflation",
    "NI": "Natural Inflation",
    "ESI": "Exponential SUSY Inflation",
    "PLI": "Power Law Inflation",
    "KMII": "Kähler Moduli Inflation I",
    "HF1I": "Horizon Flow Inflation at first order",
    "CWI": "Coleman–Weinberg Inflation",
    "LI": "Loop Inflation",
    "RpI": "(R + R^{2p}) Inflation",
    "DWI": "Double-Well Inflation",
    "MHI": "Mutated Hilltop Inflation",
    "RGI": "Radion Gauge Inflation",
    "MSSMI": "MSSM Inflation",
    "RIPI": "Renormalizable Inflection Point Inflation",
    "AI": "Arctan Inflation",
    "CNAI": "Constant n_s A Inflation",
    "CNBI": "Constant n_s B Inflation",
    "OSTI": "Open String Tachyonic Inflation",
    "WRI": "Witten–O’Raifeartaigh Inflation",
    "CCSI": "Cubicly Corrected Starobinsky Inflation",
    "SBKI": "Symmetry Breaking Kähler Inflation",
    "AHI": "Axion Hilltop Inflation",
    "PAI": "Pure Arctan Inflation",
    "SAAI": "Superconformal α-attractor A Inflation",
    "TMI": "T-Model Inflation",
    "SFI": "Small Field Inflation",
    "II": "Intermediate Inflation",
    "KMIII": "Kähler Moduli Inflation II",
    "LMI": "Logamediate Inflation",
    "TWI": "Twisted Inflation",
    "GMSSMI": "Generalized MSSM Inflation",
    "GRIPI": "Generalized Renormalizable Point Inflation",
    "BSUSYBI": "Brane SUSY Breaking Inflation",
    "TI": "Tip Inflation",
    "PSNI": "Pseudo Natural Inflation",
    "NCKI": "Non-Canonical Kähler Inflation",
    "CSI": "Constant Spectrum Inflation",
    "OI": "Orientifold Inflation",
    "CNCI": "Constant n_s C Inflation",
    "SBI": "Supergravity Brane Inflation",
    "SSBI": "Spontaneous Symmetry Breaking Inflation",
    "IMI": "Inverse Monomial Inflation",
    "BI": "Brane Inflation",
    "SAII": "String Axion Inflation I",
    "VFMI": "Mukhanov Inflation",
    "FI": "Fibre Inflation",
    "HBI": "Hyperbolic Inflation",
    "SHI": "Smeared Higgs Inflation",
    "DEI": "Double Exponential Inflation",
    "SDI": "S-Dual Inflation",
    "GDWI": "Generalized Double Well Inflation",
    "SABI": "Superconformal α-attractor B Inflation",
    "SATI": "Superconformal α-attractor T Inflation",
    "RMI": "Running-mass Inflation",
    "VHI": "Valley Hybrid Inflation",
    "DSI": "Dynamical Supersymmetric Inflation",
    "GMLFI": "Generalized Mixed Large Field Inflation",
    "LPI": "Logarithmic Potential Inflation",
    "CNDI": "Constant n_s D Inflation",
    "SAIII": "String Axion Inflation II",
    "RCLFI": "Radiatively Corrected Large Field Inflation",
    "NCLI": "Non-Renormalizable Corrected Loop Inflation",
    "HNI": "Hybrid Natural Inflation",
    "NFI": "N-Formalism Inflation",
    "RCIPI": "Radiatively Corrected Inflection Point Inflation",
}

# Catalog vs CSV tally, resolved against Martin, Ringeval & Vennin
# arXiv:1303.3787v5. Constants are M_Pl, π, e, plus RCHI's A (CSV
# agreement). Other model parameters (p, μ, …) stay variables. Notes
# are appended in seed_inflation.
_EI_COUNT_NOTES = {
    "HF1I": (
        "Formula matches Martin Eq. (5.142). CSV complexity 27, tree 28 "
        "(hand-count off by one)."
    ),
    "CNAI": (
        "Formula matches Martin Eq. (5.291) after the bracket repair. "
        "CSV complexity 19, tree 20."
    ),
    "SFI": (
        "CSV counted 3 variables. Martin Eq. (6.1) has two parameters "
        "μ, p besides M and φ; catalog has 4 (a0=μ, a1=p). Complexity "
        "still 10 — they dropped the exponent leaf."
    ),
    "II": (
        "Martin Eq. (6.23) is two M^4 terms. CSV counted 8 variable "
        "instances; the tree has 9 (the second M). Complexity still 30."
    ),
    "SAII": (
        "Formula matches Martin Eq. (6.351). CSV complexity 19, tree 20."
    ),
    "VFMI": (
        "CSV 14/19/67 vs tree (see tally). Catalog follows the table: "
        "both (1+…) factors are 1 + (2-α)/(2 √3 β)·(φ/M_Pl), so x/M_Pl "
        "sits inside the sum (CSV had it multiplying the first (1+…) "
        "from the outside). original_formula is the CSV string."
    ),
    "FI": (
        "CSV complexity 64, tree 66. Nullary split matches. Formula left "
        "as CSV."
    ),
    "DEI": (
        "Formula matches Martin Eq. (6.419). CSV complexity 18, tree 19."
    ),
    "SABI": (
        "CSV counted 3 numbers, tree 4 (the extra 2 in ( )**(2 a1) or in "
        "sqrt(2/(3 a0))). Complexity still 21. Formula left as CSV."
    ),
    "VHI": (
        "Same missed-exponent slip as SFI. Martin Eq. (7.21) is "
        "M^4 [1+(φ/μ)^p]; catalog 4 variables, CSV 3. Complexity still 10."
    ),
    "DSI": (
        "Same missed-exponent slip as SFI. Martin Eq. (7.46) is "
        "M^4 [1+(φ/μ)^(-p)]; catalog 4 variables, CSV 3. Complexity still 11."
    ),
}

# CSV Filename blank, Number blank, between SAII (#48) and FI (#50).
_VFMI_CSV = (
    "(M**4)*(1-b*(2*(1+(2-a)/(2*sqrt(3)*b))*x/c)**(-1))*exp("
    "3*b/(1-a)*((1+((2-a)/2*sqrt(3)*b)*x/c)**(2*(1-a)/(2-a))-1))"
)
# Martin table / ASPIC, CSV `c` spelling. Applied before rename_ei_c.
# Status stays original; original_formula is the CSV string. Not used by
# metalaw.load_corpus (Result 1 still reads the released CSV).
_EI_FORMULA_PATCHES = {
    # Table + ASPIC rchisr: AI x / (16 π² √6), not A/π²·16·√6.
    "RCHI": (
        "(M**4)*(1-2*exp(-sqrt(2/3)*x/c)"
        "+(A/(16*pi**2*sqrt(6)))*(x/c))"
    ),
    # Both (1+…) factors are 1 + (2-α)/(2 √3 β)·(φ/M_Pl); CSV multiplied
    # the first (1+…) by x/c from the outside.
    "VFMI": (
        "(M**4)*(1-b*(2*(1+(2-a)/(2*sqrt(3)*b)*(x/c)))**(-1))*exp("
        "3*b/(1-a)*((1+(2-a)/(2*sqrt(3)*b)*(x/c))**(2*(1-a)/(2-a))-1))"
    ),
}
# Per-row catalog fields (kinds, not the Python string). Applied in
# seed_inflation before apply_constants.
_EI_RECORD_PATCHES = {
    "RCHI": {"constant_names": ["A"]},
}
_EI_C_NAME = re.compile(r"(?<![A-Za-z0-9_])c(?![A-Za-z0-9_])")


def rename_ei_c(formula: str) -> str:
    """CSV `c` is Martin's M_Pl, not the speed of light."""
    return _EI_C_NAME.sub("M_Pl", formula)


def seed_inflation(csv_path: Path | None = None) -> list[FormulaRecord]:
    """Seed Encyclopaedia Inflationaris from Inflation.csv. CSV is not edited.

    Catalog id is the CSV Filename. The blank Filename is VFMI. Six
    unbalanced-bracket rows get `metalaw.EI_BRACKET_REPAIRS` as `formula`;
    RCHI and VFMI get `_EI_FORMULA_PATCHES` (Martin's 16/√6; VFMI
    x/M_Pl inside both 1+(…) factors). RCHI's A is a constant.
    `original_formula` stays the CSV string. Status stays original.
    CSV `c` is renamed `M_Pl`.
    """
    path = csv_path or RELEASED_DIR / "Inflation.csv"
    with path.open(newline="") as fh:
        rows = [r for r in csv.reader(fh, delimiter=";")
                if any(c.strip() for c in r)][1:]
    records = []
    for row in rows:
        csv_id = row[0].strip()
        formula = row[2]
        if not formula.strip():
            continue
        notes = []
        if not csv_id:
            if formula != _VFMI_CSV:
                raise ValueError(
                    "blank Filename in Inflation.csv is not the VFMI row")
            rec_id = "VFMI"
            notes.append(
                "CSV Filename and Number blank between SAII (#48) and "
                "FI (#50); this is Mukhanov Inflation (VFMI).")
        else:
            rec_id = csv_id
        if rec_id == "PAI" and not row[1].strip():
            notes.append("CSV Number blank.")
        if rec_id in _EI_FORMULA_PATCHES:
            catalog_formula = rename_ei_c(_EI_FORMULA_PATCHES[rec_id])
            notes.append(
                "Catalog formula follows Martin arXiv:1303.3787; "
                "original_formula is the CSV string. Status original.")
            if rec_id == "RCHI":
                notes.append("A is a constant, agreeing with the CSV.")
        else:
            catalog_formula = rename_ei_c(
                M.EI_BRACKET_REPAIRS.get(rec_id, formula))
            if rec_id in M.EI_BRACKET_REPAIRS:
                notes.append(
                    "CSV has unbalanced brackets; catalog formula is "
                    "metalaw.EI_BRACKET_REPAIRS. Inflation.csv is not edited.")
        n_var, n_num, n_const, k = (int(row[3]), int(row[4]),
                                    int(row[5]), int(row[6]))
        csv_number = row[1].strip()
        locator = rec_id if not csv_number else f"{rec_id} (CSV #{csv_number})"
        rec = FormulaRecord(
            corpus="inflation",
            original_id=csv_id,
            id=rec_id,
            formula=catalog_formula,
            original_formula=formula,
            output="V",
            n_variables=n_var, n_numbers=n_num, n_constants=n_const,
            complexity=k,
            original_n_variables=n_var, original_n_numbers=n_num,
            original_n_constants=n_const, original_complexity=k,
            status="original",
            family=None, family_role=None,
            representative=True, in_released=True,
            source_kind="pdf",
            source_url=INFLATION_PDF_URL,
            source_locator=locator,
            notes=" ".join(notes),
            name=_EI_NAMES[rec_id],
        )
        for key, value in _EI_RECORD_PATCHES.get(rec_id, {}).items():
            setattr(rec, key, value)
        rec = attach_counts(apply_constants(rec))
        extra = _EI_COUNT_NOTES.get(rec.id)
        if extra:
            rec.notes = f"{rec.notes} {extra}".strip() if rec.notes else extra
        records.append(rec)
    return records


def _extra_feynman() -> list[FormulaRecord]:
    """Catalog-only rows that are not in Constantin's CSV."""
    return [
        FormulaRecord(
            corpus="feynman",
            original_id="I.9.18b",
            id="I.9.18b",
            formula="G*m1*m2/r**2",
            original_formula="G*m1*m2/r**2",
            output="F",
            n_variables=3, n_numbers=0, n_constants=1, complexity=8,
            original_n_variables=3, original_n_numbers=0,
            original_n_constants=1, original_complexity=8,
            status="rewritten",
            family="I.9.18",
            family_role="convention",
            representative=True,
            in_released=False,
            source_kind="html",
            source_url=feynman_chapter_url("I.9.18"),
            source_locator="G m1 m2 / r²; r not expanded",
            notes=("Rewritten I.9.18 to the same r² convention as the "
                   "other inverse-square rows, instead of "
                   "(x2-x1)²+(y2-y1)²+(z2-z1)². Working default for the "
                   "family (k=8). G is Newton's constant."),
            constant_names=["G"],
        ),
        FormulaRecord(
            corpus="feynman",
            original_id="I.12.2b",
            id="I.12.2b",
            formula="q1*q2/(4*pi*epsilon*r**2)",
            original_formula="q1*q2/(4*pi*epsilon*r**2)",
            output="F",
            n_variables=3, n_numbers=1, n_constants=2, complexity=12,
            original_n_variables=3, original_n_numbers=1,
            original_n_constants=2, original_complexity=12,
            status="rewritten",
            family="I.12.2",
            family_role="convention",
            representative=True,
            in_released=False,
            source_kind="html",
            source_url=feynman_chapter_url("I.12.2"),
            source_locator="scalar 1/r²; lecture prose, not the (12.3) vector",
            notes=("Rewritten Coulomb force to the scalar 1/r² convention "
                   "(same as Newton I.9.18b), instead of the CSV r/r³ "
                   "vector. Working default for the family (k=12). "
                   "ε is a constant."),
            constant_names=["epsilon"],
        ),
        FormulaRecord(
            corpus="feynman",
            original_id="I.12.4b",
            id="I.12.4b",
            formula="q1/(4*pi*epsilon*r**2)",
            original_formula="q1/(4*pi*epsilon*r**2)",
            output="Ef",
            n_variables=2, n_numbers=1, n_constants=2, complexity=10,
            original_n_variables=2, original_n_numbers=1,
            original_n_constants=2, original_complexity=10,
            status="rewritten",
            family="I.12.4",
            family_role="convention",
            representative=True,
            in_released=False,
            source_kind="html",
            source_url=feynman_chapter_url("I.12.4"),
            source_locator="scalar 1/r²; not the r/r³ vector",
            notes=("Rewritten electric-field magnitude to the same scalar "
                   "1/r² convention as I.12.2b, instead of the CSV r/r³ "
                   "vector. Working default for the family (k=10). "
                   "ε is a constant."),
            constant_names=["epsilon"],
        ),
    ]


def save_catalog(path: Path, records: list[FormulaRecord],
                 corpus: str, extra_meta: dict | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    meta = {
        "corpus": corpus,
        "n_records": len(records),
        "policies": sorted(POLICIES),
        "records": [r.to_json() for r in records],
    }
    if extra_meta:
        meta = {**extra_meta, **meta}
    path.write_text(json.dumps(meta, indent=2, ensure_ascii=False) + "\n")


def load_catalog(path: Path) -> list[FormulaRecord]:
    payload = json.loads(path.read_text())
    records = []
    for raw in payload["records"]:
        raw = dict(raw)
        raw.pop("tree", None)
        raw.pop("latex", None)
        raw.pop("c1", None)
        raw.pop("hbar1", None)
        raw.pop("c1_hbar1", None)
        raw.pop("csv_disagreements", None)
        for extra in ("chapter", "chapter_name", "section", "book_page", "pdf_page",
                      "equation", "handbook_latex", "operator_counts"):
            raw.pop(extra, None)
        known = FormulaRecord.__dataclass_fields__
        raw = {k: v for k, v in raw.items() if k in known}  # source-specific extras
        rec = FormulaRecord(**raw)
        records.append(attach_counts(apply_constants(rec)))
    return records


def feynman_path() -> Path:
    return CATALOG_DIR / "feynman.json"


def wiki_path() -> Path:
    return CATALOG_DIR / "wiki.json"


def inflation_path() -> Path:
    return CATALOG_DIR / "inflation.json"


def load_feynman() -> list[FormulaRecord]:
    path = feynman_path()
    if not path.exists():
        raise FileNotFoundError(
            f"{path} missing; run `python catalog.py --seed` with the "
            "released CSV cloned under .scratch/")
    return load_catalog(path)


def load_wiki() -> list[FormulaRecord]:
    path = wiki_path()
    if not path.exists():
        raise FileNotFoundError(
            f"{path} missing; run `python catalog.py --seed` with the "
            "released CSV cloned under .scratch/")
    return load_catalog(path)


def load_inflation() -> list[FormulaRecord]:
    path = inflation_path()
    if not path.exists():
        raise FileNotFoundError(
            f"{path} missing; run `python catalog.py --seed` with the "
            "released CSV cloned under .scratch/")
    return load_catalog(path)


SCATTER_CORPORA = (
    ("feynman", "Feynman"),
    ("wiki", "Wikipedia"),
    ("inflation", "Encyclopaedia Inflationaris"),
    ("leemis", "Leemis Univariate"),
    ("fr", "f(R)"),
)


# Feynman Lectures Vol. I-III catalogs. Shown on the Formulas page of the ledger
# (corpus button "Feynman Lectures"); not in the per-expression scatter, which keeps
# the five corpora above. The branch test (scripts/branch_test.py) pools them with the
# 100 benchmark expressions after deduplication.
FEYNMAN_LECTURES = ("feynman_lectures_v1", "feynman_lectures_v2", "feynman_lectures_v3")


def load_feynman_lectures(volume: int | None = None) -> list[FormulaRecord]:
    """Records of the Feynman Lectures catalogs (one volume, or all present)."""
    keys = FEYNMAN_LECTURES if volume is None else (f"feynman_lectures_v{volume}",)
    return [r for key in keys for r in scatter_records(key)]


def scatter_records(corpus: str) -> list[FormulaRecord]:
    """Expressions drawn on the ledger's β plots.

    Feynman and Inflationaris are the released slice. Wikipedia is the
    released 41 plus the algebraic extension after Doppler. Leemis is the
    laws under consideration. f(R) is the expressions held in that corpus.
    """
    if corpus in FEYNMAN_LECTURES:
        path = CATALOG_DIR / f"{corpus}.json"
        if not path.exists():
            return []
        return [r for r in load_catalog(path) if r.consideration == "in" and r.formula.strip()]
    recs = load_catalog(CATALOG_DIR / f"{corpus}.json")
    if corpus in {"leemis", "fr"}:
        return [r for r in recs if r.consideration == "in" and r.formula.strip()]
    if corpus == "wiki":
        return [r for r in recs
                if r.formula.strip() and (r.in_released or r.status == "extension")]
    return select(recs, "released")


def expression_scatter() -> dict:
    """One point per expression: complexity, and β fit on its own ranks.

    Each corpus carries the multinomial-MLE β of its pooled ranks (with
    bootstrap and Fisher errors), drawn as a horizontal line. That fit uses every formula in the corpus, including
    one that has too few ranks for a dot of its own.
    """
    from .beta import expression_beta

    corpora = []
    for key, label in SCATTER_CORPORA:
        records = scatter_records(key)
        if not records:
            continue
        points = []
        omitted = []
        entries = []
        for record in records:
            entry = record_as_entry(record, "as_written")
            entries.append(entry)
            beta, n_ranks = expression_beta(entry)
            row = {
                "id": record.id,
                "name": record.name,
                "k": int(entry.complexity_declared),
                "n_ranks": n_ranks,
            }
            if beta is None:
                omitted.append(row)
            else:
                row["beta"] = round(float(beta), 4)
                points.append(row)
        fit = M.beta_with_errors(entries)
        corpus_beta = fit["beta"]
        corpora.append({
            "id": key,
            "label": label,
            "n": len(points) + len(omitted),
            "n_fit": len(points),
            "omitted": omitted,
            "points": points,
            "corpus_beta": round(float(corpus_beta), 4),
            "corpus_beta_estimator": "multinomial MLE",
            "corpus_beta_boot_sigma": round(fit["boot_sigma"], 4),
            "corpus_beta_p16": round(fit["boot_p16"], 4),
            "corpus_beta_p84": round(fit["boot_p84"], 4),
            "corpus_beta_fisher_sigma": round(fit["fisher_sigma"], 4),
        })
    return {
        "note": (
            "Each dot is one expression. The x-axis is its complexity. "
            "β is the multinomial maximum-likelihood estimate: each operator "
            "instance is a draw that lands on rank r with probability "
            "p_r = exp(−β r) / sum_s exp(−β s), and β maximises "
            "sum_r n_r log p_r, where n_r is the number of instances at rank r, "
            "over the ranks that appear. The horizontal line is the same "
            "fit after every expression's operator counts in that corpus "
            "are added together and ranked once. Feynman and Inflationaris "
            "are the released formulas. Wikipedia adds the algebraic "
            "extension after Doppler. f(R) is the expressions held in that "
            "corpus."
        ),
        "corpora": corpora,
    }


_SCATTER_CACHE: dict | None = None


def expression_scatter_cached() -> dict:
    global _SCATTER_CACHE
    if _SCATTER_CACHE is None:
        _SCATTER_CACHE = expression_scatter()
    return _SCATTER_CACHE


def serve_catalog(host: str = "0.0.0.0", port: int = 8765) -> None:
    """Serve the HTML ledger. Bound on all interfaces so a Windows browser
    can use the machine IP if Cursor's 8765 forward is down. Loopback still
    works: `ssh -L 8765:127.0.0.1:8765 USER@HOST` then
    http://127.0.0.1:8765/
    """
    from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer

    class Handler(SimpleHTTPRequestHandler):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, directory=str(CATALOG_DIR), **kwargs)

        def do_GET(self):
            path = self.path.split("?", 1)[0].split("#", 1)[0]
            if path == "/expression-beta.json":
                body = json.dumps(expression_scatter_cached()).encode()
                self.send_response(200)
                self.send_header("Content-Type", "application/json; charset=utf-8")
                self.send_header("Content-Length", str(len(body)))
                self.end_headers()
                self.wfile.write(body)
                return
            if path in ("/formula-ledger.html", "/runs/formula-ledger.html"):
                path = "/index.html"
            self.path = path
            return super().do_GET()

        def end_headers(self):
            self.send_header("Cache-Control", "no-store")
            super().end_headers()

    httpd = ThreadingHTTPServer((host, port), Handler)
    print(f"catalog ledger  http://127.0.0.1:{port}/", flush=True)
    print(f"                bound {host}:{port}", flush=True)
    print("Windows (if Cursor is not forwarding the port):", flush=True)
    print(f"  ssh -L {port}:127.0.0.1:{port} USER@HOST", flush=True)
    print(f"  then http://127.0.0.1:{port}/", flush=True)
    httpd.serve_forever()


if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument("--seed", action="store_true",
                   help="rebuild data/catalog/{feynman,wiki,inflation}.json from the released CSVs")
    p.add_argument("--serve", action="store_true",
                   help="serve data/catalog/ on 0.0.0.0:8765 for the HTML ledger")
    p.add_argument("--port", type=int, default=8765)
    p.add_argument("--host", default="0.0.0.0",
                   help="bind address for --serve (default 0.0.0.0)")
    args = p.parse_args()
    if args.seed:
        seeds = [
            ("feynman", seed_feynman(), feynman_path(), {
                **released_csv_meta("Feynman.csv"),
                "source_kind": "html",
                "source_url_template": FEYNMAN_CHAPTER_URL,
                "board_url": FEYNMAN_HOME_URL,
                "board_label": "Feynman Lectures",
            }),
            ("wiki", seed_wiki(), wiki_path(), {
                **released_csv_meta("Wiki.csv"),
                "source_kind": "html",
                "source_url_note": "per-record Wikipedia article, not the named-equations list",
                "board_url": WIKI_LIST_URL,
                "board_label": "named-equations list",
            }),
            ("inflation", seed_inflation(), inflation_path(), {
                **released_csv_meta("Inflation.csv"),
                "source_kind": "pdf",
                "source_url": INFLATION_PDF_URL,
                "board_url": INFLATION_ABS_URL,
                "board_label": "arXiv:1303.3787",
                "M_Pl_is": "reduced Planck mass (CSV wrote c)",
                "M_is": "amplitude parameter (variable)",
            }),
        ]
        for name, recs, path, meta in seeds:
            save_catalog(path, recs, name, extra_meta=meta)
            print(f"wrote {path} ({len(recs)} records)")
            for pol in POLICIES:
                print(f"  {pol:10s} n={len(select(recs, pol)):3d}")
    if args.serve:
        serve_catalog(host=args.host, port=args.port)
    elif not args.seed:
        for name, loader in (("feynman", load_feynman),
                             ("wiki", load_wiki),
                             ("inflation", load_inflation)):
            try:
                recs = loader()
            except FileNotFoundError:
                continue
            print(name)
            for pol in POLICIES:
                sel = select(recs, pol)
                print(f"  {pol:10s} n={len(sel):3d}")
