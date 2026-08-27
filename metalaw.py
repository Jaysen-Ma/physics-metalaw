"""Operator-frequency statistics, following Constantin, Bartlett, Desmond &
Ferreira 2024 (arXiv:2408.11065).

The counting convention is not stated exhaustively in the paper, so it was
reverse-engineered from the authors' own released corpus
(github.com/andrei-const/physics-formulae-patterns) and validated against the
per-formula complexity column they publish. See `validate_ei()`.

Convention, as recovered:
  * An expression is a tree. Complexity is the total node count.
  * Nullary nodes split into three separate types for ranking purposes:
    variables, numbers (literals) and constants.
  * Binary: Add, Subtract, Mult, Div, Pow.
  * Unary: named functions (exp, log, sqrt, sin, cos, tanh, ...), unary minus,
    and Pow2/Pow3/Pow4 -- an integer literal exponent of 2, 3 or 4 collapses
    into a single unary operator and the literal is NOT counted separately.
"""
from __future__ import annotations

import ast
import csv
import re
from collections import Counter
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.optimize import curve_fit

BINARY = {ast.Add: "Add", ast.Sub: "Subtract", ast.Mult: "Mult",
          ast.Div: "Div", ast.Pow: "Pow"}
COLLAPSED_POWERS = (2, 3, 4)

NULLARY_VARIABLE = "variable"
NULLARY_NUMBER = "number"
NULLARY_CONSTANT = "constant"


def count_operators(expr: str, variables: set[str] | None = None) -> Counter:
    """Count operator instances in `expr`, written in Python/sympy syntax.

    Symbolic names land in `variable` if listed in `variables`, else `constant`.
    Pass `variables=None` to keep every symbolic name under its own key, which
    is what `validate_ei` needs since Encyclopaedia Inflationaris supplies its
    own nullary classification.
    """
    counts: Counter = Counter()

    def walk(node):
        if isinstance(node, ast.Expression):
            return walk(node.body)
        if isinstance(node, ast.BinOp):
            if (isinstance(node.op, ast.Pow)
                    and isinstance(node.right, ast.Constant)
                    and isinstance(node.right.value, int)
                    and node.right.value in COLLAPSED_POWERS):
                counts[f"Pow{node.right.value}"] += 1
                return walk(node.left)
            counts[BINARY[type(node.op)]] += 1
            walk(node.left)
            return walk(node.right)
        if isinstance(node, ast.UnaryOp):
            counts["Neg" if isinstance(node.op, ast.USub) else "Pos"] += 1
            return walk(node.operand)
        if isinstance(node, ast.Call):
            counts[node.func.id] += 1
            for arg in node.args:
                walk(arg)
            return None
        if isinstance(node, ast.Constant):
            counts[NULLARY_NUMBER] += 1
            return None
        if isinstance(node, ast.Name):
            if variables is None:
                counts["name:" + node.id] += 1
            else:
                counts[NULLARY_VARIABLE if node.id in variables
                       else NULLARY_CONSTANT] += 1
            return None
        raise TypeError(f"unhandled node: {ast.dump(node)}")

    walk(ast.parse(_sanitise(expr), mode="eval"))
    return counts


_RESERVED = ("lambda", "in", "is", "not", "and", "or", "if", "else")


def _sanitise(expr: str) -> str:
    """Rename symbol names that collide with Python keywords (Wikipedia's
    corpus uses `lambda` as a wavelength symbol)."""
    for kw in _RESERVED:
        expr = re.sub(rf"\b{kw}\b", kw + "_sym", expr)
    return expr


def complexity(expr: str) -> int:
    return sum(count_operators(expr).values())


def fit_beta(counts) -> tuple[float, float]:
    """Fit the normalised exponential f(r) = e^{-beta r} / sum_r' e^{-beta r'}
    to a rank-ordered frequency list. Returns (beta, sigma).

    The normalisation matters and is easy to get wrong: the authors' released
    notebook defines `exponential(x, beta) = np.exp(-beta*x)` without the
    denominator, while its Zipf model does include one. Fitting the
    unnormalised form to normalised frequencies gives beta ~ 1.07 for
    Encyclopaedia Inflationaris; the normalised form reproduces their published
    0.266 +/- 0.014. So the published beta is the normalised-model parameter.
    """
    obs = np.sort(np.asarray(list(counts), dtype=float))[::-1]
    freq = obs / obs.sum()
    rank = np.arange(len(freq)) + 1

    def model(r, beta):
        return np.exp(-beta * r) / np.sum(np.exp(-beta * rank))

    popt, pcov = curve_fit(model, rank, freq, p0=[0.3])
    return float(popt[0]), float(np.sqrt(pcov[0, 0]))


# Six Encyclopaedia Inflationaris entries in the authors' CSV have unbalanced
# brackets and cannot be parsed. The repairs below are the minimal bracket
# insertions/deletions that make each parse; none changes the intended model.
EI_BRACKET_REPAIRS = {
    "CNAI": "(M**4)*(3-(3+a**2)*(tanh(a*x/(c*sqrt(2)))**2))",
    "CNBI": "(M**4)*((3-a**2)*(tan(a*x/(c*sqrt(2)))**2-3))",
    "TMI": "(M**4)*(tanh(x/c*sqrt(6))**(4*a0))",
    "CNCI": "(M**4)*((3+a0**2)*(coth(a0*x/(c*sqrt(2)))**2-3))",
    "HBI": "(M**4)*((sinh(x/a0))**a1)",
    "SATI": "(M**4)*((tanh(x/(sqrt(6*a0)*c)))**(2*a1))",
}


@dataclass
class EIEntry:
    name: str
    formula: str
    n_variables: int
    n_numbers: int
    n_constants: int
    complexity_declared: int
    operator_counts: Counter   # non-nullary only
    repaired: bool

    @property
    def complexity_recomputed(self) -> int:
        return (sum(self.operator_counts.values())
                + self.n_variables + self.n_numbers + self.n_constants)


# column index of (name, formula, n_variables, n_numbers, n_constants, complexity)
CORPUS_COLUMNS = {
    "inflation": (0, 2, 3, 4, 5, 6),
    "feynman": (0, 3, 4, 5, 6, 7),
    "wiki": (0, 2, 3, 4, 5, 6),
}


def load_corpus(csv_path: str | Path, corpus: str) -> list[EIEntry]:
    """Load any of the three released corpora. Only `inflation` needs repairs."""
    ci_name, ci_formula, ci_var, ci_num, ci_const, ci_k = CORPUS_COLUMNS[corpus]
    rows = [r for r in csv.reader(open(csv_path), delimiter=";")
            if any(c.strip() for c in r)][1:]
    entries = []
    for row in rows:
        name = row[ci_name].strip() or "(unnamed)"
        formula = row[ci_formula]
        if not formula.strip():
            continue  # trailing column-totals row (Wikipedia corpus has one)
        if corpus == "inflation":
            formula = EI_BRACKET_REPAIRS.get(name, formula)
        counts = count_operators(formula)
        operators = Counter({k: v for k, v in counts.items()
                             if not k.startswith("name:") and k != NULLARY_NUMBER})
        entries.append(EIEntry(
            name=name, formula=formula,
            n_variables=int(row[ci_var]), n_numbers=int(row[ci_num]),
            n_constants=int(row[ci_const]), complexity_declared=int(row[ci_k]),
            operator_counts=operators,
            repaired=corpus == "inflation" and name in EI_BRACKET_REPAIRS))
    return entries


def load_ei(csv_path: str | Path) -> list[EIEntry]:
    """Load the Encyclopaedia Inflationaris corpus from the authors' CSV.

    Nullary counts come from the CSV's own three columns rather than from the
    formula text: the authors' split of symbolic names into `variables` and
    `constants` is not recoverable from the formulas alone (`p` in LFI is a
    variable, `A` in RCHI is a constant), so their classification is taken as
    given and only the operators are recounted.
    """
    rows = [r for r in csv.reader(open(csv_path), delimiter=";")
            if any(c.strip() for c in r)][1:]
    entries = []
    for row in rows:
        name = row[0].strip() or "(unnamed)"
        formula = EI_BRACKET_REPAIRS.get(name, row[2])
        counts = count_operators(formula)
        operators = Counter({k: v for k, v in counts.items()
                             if not k.startswith("name:") and k != NULLARY_NUMBER})
        entries.append(EIEntry(
            name=name, formula=formula,
            n_variables=int(row[3]), n_numbers=int(row[4]),
            n_constants=int(row[5]), complexity_declared=int(row[6]),
            operator_counts=operators, repaired=name in EI_BRACKET_REPAIRS))
    return entries


def arity_shares(entries: list[EIEntry]) -> dict[str, float]:
    """Fraction of operator instances that are nullary, unary and binary."""
    total = rank_counts(entries)
    nullary = sum(total[k] for k in
                  (NULLARY_VARIABLE, NULLARY_NUMBER, NULLARY_CONSTANT))
    binary = sum(v for k, v in total.items() if k in BINARY.values())
    unary = sum(total.values()) - nullary - binary
    n = sum(total.values())
    return {"nullary": nullary / n, "unary": unary / n, "binary": binary / n,
            "n_instances": n}


def rank_counts(entries: list[EIEntry]) -> Counter:
    """Aggregate to the ranked operator-instance counts the paper fits beta to."""
    total = Counter()
    for e in entries:
        total += e.operator_counts
        total[NULLARY_VARIABLE] += e.n_variables
        total[NULLARY_NUMBER] += e.n_numbers
        total[NULLARY_CONSTANT] += e.n_constants
    return total


ei_rank_counts = rank_counts
