"""Operator-frequency statistics, following Constantin, Bartlett, Desmond &
Ferreira 2024 (arXiv:2408.11065).

The counting convention is not stated exhaustively in the paper, so it was
reverse-engineered from the authors' own released corpus
(github.com/andrei-const/physics-formulae-patterns) and validated against the
per-formula complexity column they publish. `load_corpus` takes their
nullary columns as given; only operators are recounted.

Convention, as recovered:
  * An expression is a tree. Complexity is the total node count.
  * Nullary nodes split into three separate types for ranking purposes:
    variables, numbers (literals) and constants.
  * Binary: Add, Subtract, Mult, Div, Pow.
  * Unary: named functions (exp, log, sqrt, sin, cos, tanh, ...), unary minus,
    and Pow2/Pow3/Pow4 -- an integer literal exponent of 2, 3 or 4 collapses
    into a single unary operator and the literal is NOT counted separately.
  * An equation is two expressions under one equality. The equality is the
    tree root and is neither an operator nor a complexity node.
"""
from __future__ import annotations

import ast
import csv
import math
import re
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np

BINARY = {ast.Add: "Add", ast.Sub: "Subtract", ast.Mult: "Mult",
          ast.Div: "Div", ast.Pow: "Pow"}
COLLAPSED_POWERS = (2, 3, 4)

NULLARY_VARIABLE = "variable"
NULLARY_NUMBER = "number"
NULLARY_CONSTANT = "constant"
NULLARY_PARAMETER = "parameter"
NULLARY = (NULLARY_VARIABLE, NULLARY_NUMBER, NULLARY_CONSTANT, NULLARY_PARAMETER)


def _split_equation(expr: str) -> list[str]:
    """Split one equation on its top-level equals sign.

    A chain of equalities is rejected. `==` and comparison signs are not
    separators. An expression with no equals sign comes back as one piece.
    """
    cuts = []
    depth = 0
    for i, ch in enumerate(expr):
        if ch in "([{":
            depth += 1
        elif ch in ")]}":
            depth -= 1
        elif ch == "=" and depth == 0:
            prev = expr[i - 1] if i else ""
            nxt = expr[i + 1] if i + 1 < len(expr) else ""
            if prev not in "<>!" and nxt != "=":
                cuts.append(i)
    if not cuts:
        return [expr]
    if len(cuts) != 1:
        raise ValueError(f"equation chain, not one equality: {expr}")
    cut = cuts[0]
    return [expr[:cut].strip(), expr[cut + 1:].strip()]


def count_operators(expr: str, variables: set[str] | None = None) -> Counter:
    """Count operator instances in `expr`, written in Python/sympy syntax.

    Symbolic names land in `variable` if listed in `variables`, else `constant`.
    Pass `variables=None` to keep every symbolic name under its own key.
    Encyclopaedia Inflationaris supplies its own nullary classification, so
    `load_corpus` does not guess it from the formula text.
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

    for side in _split_equation(_sanitise(expr)):
        walk(ast.parse(side, mode="eval"))
    return counts


def expression_tree(expr: str, constants: set[str] | None = None) -> dict:
    """Same nodes as `count_operators`, nested for the HTML ledger.

    Each node is `{label, arity, kind, children}`. `kind` is `op` for
    operators, else `variable` / `number` / `constant`. Mathematical `pi`
    and `e` are constants by default; extra names (Newton's `G`, …) are
    passed in `constants`.
    """

    def leaf(label: str, kind: str) -> dict:
        return {"label": label, "arity": "nullary", "kind": kind, "children": []}

    constants = set(constants) if constants is not None else {"pi", "e"}

    def walk(node) -> dict:
        if isinstance(node, ast.Expression):
            return walk(node.body)
        if isinstance(node, ast.BinOp):
            if (isinstance(node.op, ast.Pow)
                    and isinstance(node.right, ast.Constant)
                    and isinstance(node.right.value, int)
                    and node.right.value in COLLAPSED_POWERS):
                return {
                    "label": f"Pow{node.right.value}",
                    "arity": "unary",
                    "kind": "op",
                    "children": [walk(node.left)],
                }
            return {
                "label": BINARY[type(node.op)],
                "arity": "binary",
                "kind": "op",
                "children": [walk(node.left), walk(node.right)],
            }
        if isinstance(node, ast.UnaryOp):
            return {
                "label": "Neg" if isinstance(node.op, ast.USub) else "Pos",
                "arity": "unary",
                "kind": "op",
                "children": [walk(node.operand)],
            }
        if isinstance(node, ast.Call):
            args = [walk(a) for a in node.args]
            return {
                "label": node.func.id,
                "arity": "unary" if len(args) == 1 else "binary",
                "kind": "op",
                "children": args,
            }
        if isinstance(node, ast.Constant):
            return leaf(str(node.value), "number")
        if isinstance(node, ast.Name):
            label = node.id[:-4] if node.id.endswith("_sym") else node.id
            kind = "constant" if label in constants else "variable"
            return leaf(label, kind)
        raise TypeError(f"unhandled node: {ast.dump(node)}")

    sides = [walk(ast.parse(side, mode="eval"))
             for side in _split_equation(_sanitise(expr))]
    if len(sides) == 1:
        return sides[0]
    return {"label": "=", "arity": "binary", "kind": "relation",
            "children": sides}


def tree_size(node: dict) -> int:
    """Counted nodes. An equality root is not an operator and adds nothing."""
    children = node.get("children") or ()
    here = 0 if node.get("kind") == "relation" else 1
    return here + sum(tree_size(c) for c in children)


_GREEK = (
    "alpha", "beta", "gamma", "delta", "epsilon", "theta", "lambda",
    "mu", "nu", "pi", "rho", "sigma", "tau", "phi", "chi", "psi",
    "omega", "kappa", "eta",
)
# Concatenated CSV identifiers that are one letter plus a subscript,
# or a conventional physics symbol the CSV spelled out in ASCII.
_SYMBOL_TEX_ALIASES = {
    "Ef": r"E_{f}",
    "hbar": r"\hbar",
    "Psi": r"\Psi",
    "lambd": r"\lambda",
    "kb": r"k_{b}",
    "kB": r"k_{b}",
    "pr": r"P",
    "Volt": r"V_{e}",
    "mob": r"\mu",
    "mom": r"\mu",
    "Pwr": r"P",
    "A_vec": r"\vec{A}",
    "epsilon0": r"\varepsilon_{0}",
    "epsilonr": r"\varepsilon_{r}",
    "M_Pl": r"M_{\mathrm{Pl}}",
    "Lambda1": r"\Lambda_{1}",
    "Lambda2": r"\Lambda_{2}",
}
_XYZ_COMP = re.compile(r"^([A-Za-z])([xyz])$")
_PREC_ADD, _PREC_MUL, _PREC_NEG, _PREC_POW, _PREC_ATOM = 1, 2, 3, 4, 5
_CALL_TEX = {
    "exp": lambda a: rf"e^{{{a}}}",
    "sqrt": lambda a: rf"\sqrt{{{a}}}",
    "sin": lambda a: rf"\sin\left({a}\right)",
    "cos": lambda a: rf"\cos\left({a}\right)",
    "tan": lambda a: rf"\tan\left({a}\right)",
    "tanh": lambda a: rf"\tanh\left({a}\right)",
    "log": lambda a: rf"\log\left({a}\right)",
    "ln": lambda a: rf"\ln\left({a}\right)",
    "arcsin": lambda a: rf"\arcsin\left({a}\right)",
    "arccos": lambda a: rf"\arccos\left({a}\right)",
    "arctan": lambda a: rf"\arctan\left({a}\right)",
    "abs": lambda a: rf"\left|{a}\right|",
    "sinh": lambda a: rf"\sinh\left({a}\right)",
    "cosh": lambda a: rf"\cosh\left({a}\right)",
    "sech": lambda a: rf"\operatorname{{sech}}\left({a}\right)",
    "coth": lambda a: rf"\coth\left({a}\right)",
    "zeta": lambda a: rf"\zeta\left({a}\right)",
    "Abs": lambda a: rf"\left|{a}\right|",
    "Gamma": lambda a: rf"\Gamma\left({a}\right)",
    "erf": lambda a: rf"\operatorname{{erf}}\left({a}\right)",
    "erfc": lambda a: rf"\operatorname{{erfc}}\left({a}\right)",
}


def formula_latex(expr: str, symbol_tex=None) -> str:
    """LaTeX for a Python/sympy formula, preserving written order.

    ``symbol_tex`` optionally maps a symbol name to LaTeX (return None to fall
    back to the default); corpora with their own naming conventions use it.
    """
    def name_tex(name: str) -> str:
        if symbol_tex is not None:
            tex = symbol_tex(name[:-4] if name.endswith("_sym") else name)
            if tex is not None:
                return tex
        return _symbol_tex(name)

    def paren(tex: str, inner: int, outer: int) -> str:
        return tex if inner >= outer else rf"\left({tex}\right)"

    def walk(node, outer: int) -> str:
        tex, inner = emit(node)
        return paren(tex, inner, outer)

    def emit(node) -> tuple[str, int]:
        if isinstance(node, ast.Expression):
            return emit(node.body)
        if isinstance(node, ast.BinOp):
            if (isinstance(node.op, ast.Pow)
                    and isinstance(node.right, ast.Constant)
                    and isinstance(node.right.value, int)
                    and node.right.value in COLLAPSED_POWERS):
                base_tex, inner = emit(node.left)
                return _pow_tex(base_tex, str(node.right.value), inner), _PREC_POW
            left_outer = {
                ast.Add: _PREC_ADD, ast.Sub: _PREC_ADD,
                ast.Mult: _PREC_MUL, ast.Div: 0, ast.Pow: _PREC_POW,
            }[type(node.op)]
            right_outer = {
                ast.Add: _PREC_ADD, ast.Sub: _PREC_MUL,
                ast.Mult: _PREC_MUL, ast.Div: 0, ast.Pow: _PREC_NEG,
            }[type(node.op)]
            if isinstance(node.op, ast.Pow):
                base_tex, inner = emit(node.left)
                exp_tex = walk(node.right, _PREC_NEG)
                return _pow_tex(base_tex, exp_tex, inner), _PREC_POW
            left, right = walk(node.left, left_outer), walk(node.right, right_outer)
            if isinstance(node.op, ast.Add):
                return rf"{left} + {right}", _PREC_ADD
            if isinstance(node.op, ast.Sub):
                return rf"{left} - {right}", _PREC_ADD
            if isinstance(node.op, ast.Mult):
                return _mul_tex(left, right), _PREC_MUL
            if isinstance(node.op, ast.Div):
                return rf"\frac{{{left}}}{{{right}}}", _PREC_ATOM
            raise TypeError(f"unhandled binop {type(node.op)}")
        if isinstance(node, ast.UnaryOp):
            if isinstance(node.op, ast.UAdd):
                return emit(node.operand)
            return rf"-{walk(node.operand, _PREC_NEG)}", _PREC_NEG
        if isinstance(node, ast.Call):
            tex = _call_tex(node, walk)
            if tex is not None:
                return tex, _PREC_ATOM
            args = ", ".join(walk(a, 0) for a in node.args)
            fmt = _CALL_TEX.get(node.func.id)
            tex = fmt(args) if fmt and len(node.args) == 1 else (
                rf"\operatorname{{{node.func.id}}}\left({args}\right)")
            return tex, _PREC_ATOM
        if isinstance(node, ast.Constant):
            return str(node.value), _PREC_ATOM
        if isinstance(node, ast.Name):
            return name_tex(node.id), _PREC_ATOM
        raise TypeError(f"unhandled node: {ast.dump(node)}")

    parts = [emit(ast.parse(side, mode="eval"))[0]
             for side in _split_equation(_sanitise(expr))]
    return " = ".join(parts)


def _call_tex(node, walk) -> str | None:
    """Rendered calculus and bracket operators. None leaves the usual call."""
    name = node.func.id
    args = [walk(a, _PREC_ATOM) for a in node.args]
    if name in {"Partial", "D", "Differential", "Grad"} and len(args) == 1:
        glyph = {"Partial": r"\partial ", "D": r"\mathrm{d} ",
                 "Differential": r"\mathrm{d} ", "Grad": r"\nabla "}[name]
        return glyph + args[0]
    if name == "Partial2" and len(args) == 1:
        return r"\partial^{2} " + args[0]
    if name == "D2" and len(args) == 1:
        return r"\mathrm{d}^{2} " + args[0]
    if name == "Laplacian" and len(args) == 1:
        return r"\nabla^{2} " + args[0]
    if name == "Divergence" and len(args) == 1:
        return r"\nabla\cdot " + args[0]
    if name == "Curl" and len(args) == 1:
        return r"\nabla\times " + args[0]
    if name == "conj" and len(args) == 1:
        return args[0] + "^{*}"
    if name == "Expectation" and len(args) == 1:
        return r"\langle " + args[0] + r" \rangle"
    if name == "Commutator" and len(args) == 2:
        return r"\left[" + args[0] + ", " + args[1] + r"\right]"
    if name == "Dot" and len(args) == 2:
        return args[0] + r" \cdot " + args[1]
    if name == "Cross" and len(args) == 2:
        return args[0] + r" \times " + args[1]
    if name == "Integral" and len(args) == 1:
        return r"\int " + args[0]
    if name == "Integral" and len(args) == 3:
        return (r"\int_{" + args[1] + "}^{" + args[2] + "} " + args[0])
    if name == "Sum" and len(args) == 1:
        return r"\sum " + args[0]
    if name == "Sum" and len(args) == 3:
        return (r"\sum_{" + args[1] + "}^{" + args[2] + "} " + args[0])
    if name == "delta" and len(args) == 1:
        return r"\delta\left(" + args[0] + r"\right)"
    return None


def _symbol_tex(name: str) -> str:
    if name.endswith("_sym"):
        name = name[:-4]
    if name in _SYMBOL_TEX_ALIASES:
        return _SYMBOL_TEX_ALIASES[name]
    bits = [b for b in name.split("_") if b]
    if not bits:
        return name

    def atom(part: str) -> str:
        if part in _SYMBOL_TEX_ALIASES:
            return _SYMBOL_TEX_ALIASES[part]
        m = re.match(r"^([A-Za-z]+)(\d+)$", part)
        if m:
            return atom(m.group(1)) + "_{" + m.group(2) + "}"
        xyz = _XYZ_COMP.match(part)
        if xyz:
            return xyz.group(1) + "_{" + xyz.group(2) + "}"
        low = part.lower()
        if low in _GREEK:
            return "\\" + low
        return part

    def sub(part: str) -> str:
        low = part.lower()
        if low in _GREEK:
            return "\\" + low
        if part.isdigit() or re.fullmatch(r"[A-Za-z]", part):
            return part
        if re.fullmatch(r"[A-Za-z]+\d+", part):
            return atom(part)
        return rf"\text{{{part}}}"

    if len(bits) == 1:
        return atom(bits[0])
    return atom(bits[0]) + "_{" + ",".join(sub(b) for b in bits[1:]) + "}"


def _pow_tex(base: str, exp: str, inner: int = _PREC_ATOM) -> str:
    if inner < _PREC_POW or base.startswith(r"\frac"):
        return rf"\left({base}\right)^{{{exp}}}"
    return rf"{base}^{{{exp}}}"


def _mul_tex(left: str, right: str) -> str:
    if left[-1:].isdigit() and right[:1].isdigit():
        return rf"{left} \cdot {right}"
    if right.startswith("-"):
        return rf"{left} \cdot \left({right}\right)"
    return rf"{left} {right}"


_RESERVED = ("lambda", "in", "is", "not", "and", "or", "if", "else")


def _sanitise(expr: str) -> str:
    """Rename symbol names that collide with Python keywords (Wikipedia's
    corpus uses `lambda` as a wavelength symbol)."""
    for kw in _RESERVED:
        expr = re.sub(rf"\b{kw}\b", kw + "_sym", expr)
    return expr


def complexity(expr: str) -> int:
    return sum(count_operators(expr).values())


def set_names_to_one(expr: str, names: Iterable[str] = ("c",)) -> str:
    """Replace named constants by 1 and drop the resulting trivial factors.

    `omega/c` → `omega`, `v**2/c**2` → `v**2`, `m*c**2` → `m`. The rest of
    the expression is not rearranged. For ħ=1 on unreduced `h`, use
    `set_hbar_to_one`.
    """
    names = set(names)
    if not names:
        return expr
    sides = _split_equation(expr)
    if len(sides) == 2:
        return "=".join(set_names_to_one(side, names) for side in sides)
    tree = ast.parse(_sanitise(expr), mode="eval")

    def is_one(node):
        return isinstance(node, ast.Constant) and node.value == 1

    def walk(node):
        if isinstance(node, ast.Expression):
            return walk(node.body)
        if isinstance(node, ast.Name) and node.id in names:
            return ast.Constant(value=1)
        if isinstance(node, ast.BinOp):
            left, right = walk(node.left), walk(node.right)
            if isinstance(node.op, ast.Mult):
                if is_one(left):
                    return right
                if is_one(right):
                    return left
            elif isinstance(node.op, ast.Div) and is_one(right):
                return left
            elif isinstance(node.op, ast.Pow):
                if is_one(left):
                    return ast.Constant(value=1)
                if is_one(right):
                    return left
            return ast.BinOp(left=left, op=node.op, right=right)
        if isinstance(node, ast.UnaryOp):
            operand = walk(node.operand)
            if isinstance(node.op, ast.UAdd):
                return operand
            return ast.UnaryOp(op=node.op, operand=operand)
        if isinstance(node, ast.Call):
            return ast.Call(
                func=node.func,
                args=[walk(a) for a in node.args],
                keywords=node.keywords,
            )
        return node

    return ast.unparse(walk(tree)).replace(" ", "")


def _is_one(node) -> bool:
    return isinstance(node, ast.Constant) and node.value == 1


def _const_int(node):
    if isinstance(node, ast.Constant) and isinstance(node.value, int):
        return node.value
    if isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.USub):
        v = _const_int(node.operand)
        return None if v is None else -v
    return None


def _flatten_mult(node):
    if isinstance(node, ast.BinOp) and isinstance(node.op, ast.Mult):
        return _flatten_mult(node.left) + _flatten_mult(node.right)
    if isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.UAdd):
        return _flatten_mult(node.operand)
    return [node]


def _product(factors, coeff=1):
    sign = 1
    if coeff < 0:
        sign = -1
        coeff = -coeff
    parts = []
    if coeff != 1:
        parts.append(ast.Constant(value=coeff))
    parts.extend(factors)
    if not parts:
        node = ast.Constant(value=1)
    else:
        node = parts[0]
        for part in parts[1:]:
            node = ast.BinOp(left=node, op=ast.Mult(), right=part)
    if sign < 0:
        node = ast.UnaryOp(op=ast.USub(), operand=node)
    return node


def _cancel_div(num, den):
    """Cancel integer coefficients, `pi`, and identical factors in a quotient."""
    n_coeff, d_coeff = 1, 1
    n_pi, d_pi = 0, 0
    n_rest, d_rest = [], []
    for factor in _flatten_mult(num):
        value = _const_int(factor)
        if value is not None:
            n_coeff *= value
        elif isinstance(factor, ast.Name) and factor.id == "pi":
            n_pi += 1
        else:
            n_rest.append(factor)
    for factor in _flatten_mult(den):
        value = _const_int(factor)
        if value is not None:
            d_coeff *= value
        elif isinstance(factor, ast.Name) and factor.id == "pi":
            d_pi += 1
        else:
            d_rest.append(factor)
    g = math.gcd(abs(n_coeff), abs(d_coeff)) or 1
    n_coeff //= g
    d_coeff //= g
    both = min(n_pi, d_pi)
    n_pi -= both
    d_pi -= both
    used = [False] * len(d_rest)
    n_keep = []
    for factor in n_rest:
        key = ast.dump(factor)
        hit = next((i for i, other in enumerate(d_rest)
                    if not used[i] and ast.dump(other) == key), None)
        if hit is None:
            n_keep.append(factor)
        else:
            used[hit] = True
    d_keep = [other for i, other in enumerate(d_rest) if not used[i]]
    n_keep = [ast.Name(id="pi", ctx=ast.Load())] * n_pi + n_keep
    d_keep = [ast.Name(id="pi", ctx=ast.Load())] * d_pi + d_keep
    left = _product(n_keep, n_coeff)
    if not d_keep and d_coeff == 1:
        return left
    return ast.BinOp(left=left, op=ast.Div(), right=_product(d_keep, d_coeff))


def _drop_ones(node):
    if isinstance(node, ast.BinOp):
        left, right = _drop_ones(node.left), _drop_ones(node.right)
        if isinstance(node.op, ast.Mult):
            if _is_one(left):
                return right
            if _is_one(right):
                return left
        elif isinstance(node.op, ast.Div) and _is_one(right):
            return left
        elif isinstance(node.op, ast.Pow):
            if _is_one(left):
                return ast.Constant(value=1)
            if _is_one(right):
                return left
        return ast.BinOp(left=left, op=node.op, right=right)
    if isinstance(node, ast.UnaryOp):
        operand = _drop_ones(node.operand)
        if isinstance(node.op, ast.UAdd):
            return operand
        return ast.UnaryOp(op=node.op, operand=operand)
    if isinstance(node, ast.Call):
        return ast.Call(func=node.func,
                        args=[_drop_ones(a) for a in node.args],
                        keywords=node.keywords)
    return node


def set_hbar_to_one(expr: str) -> str:
    """ħ = 1 when the formula writes unreduced h (ħ = h/2π).

    Substitutes `h → 2π` and cancels the resulting 2 and π factors:
    `h/(2*pi)*omega` → `omega`, `q*h/(4*pi*m)` → `q/(2*m)`.
    """
    sides = _split_equation(expr)
    if len(sides) == 2:
        return "=".join(set_hbar_to_one(side) for side in sides)
    tree = ast.parse(_sanitise(expr), mode="eval")

    def walk(node):
        if isinstance(node, ast.Expression):
            return walk(node.body)
        if isinstance(node, ast.Name) and node.id == "h":
            return ast.BinOp(left=ast.Constant(value=2), op=ast.Mult(),
                             right=ast.Name(id="pi", ctx=ast.Load()))
        if isinstance(node, ast.BinOp):
            left, right = walk(node.left), walk(node.right)
            if isinstance(node.op, ast.Div):
                return _cancel_div(left, right)
            return ast.BinOp(left=left, op=node.op, right=right)
        if isinstance(node, ast.UnaryOp):
            operand = walk(node.operand)
            if isinstance(node.op, ast.UAdd):
                return operand
            return ast.UnaryOp(op=node.op, operand=operand)
        if isinstance(node, ast.Call):
            return ast.Call(func=node.func,
                            args=[walk(a) for a in node.args],
                            keywords=node.keywords)
        return node

    return ast.unparse(_drop_ones(walk(tree))).replace(" ", "")


# ili_zipf.py: beta ~ Uniform(beta_min, beta_max) for the exponential model.
BETA_PRIOR_LO = 0.1
BETA_PRIOR_HI = 1.0


def _positive_rank_counts(counts):
    obs = np.sort(np.asarray(list(counts), dtype=float))[::-1]
    obs = obs[obs > 0]
    if obs.size < 2 or not np.isfinite(obs.sum()) or obs.sum() <= 0:
        raise ValueError("beta needs at least two positive ranks")
    return obs


def beta_posterior(counts, n_grid: int = 2001):
    """Posterior of β for the normalised exponential under a uniform prior.

    The prior is Uniform(0.1, 1), the range in the authors' `ili_zipf.py`.
    The likelihood is multinomial: rank r has probability proportional to
    e^{-β r}, over the ranks with a positive count. Returns
    (mean, std, 16th percentile, 84th percentile).
    """
    obs = _positive_rank_counts(counts)
    freq = obs / obs.sum()
    rank = np.arange(1, obs.size + 1)
    beta = np.linspace(BETA_PRIOR_LO, BETA_PRIOR_HI, n_grid)
    logp = -beta[:, None] * rank[None, :]
    logp -= np.logaddexp.reduce(logp, axis=1, keepdims=True)
    loglik = obs.sum() * (freq[None, :] * logp).sum(axis=1)
    loglik -= loglik.max()
    weight = np.exp(loglik)
    weight /= weight.sum()
    mean = float((weight * beta).sum())
    var = float((weight * (beta - mean) ** 2).sum())
    cdf = np.cumsum(weight)
    lo = float(beta[min(np.searchsorted(cdf, 0.16), beta.size - 1)])
    hi = float(beta[min(np.searchsorted(cdf, 0.84), beta.size - 1)])
    return mean, float(np.sqrt(var)), lo, hi


def fit_beta(counts) -> tuple[float, float]:
    """β as the posterior mean under Uniform(0.1, 1).

    Returns (mean, posterior standard deviation). The 16% and 84% points
    of the same posterior are `beta_posterior`.
    """
    mean, std, _, _ = beta_posterior(counts)
    return mean, std


def fit_beta_logls(counts) -> tuple[float, float]:
    """β by least squares of log frequency to the normalised exponential.

    Each observed rank is one residual of log(f) − log(p), so a rare rank
    pulls as hard as the most common rank. Returns (β, residual standard error).
    """
    from scipy.optimize import least_squares

    obs = _positive_rank_counts(counts)
    rank = np.arange(1, obs.size + 1, dtype=float)
    logf = np.log(obs / obs.sum())

    def resid(theta):
        logp = -theta[0] * rank
        logp = logp - np.logaddexp.reduce(logp)
        return logf - logp

    sol = least_squares(resid, np.array([0.3]))
    beta = float(sol.x[0])
    jac = np.asarray(sol.jac, dtype=float)
    dof = max(int(obs.size) - 1, 1)
    scale = float(np.sum(sol.fun ** 2) / dof)
    gram = float((jac.T @ jac)[0, 0])
    sigma = float(np.sqrt(scale / gram)) if gram > 0 else float("nan")
    return beta, sigma


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
    """Load any of the three released corpora. Only `inflation` needs repairs.

    Nullary counts are the CSV columns, not a guess from the formula
    (`p` in LFI is a variable; `A` in RCHI is a constant).
    """
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


def _nullary_rank(name: str) -> bool:
    return name in NULLARY or name.startswith("parameter_")


def arity_shares(entries: list[EIEntry]) -> dict[str, float]:
    """Fraction of operator instances that are nullary, unary and binary."""
    total = rank_counts(entries)
    nullary = sum(v for k, v in total.items() if _nullary_rank(k))
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


# --- Default estimator: multinomial maximum likelihood --------------------
# Each operator instance is one draw from p_r = exp(-beta r) / sum_s exp(-beta s)
# over the R observed ranks (counts sorted descending). This is the likelihood
# the shared-vs-separate tests use. curve_fit / log least squares stay available
# (fit_beta_logls, scripts/make_plots.fit_exp_curvefit) but are not the default.
MLE_BOUNDS = (0.0, 5.0)
BOOTSTRAP_REPS = 2000
BOOTSTRAP_SEED = 20261013


def _mle_from_sorted(obs: np.ndarray) -> float:
    from scipy.optimize import minimize_scalar
    rank = np.arange(1, obs.size + 1, dtype=float)
    n = obs.sum()
    rbar = float((obs * rank).sum() / n)

    def nll(beta):
        lp = -beta * rank
        return beta * rbar + np.logaddexp.reduce(lp)  # per-instance −loglik
    res = minimize_scalar(nll, bounds=MLE_BOUNDS, method="bounded",
                          options={"xatol": 1e-10})
    return float(res.x)


def fit_beta_mle(counts) -> float:
    """Multinomial MLE of β on the ranked positive counts (any order in)."""
    return _mle_from_sorted(_positive_rank_counts(counts))


def fisher_sigma(beta: float, n_ranks: int, n_instances: float) -> float:
    """1/sqrt(Fisher information): N Var_p(r) for the truncated exponential."""
    rank = np.arange(1, n_ranks + 1, dtype=float)
    lp = -beta * rank
    p = np.exp(lp - np.logaddexp.reduce(lp))
    var = float((p * rank ** 2).sum() - (p * rank).sum() ** 2)
    return float(1.0 / np.sqrt(n_instances * var)) if var > 0 else float("nan")


def _count_matrix(entries: list) -> np.ndarray:
    """Rows = expressions, columns = operator types (incl. x, a, c)."""
    per = []
    for e in entries:
        per.append(rank_counts([e]))
    keys = sorted({k for c in per for k in c})
    return np.array([[c.get(k, 0) for k in keys] for c in per], dtype=float)


def beta_with_errors(entries: list, reps: int = BOOTSTRAP_REPS,
                     seed: int = BOOTSTRAP_SEED) -> dict:
    """β (multinomial MLE) with a bootstrap over expressions and the Fisher error.

    Each bootstrap replicate resamples the expressions with replacement, sums
    their operator counts, re-ranks, and refits. Returns n, beta, boot_sigma,
    boot_p16, boot_p84, fisher_sigma, n_instances, n_ranks.
    """
    mat = _count_matrix(entries)
    tot = mat.sum(axis=0)
    obs = np.sort(tot[tot > 0])[::-1]
    beta = _mle_from_sorted(obs)
    out = {"n": len(entries), "beta": beta, "n_instances": int(obs.sum()),
           "n_ranks": int(obs.size),
           "fisher_sigma": fisher_sigma(beta, obs.size, obs.sum())}
    if reps and len(entries) > 1:
        rng = np.random.default_rng(seed)
        bs = np.empty(reps)
        n = len(entries)
        for i in range(reps):
            w = np.bincount(rng.integers(0, n, n), minlength=n)
            t = w @ mat
            o = np.sort(t[t > 0])[::-1]
            bs[i] = _mle_from_sorted(o) if o.size >= 2 else np.nan
        bs = bs[np.isfinite(bs)]
        out.update(boot_sigma=float(bs.std(ddof=1)),
                   boot_p16=float(np.percentile(bs, 16)),
                   boot_p84=float(np.percentile(bs, 84)), reps=int(bs.size), seed=seed)
    return out


def counts_with_fisher(rank_ordered_counts) -> dict:
    """For aggregate-only data (e.g. published rank counts): MLE and Fisher error only."""
    obs = _positive_rank_counts(rank_ordered_counts)
    beta = _mle_from_sorted(obs)
    return {"beta": beta, "fisher_sigma": fisher_sigma(beta, obs.size, obs.sum()),
            "n_instances": int(obs.sum()), "n_ranks": int(obs.size)}
