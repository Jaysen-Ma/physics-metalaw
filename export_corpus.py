"""Export the three released corpora to one self-describing JSON file.

The CSVs the meta-law authors publish are spreadsheet exports: a formula as a
flat string, plus three hand-maintained nullary tallies. The object actually
being studied is the *expression tree*, which a CSV cell cannot hold, so every
consumer re-parses the string and re-derives the tree. This script does that
parse once, under `metalaw`'s validated convention, and writes the result --
tree included, and with each node tagged by the rank bucket it contributes to.

Provenance is kept alongside: the authors' declared complexity, our recomputed
complexity, and a flag for every formula where the two disagree or where the
source text had to be repaired to parse. Nothing is silently corrected.

Usage:  .venv/bin/python export_corpus.py [-o runs/corpora.json]
The corpus clone in .scratch/ is CC0 (Constantin, Bartlett, Desmond & Ferreira
2024, arXiv:2408.11065); it is re-cloned rather than vendored, so this script's
output is generated, not tracked.
"""
from __future__ import annotations

import argparse
import ast
import csv
import json
from pathlib import Path

import metalaw as M

DATA = Path(__file__).parent / ".scratch/physics-formulae-patterns/data"
CORPORA = [("inflation", "Inflation.csv"), ("feynman", "Feynman.csv"),
           ("wiki", "Wiki.csv")]

BUCKET = {"variable": "nullary", "number": "nullary", "constant": "nullary"}


def to_tree(expr: str) -> dict:
    """Parse to the same tree `metalaw.count_operators` walks.

    Each node carries the rank key it contributes to, so the tree and the
    operator counts can never drift apart.
    """
    def walk(node):
        if isinstance(node, ast.Expression):
            return walk(node.body)
        if isinstance(node, ast.BinOp):
            if (isinstance(node.op, ast.Pow)
                    and isinstance(node.right, ast.Constant)
                    and isinstance(node.right.value, int)
                    and node.right.value in M.COLLAPSED_POWERS):
                # the literal exponent collapses into the operator and is not
                # counted as a node of its own
                return {"key": f"Pow{node.right.value}", "arity": "unary",
                        "children": [walk(node.left)]}
            return {"key": M.BINARY[type(node.op)], "arity": "binary",
                    "children": [walk(node.left), walk(node.right)]}
        if isinstance(node, ast.UnaryOp):
            return {"key": "Neg" if isinstance(node.op, ast.USub) else "Pos",
                    "arity": "unary", "children": [walk(node.operand)]}
        if isinstance(node, ast.Call):
            return {"key": node.func.id, "arity": "unary",
                    "children": [walk(a) for a in node.args]}
        if isinstance(node, ast.Constant):
            return {"key": "number", "arity": "nullary", "text": str(node.value)}
        if isinstance(node, ast.Name):
            return {"key": "symbol", "arity": "nullary", "text": node.id}
        raise TypeError(f"unhandled node: {ast.dump(node)}")

    return walk(ast.parse(M._sanitise(expr), mode="eval"))


def symbol_roles(entry: M.EIEntry, symbols: set[str]) -> dict[str, str]:
    """The authors publish only *counts* of variables and constants, never which
    symbol is which, and the split is a per-formula physics judgement (`c` is a
    constant in I.15.3x and a variable elsewhere). So it is only recoverable at
    the ends of the range; anything else is honestly marked unknown.
    """
    if entry.n_constants == 0:
        return {s: "variable" for s in symbols}
    if entry.n_variables == 0:
        return {s: "constant" for s in symbols}
    return {s: "unknown" for s in symbols}


def build(corpus: str, fname: str) -> list[dict]:
    raw = {}
    ci = M.CORPUS_COLUMNS[corpus]
    rows = [r for r in csv.reader(open(DATA / fname), delimiter=";")
            if any(c.strip() for c in r)][1:]
    for i, r in enumerate(rows):
        if r[ci[1]].strip():
            raw.setdefault(r[0].strip() or "(unnamed)", (i + 2, r[ci[1]]))

    out = []
    for e in M.load_corpus(DATA / fname, corpus):
        line, original = raw.get(e.name, (None, e.formula))
        counts = M.count_operators(e.formula)
        symbols = {k[5:] for k in counts if k.startswith("name:")}
        flags = []
        if e.repaired:
            flags.append("bracket_repaired")
        if e.complexity_declared != e.complexity_recomputed:
            flags.append("complexity_mismatch")
        n_lit = counts.get(M.NULLARY_NUMBER, 0)
        if n_lit != e.n_numbers:
            flags.append("number_count_mismatch")
        out.append({
            "id": f"{corpus}.{e.name}",
            "corpus": corpus,
            "name": e.name,
            "expr": e.formula,
            "expr_original": original if e.repaired else None,
            "tree": to_tree(e.formula),
            "operators": dict(sorted(e.operator_counts.items())),
            "nullary": {"variables": e.n_variables, "numbers": e.n_numbers,
                        "constants": e.n_constants},
            "symbol_roles": symbol_roles(e, symbols),
            "complexity": {"declared": e.complexity_declared,
                           "recomputed": e.complexity_recomputed},
            "flags": flags,
            "source": {"file": fname, "line": line},
        })
    return out


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("-o", "--out", default="runs/corpora.json")
    args = p.parse_args()

    payload = {
        "convention": (
            "Complexity is total node count. Nullary nodes rank as three "
            "separate types (variable/number/constant). An integer exponent of "
            "2, 3 or 4 collapses into a single unary Pow2/Pow3/Pow4 and the "
            "literal is not counted. Unary minus is its own operator."),
        "source": "github.com/andrei-const/physics-formulae-patterns (CC0); "
                  "Constantin, Bartlett, Desmond & Ferreira 2024, arXiv:2408.11065",
        "corpora": {},
    }
    for corpus, fname in CORPORA:
        formulas = build(corpus, fname)
        ranks = sorted(M.rank_counts(M.load_corpus(DATA / fname, corpus)).items(),
                       key=lambda kv: -kv[1])
        beta, sigma = M.fit_beta([v for _, v in ranks])
        published = DATA / f"{corpus}_data.txt"
        payload["corpora"][corpus] = {
            "formulas": formulas,
            "ranks": [{"key": k, "count": v} for k, v in ranks],
            "ranks_published": ([int(x) for x in open(published)]
                                if published.exists() else None),
            "beta": beta, "beta_sigma": sigma,
        }

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, indent=1))
    total = sum(len(c["formulas"]) for c in payload["corpora"].values())
    print(f"wrote {out}  ({total} formulas, {out.stat().st_size/1024:.0f} kB)")
    for name, c in payload["corpora"].items():
        flagged = sum(1 for f in c["formulas"] if f["flags"])
        print(f"  {name:10s} n={len(c['formulas']):3d}  beta={c['beta']:.3f}  "
              f"flagged={flagged}")


if __name__ == "__main__":
    main()
