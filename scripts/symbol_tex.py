"""LaTeX for the identifiers produced by the LaTeX-to-sympy extraction (x' -> xp, <c> -> avgc, ...).

Used by scripts/build_feynman_lectures_catalog.py to render the `latex` field of each record.
"""
from __future__ import annotations

import re

# The extraction encodes LaTeX decorations in identifiers (x' -> xp, <c> -> avgc,
# \bar{x} -> xbar, \varkappa, multi-letter subscripts, zeta(3/2) -> zeta3o2).
# metalaw's default renderer does not know these, so the page showed them as
# italic letter strings (e.g. "varkappa", "lHopitalsrule"). This maps them back.
GREEK = ("alpha beta gamma delta epsilon varepsilon zeta eta theta vartheta iota kappa varkappa lambda mu nu xi "
         "pi varpi rho varrho sigma varsigma tau upsilon phi varphi chi psi omega").split()
GREEK_CAP = "Gamma Delta Theta Lambda Xi Pi Sigma Upsilon Phi Psi Omega".split()
SUB_TEX = {"perp": r"\perp", "parallel": r"\parallel", "odot": r"\odot", "infty": r"\infty", "sun": r"\odot"}
NAME_TEX = {
    "lambd": r"\lambda", "zeta3o2": r"\zeta(3/2)", "zeta3": r"\zeta(3)", "zeta5": r"\zeta(5)",
    "avgxixj": r"\langle x_i x_j \rangle", "avgxi": r"\langle x_i \rangle", "avgxj": r"\langle x_j \rangle",
    "avgc2": r"\langle c^2 \rangle", "avgr2": r"\langle r^2 \rangle", "avgExEy": r"\langle E_x E_y \rangle",
    "x_n___1": r"x_{n+1}", "Z_i___1": r"Z_{i+1}", "g_i___1": r"g_{i+1}",
    "lambda_i__j__k_l__": r"\lambda_{ijkl}", "R_alpha__beta__delta_gamma__": r"R_{\alpha\beta\gamma\delta}",
    "G_h__k_l_": r"G_{hkl}", "E_l__c_o_": r"E_{\mathrm{loc}}", "u_r__a_d_": r"u_{\mathrm{rad}}",
    "u_m__a_g_": r"u_{\mathrm{mag}}", "m_b__l_o_": r"m_{\mathrm{bol}}", "M_b__l_o_": r"M_{\mathrm{bol}}",
    "n_0_j1": r"n_{0,j+1}", "g_0_j1": r"g_{0,j+1}", "BV_0": r"(B-V)_0", "Gamma_12_tau": r"\Gamma_{12}(\tau)",
    "Gamma_11_0": r"\Gamma_{11}(0)", "Gamma_22_0": r"\Gamma_{22}(0)", "gamma_12_tau": r"\gamma_{12}(\tau)",
    "Y_1_m1": r"Y_1^{-1}", "Y_1_1": r"Y_1^{1}", "Y_2_m1": r"Y_2^{-1}", "Y_2_1": r"Y_2^{1}",
    "a_star": r"a^{*}", "b_star": r"b^{*}", "c_star": r"c^{*}", "hbar": r"\hbar", "chi_HP": r"\chi_{\mathrm{HP}}",
    "Omega_m0": r"\Omega_{\mathrm{m}0}", "Omega_Lambda0": r"\Omega_{\Lambda 0}", "omega_f": r"\omega_{\mathrm{f}}",
    "nu_ch": r"\nu_{\mathrm{ch}}",
}


def _atom_tex(part: str) -> str:
    if part in GREEK or part in GREEK_CAP:
        return "\\" + part
    m = re.fullmatch(r"avg(.+)", part)
    if m:
        return r"\langle " + _atom_tex(m.group(1)) + r" \rangle"
    m = re.fullmatch(r"(.+?)bar", part)
    if m:
        return r"\bar{" + _atom_tex(m.group(1)) + "}"
    m = re.fullmatch(r"([A-Za-z]|" + "|".join(GREEK + GREEK_CAP) + r")p", part)
    if m:
        return _atom_tex(m.group(1)) + "'"
    m = re.fullmatch(r"([A-Za-z]|" + "|".join(GREEK + GREEK_CAP) + r")(\d+)", part)
    if m:
        return _atom_tex(m.group(1)) + "_{" + m.group(2) + "}"
    if len(part) > 1 and part.isalpha():
        return r"\mathrm{" + part + "}"
    return part


def symbol_tex(name: str):
    """LaTeX for an extracted identifier, or None to use metalaw's default."""
    if name in NAME_TEX:
        return NAME_TEX[name]
    base, *subs = name.split("_")
    if "" in subs:
        return None
    sub = ",".join(SUB_TEX.get(b) or (_atom_tex(b) if (b in GREEK or b in GREEK_CAP) else (b if len(b) == 1 or b.isdigit() else r"\mathrm{" + b + "}")) for b in subs)
    tex = _atom_tex(base)
    return tex + ("_{" + sub + "}" if subs else "")
