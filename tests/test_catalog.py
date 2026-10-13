"""Sanity checks for the working Feynman catalog."""
from __future__ import annotations

import unittest

from physics_metalaw import catalog as C
from physics_metalaw import metalaw as M


class FeynmanCatalogTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.recs = C.seed_feynman()
        cls.by_orig = {r.original_id: r for r in cls.recs}

    def test_released_keeps_all_hundred(self):
        self.assertEqual(len(self.recs), 103)
        self.assertEqual(len(C.select(self.recs, "released")), 100)

    def test_primary_uses_gaussian_generalization(self):
        primary = C.select(self.recs, "primary")
        ids = {r.id for r in primary}
        self.assertEqual(len(primary), 97)
        self.assertIn("I.6.20b", ids)
        self.assertIn("I.9.18b", ids)
        self.assertIn("I.12.2b", ids)
        self.assertIn("I.12.4b", ids)
        self.assertIn("II.6.15b", ids)
        self.assertNotIn("I.6.20a", ids)
        self.assertNotIn("I.6.20", ids)
        self.assertNotIn("I.9.18", ids)
        self.assertNotIn("I.12.2", ids)
        self.assertNotIn("I.12.4", ids)
        self.assertNotIn("II.6.15a", ids)

    def test_lecture_policy_is_the_original_form(self):
        lecture = C.select(self.recs, "lecture")
        self.assertTrue(all(r.status == "original" for r in lecture))
        self.assertEqual({r.status for r in self.recs}, {"original", "rewritten"})
        self.assertIn("I.6.20", {r.id for r in lecture})
        self.assertIn("I.39.10", {r.id for r in lecture})
        self.assertEqual(len(lecture), 100)
        self.assertIn("I.6.20a", {r.id for r in lecture})
        self.assertIn("I.6.20b", {r.id for r in lecture})
        self.assertIn("II.6.15b", {r.id for r in lecture})
        self.assertEqual(
            {r.id for r in self.recs if r.status == "rewritten"},
            {"I.9.18b", "I.12.2b", "I.12.4b"},
        )
        self.assertEqual(self.by_orig["I.6.2"].status, "original")
        self.assertEqual(self.by_orig["I.39.1"].status, "original")

    def test_gaussian_family_urls_and_roles(self):
        a = self.by_orig["I.6.2a"]
        g = self.by_orig["I.6.2"]
        b = self.by_orig["I.6.2b"]
        chapter = "https://www.feynmanlectures.caltech.edu/I_06.html"
        self.assertEqual(a.id, "I.6.20a")
        self.assertEqual(g.id, "I.6.20")
        self.assertEqual(b.id, "I.6.20b")
        self.assertEqual(a.status, "original")
        self.assertEqual(g.status, "original")
        self.assertEqual(b.status, "original")
        self.assertEqual(g.family_role, "centred")
        self.assertEqual(b.family_role, "generalization")
        self.assertFalse(g.representative)
        self.assertTrue(b.representative)
        self.assertFalse(a.representative)
        for rec in (a, g, b):
            self.assertEqual(rec.family, "I.6.20")
            self.assertEqual(rec.source_kind, "html")
            self.assertEqual(rec.source_url, chapter)
        self.assertEqual(g.source_locator, "(6.20)")

    def test_counts_still_match_released_complexity(self):
        g = self.by_orig["I.6.2"]
        self.assertEqual(g.complexity, g.original_complexity)
        self.assertEqual(g.complexity, 15)
        # rank_counts duck-types FormulaRecord the same way as EIEntry
        n = sum(M.rank_counts([g]).values())
        self.assertEqual(n, 15)

    def test_tree_nodes_match_complexity(self):
        g = self.by_orig["I.6.2"]
        tree = M.expression_tree(g.formula)
        self.assertEqual(M.tree_size(tree), g.complexity)
        self.assertEqual(tree["label"], "Div")
        for rec in self.recs:
            size = M.tree_size(M.expression_tree(rec.formula))
            self.assertEqual(size, M.complexity(rec.formula), rec.original_id)

    def test_json_carries_the_tree(self):
        rec = self.by_orig["I.6.2a"]
        blob = rec.to_json()
        self.assertEqual(blob["tree"]["label"], "Div")
        self.assertEqual(M.tree_size(blob["tree"]), 11)
        self.assertIn(r"\theta", blob["latex"])

    def test_I918_G_is_a_constant(self):
        rec = self.by_orig["I.9.18"]
        self.assertEqual(rec.n_variables, 8)
        self.assertEqual(rec.n_constants, 1)
        self.assertFalse(rec.representative)
        self.assertEqual(rec.family, "I.9.18")
        self.assertIn("G", rec.constant_names)
        g = _find_label(rec.to_json()["tree"], "G")
        self.assertEqual(g["kind"], "constant")

    def test_I918_family_defaults_to_r_squared(self):
        b = self.by_orig["I.9.18b"]
        self.assertEqual(b.formula, "G*m1*m2/r**2")
        self.assertEqual(b.status, "rewritten")
        self.assertEqual(b.family_role, "convention")
        self.assertTrue(b.representative)
        self.assertFalse(b.in_released)
        self.assertEqual(b.complexity, 8)
        self.assertEqual(b.n_variables, 3)
        self.assertEqual(b.n_constants, 1)
        g = _find_label(b.to_json()["tree"], "G")
        self.assertEqual(g["kind"], "constant")
        self.assertIn("r", M.formula_latex(b.formula))

    def test_I107_c_is_a_constant(self):
        rec = self.by_orig["I.10.7"]
        self.assertEqual(rec.n_variables, 2)
        self.assertEqual(rec.n_constants, 1)
        self.assertIn("c", rec.constant_names)
        c = _find_label(rec.to_json()["tree"], "c")
        self.assertEqual(c["kind"], "constant")

    def test_II615_family_defaults_to_polar(self):
        a = self.by_orig["II.6.15a"]
        b = self.by_orig["II.6.15b"]
        self.assertEqual(a.family, "II.6.15")
        self.assertEqual(b.family, "II.6.15")
        self.assertEqual(a.status, "original")
        self.assertEqual(b.status, "original")
        self.assertFalse(a.representative)
        self.assertTrue(b.representative)
        self.assertEqual(a.family_role, "cartesian")
        self.assertEqual(b.family_role, "generalization")

    def test_I122_family_defaults_to_inverse_square(self):
        raw = self.by_orig["I.12.2"]
        b = self.by_orig["I.12.2b"]
        self.assertEqual(raw.formula, "q1*q2*r/(4*pi*epsilon*r**3)")
        self.assertFalse(raw.representative)
        self.assertEqual(raw.family_role, "vector")
        self.assertEqual(b.formula, "q1*q2/(4*pi*epsilon*r**2)")
        self.assertEqual(b.status, "rewritten")
        self.assertEqual(b.family_role, "convention")
        self.assertTrue(b.representative)
        self.assertFalse(b.in_released)
        self.assertEqual(b.complexity, 12)
        self.assertEqual(b.n_variables, 3)
        eps = _find_label(b.to_json()["tree"], "epsilon")
        self.assertEqual(eps["kind"], "constant")
        self.assertIn("r^{2}", M.formula_latex(b.formula))
        self.assertNotIn("r}^{3", M.formula_latex(b.formula))

    def test_distance_latex_has_brackets(self):
        dist = M.formula_latex("sqrt((x2-x1)**2+(y2-y1)**2)")
        grav = M.formula_latex("G*m1*m2/((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)")
        self.assertIn(r"\left(x_{2} - x_{1}\right)^{2}", dist)
        self.assertIn(r"\left(y_{2} - y_{1}\right)^{2}", dist)
        self.assertIn(r"\left(x_{2} - x_{1}\right)^{2}", grav)
        self.assertIn(r"\left(z_{2} - z_{1}\right)^{2}", grav)

    def test_I124_family_defaults_to_inverse_square(self):
        raw = self.by_orig["I.12.4"]
        b = self.by_orig["I.12.4b"]
        self.assertEqual(raw.formula, "q1*r/(4*pi*epsilon*r**3)")
        self.assertFalse(raw.representative)
        self.assertEqual(raw.family_role, "vector")
        self.assertEqual(b.formula, "q1/(4*pi*epsilon*r**2)")
        self.assertEqual(b.status, "rewritten")
        self.assertEqual(b.family_role, "convention")
        self.assertTrue(b.representative)
        self.assertFalse(b.in_released)
        self.assertEqual(b.complexity, 10)
        self.assertEqual(b.n_variables, 2)
        eps = _find_label(b.to_json()["tree"], "epsilon")
        self.assertEqual(eps["kind"], "constant")

    def test_Ef_latex_is_E_f(self):
        for orig in ("I.12.5", "I.12.11"):
            rec = self.by_orig[orig]
            tex = M.formula_latex(rec.formula)
            self.assertIn(r"E_{f}", tex, orig)
            self.assertNotIn("Ef", tex, orig)

    def test_I1312_G_is_a_constant(self):
        rec = self.by_orig["I.13.12"]
        self.assertEqual(rec.n_variables, 4)
        self.assertEqual(rec.n_constants, 1)
        self.assertIn("G", rec.constant_names)
        g = _find_label(rec.to_json()["tree"], "G")
        self.assertEqual(g["kind"], "constant")

    def test_I143_g_is_a_constant(self):
        rec = self.by_orig["I.14.3"]
        self.assertEqual(rec.n_variables, 2)
        self.assertEqual(rec.n_constants, 1)
        self.assertIn("g", rec.constant_names)
        g = _find_label(rec.to_json()["tree"], "g")
        self.assertEqual(g["kind"], "constant")

    def test_physics_constants_are_not_variables(self):
        lande_g = _find_label(self.by_orig["II.34.11"].to_json()["tree"], "g_")
        self.assertEqual(lande_g["kind"], "variable")
        rec_c = self.by_orig["I.29.4"]
        self.assertEqual(rec_c.original_n_variables, 2)
        self.assertEqual(rec_c.n_variables, 1)
        self.assertEqual(rec_c.n_constants, 1)
        self.assertEqual(rec_c.complexity, rec_c.original_complexity)
        rec_h = self.by_orig["I.34.27"]
        self.assertEqual(rec_h.original_n_variables, 2)
        self.assertEqual(rec_h.n_variables, 1)
        self.assertEqual(rec_h.n_constants, 2)
        self.assertIn("h", rec_h.constant_names)
        for rec in self.recs:
            tree = rec.to_json()["tree"]
            for name in C.PHYSICS_CONSTANTS:
                node = _find_label(tree, name)
                if node is not None:
                    self.assertEqual(node["kind"], "constant",
                                     f"{rec.original_id} {name}")
            self.assertEqual(
                rec.n_variables, _count_kind(tree, "variable"), rec.original_id)
            self.assertEqual(
                rec.n_numbers, _count_kind(tree, "number"), rec.original_id)
            self.assertEqual(
                rec.n_constants, _count_kind(tree, "constant"), rec.original_id)
            self.assertEqual(rec.complexity, M.tree_size(tree), rec.original_id)

    def test_I391_is_I3910(self):
        rec = self.by_orig["I.39.1"]
        self.assertEqual(rec.id, "I.39.10")
        self.assertEqual(rec.original_id, "I.39.1")
        self.assertEqual(rec.status, "original")
        self.assertEqual(rec.formula, "3/2*P_F*V")
        self.assertEqual(rec.original_formula, "3/2*pr*V")
        self.assertEqual(rec.source_locator, "(39.10)")
        self.assertEqual(rec.complexity, 7)
        self.assertIn(r"P_{F}", rec.to_json()["latex"])

    def test_symbol_tex_aliases(self):
        r10 = self.by_orig["I.39.1"].to_json()["latex"]
        r11 = self.by_orig["I.39.11"].to_json()["latex"]
        r22 = self.by_orig["I.39.22"]
        self.assertEqual(r10, r"\frac{3}{2} P_{F} V")
        self.assertIn(r"P_{F}", r11)
        self.assertEqual(r22.output, "P_F")
        self.assertIn(r"k_{b}", r22.to_json()["latex"])
        self.assertNotIn("kb", r22.to_json()["latex"])
        drift = self.by_orig["I.43.16"].to_json()["latex"]
        self.assertIn(r"\mu_{\text{drift}}", drift)
        self.assertIn(r"V_{e}", drift)
        mob = self.by_orig["I.43.31"].to_json()["latex"]
        self.assertIn(r"\mu", mob)
        self.assertIn(r"k_{b}", mob)
        self.assertNotIn("mob", mob)
        sound = self.by_orig["I.47.23"].to_json()["latex"]
        self.assertIn(r"P", sound)
        self.assertNotIn("pr", sound)
        for rec in self.recs:
            tex = rec.to_json()["latex"]
            if "kb" in C._name_occurrences(rec.formula):
                self.assertIn(r"k_{b}", tex, rec.id)
                self.assertNotIn("kb", tex, rec.id)
        mom = self.by_orig["II.15.4"].to_json()["latex"]
        self.assertIn(r"\mu", mom)
        self.assertNotIn("mom", mom)
        vec = self.by_orig["III.10.19"].to_json()["latex"]
        self.assertIn(r"B_{x}", vec)
        self.assertIn(r"B_{y}", vec)
        self.assertIn(r"B_{z}", vec)
        self.assertIn(r"n_{\rho}", self.by_orig["II.11.20"].to_json()["latex"])
        self.assertIn(r"\vec{A}", self.by_orig["III.21.20"].to_json()["latex"])
        self.assertEqual(M.formula_latex("g_*q*B/(2*m)"),
                         r"\frac{g q B}{2 m}")

    def test_I3217_numbers_drop_collapsed_four(self):
        rec = self.by_orig["I.32.17"]
        self.assertEqual(rec.original_n_numbers, 5)
        self.assertEqual(rec.n_numbers, 4)
        self.assertEqual(rec.n_constants, 3)
        self.assertEqual(rec.n_variables, 5)
        self.assertEqual(rec.complexity, 29)
        tree = rec.to_json()["tree"]
        self.assertIsNotNone(_find_label(tree, "Pow4"))
        self.assertEqual(_find_label(tree, "epsilon")["kind"], "constant")
        self.assertEqual(_find_label(tree, "c")["kind"], "constant")
        self.assertEqual(_find_label(tree, "pi")["kind"], "constant")

    def test_I305_lambda_has_backslash(self):
        rec = self.by_orig["I.30.5"]
        self.assertEqual(rec.formula, "arcsin(lambd/(n*d))")
        tex = rec.to_json()["latex"]
        self.assertEqual(tex, r"\arcsin\left(\frac{\lambda}{n d}\right)")

    def test_I303_intensity_is_I(self):
        rec = self.by_orig["I.30.3"]
        self.assertEqual(rec.original_formula,
                         "Int_0*sin(n*theta/2)**2/sin(theta/2)**2")
        self.assertEqual(rec.formula,
                         "I_0*sin(n*theta/2)**2/sin(theta/2)**2")
        self.assertEqual(rec.output, "I")
        self.assertEqual(rec.complexity, 15)
        tex = rec.to_json()["latex"]
        self.assertIn(r"I_{0}", tex)
        self.assertNotIn("Int", tex)

    def test_I246_half_sticks_with_x_squared(self):
        rec = self.by_orig["I.24.6"]
        self.assertEqual(rec.original_formula,
                         "1/2*m*(omega**2+omega_0**2)*1/2*x**2")
        self.assertEqual(rec.formula,
                         "1/2*m*(omega**2+omega_0**2)*(1/2*x**2)")
        self.assertEqual(rec.complexity, 18)
        self.assertEqual(
            rec.to_json()["latex"],
            r"\frac{1}{2} m \left(\omega^{2} + \omega_{0}^{2}\right) "
            r"\frac{1}{2} x^{2}",
        )

    def test_c1_drops_c_and_lowers_complexity(self):
        self.assertEqual(M.set_names_to_one("omega/c"), "omega")
        self.assertEqual(M.set_names_to_one("m_0/sqrt(1-v**2/c**2)"),
                         "m_0/sqrt(1-v**2)")
        self.assertEqual(M.set_names_to_one("m*c**2/sqrt(1-v**2/c**2)"),
                         "m/sqrt(1-v**2)")
        mass = self.by_orig["I.10.7"].to_json()
        self.assertEqual(mass["c1"]["formula"], "m_0/sqrt(1-v**2)")
        self.assertEqual(mass["c1"]["complexity"], 7)
        self.assertLess(mass["c1"]["complexity"], mass["complexity"])
        self.assertIsNone(_find_label(mass["c1"]["tree"], "c"))
        self.assertNotIn("c", mass["c1"]["constant_names"])
        k = self.by_orig["I.29.4"].to_json()["c1"]
        self.assertEqual(k["formula"], "omega")
        self.assertEqual(k["complexity"], 1)
        self.assertIsNone(self.by_orig["I.8.14"].to_json().get("c1"))
        self.assertIsNone(self.by_orig["I.47.23"].to_json().get("c1"))
        for rec in self.recs:
            blob = rec.to_json()
            if "c" not in C._name_occurrences(rec.formula):
                self.assertNotIn("c1", blob, rec.original_id)
                continue
            c1 = blob["c1"]
            self.assertLess(c1["complexity"], rec.complexity, rec.original_id)
            self.assertEqual(c1["complexity"], M.tree_size(c1["tree"]),
                             rec.original_id)
            self.assertIsNone(_find_label(c1["tree"], "c"), rec.original_id)

    def test_natural_units_drop_painted_physics_constants_only(self):
        newton = self.by_orig["I.9.18b"]
        self.assertEqual(C.rewrite_units(newton, "natural"), "m1*m2/r**2")
        gauss = self.by_orig["I.6.2"]
        self.assertEqual(C.rewrite_units(gauss, "natural"), gauss.formula)
        si = [r for r in C.seed_inflation() if r.id == "SI"][0]
        nat = C.rewrite_units(si, "natural")
        self.assertIn("M_Pl", C._name_occurrences(nat))
        self.assertEqual(nat, si.formula)

    def test_hbar1_drops_unreduced_h(self):
        self.assertEqual(M.set_hbar_to_one("(h/(2*pi))*omega"), "omega")
        self.assertEqual(M.set_hbar_to_one("q*h/(4*pi*m)"), "q/(2*m)")
        e = self.by_orig["I.34.27"].to_json()["hbar1"]
        self.assertEqual(e["formula"], "omega")
        self.assertEqual(e["complexity"], 1)
        self.assertIsNone(_find_label(e["tree"], "h"))
        bohr = self.by_orig["I.38.12"].to_json()["hbar1"]
        self.assertEqual(bohr["formula"], "4*pi*epsilon/(m*q**2)")
        magneton = self.by_orig["II.34.29a"].to_json()["hbar1"]
        self.assertEqual(magneton["formula"], "q/(2*m)")
        planck = self.by_orig["I.41.16"].to_json()
        self.assertIn("c1_hbar1", planck)
        self.assertEqual(
            planck["c1_hbar1"]["formula"],
            "omega**3/(pi**2*(exp(omega/(kb*T))-1))",
        )
        self.assertIsNone(_find_label(planck["c1_hbar1"]["tree"], "h"))
        self.assertIsNone(_find_label(planck["c1_hbar1"]["tree"], "c"))
        for rec in self.recs:
            blob = rec.to_json()
            if "h" not in rec.constant_names:
                self.assertNotIn("hbar1", blob, rec.original_id)
                continue
            hbar = blob["hbar1"]
            self.assertLess(hbar["complexity"], rec.complexity, rec.original_id)
            self.assertIsNone(_find_label(hbar["tree"], "h"), rec.original_id)

    def test_generalization_latex(self):
        b = self.by_orig["I.6.2b"]
        self.assertEqual(b.complexity, 17)
        tex = M.formula_latex(b.formula)
        self.assertIn(r"\theta_{1}", tex)
        self.assertIn(r"\sigma", tex)
        self.assertIn(r"\left(\frac", tex)


class ReleasedCsvTests(unittest.TestCase):
    def test_clone_head_is_stamped_into_catalog_meta(self):
        sha = C.released_clone_commit()
        self.assertIsNotNone(sha)
        self.assertRegex(sha, r"^[0-9a-f]{40}$")
        meta = C.released_csv_meta("Feynman.csv")
        self.assertEqual(meta["released_commit"], sha)
        self.assertEqual(meta["released_repo"], C.RELEASED_REPO_URL)
        self.assertEqual(
            meta["released_csv"], "physics-formulae-patterns/data/Feynman.csv")


class WikiCatalogTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.recs = C.seed_wiki()
        cls.by_id = {r.id: r for r in cls.recs}

    def test_released_keeps_all_forty_one(self):
        released = C.select(self.recs, "released")
        extension = [r for r in self.recs if r.status == "extension"]
        self.assertEqual(len(released), 41)
        self.assertEqual(len(extension), len(C._WIKI_EXTENSION))
        qm = [r for r in self.recs if r.status == "qm"]
        self.assertEqual(len(qm), len(C._WIKI_QM))
        self.assertTrue(all(not r.in_released for r in qm))
        self.assertEqual(len(C.select(self.recs, "primary")), len(self.recs))
        self.assertTrue(all(r.status == "original" and r.in_released for r in released))
        self.assertTrue(all(not r.in_released for r in extension))
        self.assertTrue(all(r.corpus == "wiki" for r in self.recs))
        self.assertIn("R", C.PHYSICS_CONSTANTS)

    def test_source_is_the_wikipedia_article(self):
        rec = self.by_id["Wiki_1"]
        self.assertEqual(rec.formula, "sqrt((K + (4/3)*mu)/rho)")
        self.assertEqual(rec.name, "compressional-wave speed v_p")
        self.assertIn("Adams", rec.source_locator)
        self.assertTrue(rec.source_url.endswith("Adams%E2%80%93Williamson_equation")
                        or "Adams" in rec.source_url)
        self.assertNotEqual(rec.source_url, C.WIKI_LIST_URL)
        self.assertIn("Potts", self.by_id["Wiki_11"].source_url)
        self.assertIn("Doppler", self.by_id["Wiki_41"].source_url)
        for rec in self.recs:
            self.assertTrue(rec.name, rec.id)
            self.assertIn("/wiki/", rec.source_url, rec.id)
            self.assertNotEqual(rec.source_url, C.WIKI_LIST_URL, rec.id)
        released = C.select(self.recs, "released")
        self.assertEqual(set(C._WIKI_PAGES), {r.id for r in released})
        self.assertEqual(
            len({C._WIKI_PAGES[r.id][0] for r in released}), 21)

    def test_every_formula_parses(self):
        for rec in self.recs:
            tree = rec.to_json()["tree"]
            self.assertEqual(rec.complexity, M.tree_size(tree), rec.id)
            self.assertNotIn("_sym", rec.to_json()["latex"], rec.id)

    def test_schrodinger_equation_does_not_count_equality(self):
        rec = next(r for r in self.recs if r.name.startswith("time-dependent"))
        tree = rec.to_json()["tree"]
        self.assertEqual(tree["kind"], "relation")
        self.assertEqual(tree["label"], "=")
        self.assertNotIn("=", rec.operator_counts)
        self.assertIn("Partial2", rec.operator_counts)
        self.assertIn("Partial", rec.operator_counts)
        self.assertEqual(_find_label(tree, "hbar")["kind"], "constant")
        self.assertEqual(_find_label(tree, "I")["kind"], "constant")
        self.assertIn(r"\hbar", rec.to_json()["latex"])
        self.assertIn(r"\partial^{2}", rec.to_json()["latex"])
        self.assertEqual(
            M.tree_size(tree),
            sum(M.tree_size(c) for c in tree["children"]))

    def test_Wiki33_lambda_is_wavelength(self):
        rec = self.by_id["Wiki_33"]
        self.assertIn("lambda", rec.formula)
        tex = rec.to_json()["latex"]
        self.assertIn(r"\lambda", tex)
        self.assertNotIn("lambda_sym", tex)
        leaf = _find_label(rec.to_json()["tree"], "lambda")
        self.assertIsNotNone(leaf)
        self.assertEqual(leaf["kind"], "variable")

    def test_Wiki7_R_is_the_gas_constant(self):
        rec = self.by_id["Wiki_7"]
        self.assertEqual(rec.formula, "A*exp(-E_a/(R*T))")
        self.assertEqual(rec.original_formula, "A*exp(-Ea/(R*T))")
        self.assertIn("R", rec.constant_names)
        self.assertEqual(rec.n_variables, 3)
        self.assertEqual(rec.n_constants, 1)
        self.assertEqual(rec.original_n_variables, 3)
        self.assertEqual(rec.original_n_constants, 1)
        self.assertEqual(rec.complexity, 9)
        tree = rec.to_json()["tree"]
        self.assertEqual(_find_label(tree, "R")["kind"], "constant")
        self.assertEqual(_find_label(tree, "A")["kind"], "variable")
        self.assertEqual(_find_label(tree, "E_a")["kind"], "variable")
        self.assertIsNone(_find_label(tree, "Ea"))
        self.assertEqual(_find_label(tree, "T")["kind"], "variable")
        self.assertIn(r"E_{a}", rec.to_json()["latex"])
        for wid in ("Wiki_8", "Wiki_9"):
            r = self.by_id[wid]
            self.assertIn("R", r.constant_names, wid)
            self.assertEqual(
                _find_label(r.to_json()["tree"], "R")["kind"], "constant", wid)

    def test_Wiki10_kb_h_and_R_are_constants(self):
        rec = self.by_id["Wiki_10"]
        self.assertEqual(rec.formula, "(kb*T/h)*exp(dG/(R*T))")
        self.assertEqual(rec.original_formula, "(kB*T/h)*exp(dG/(R*T))")
        self.assertIn("kb", rec.constant_names)
        self.assertNotIn("kB", rec.constant_names)
        self.assertIn("h", rec.constant_names)
        self.assertIn("R", rec.constant_names)
        self.assertEqual(rec.n_variables, 3)
        self.assertEqual(rec.n_constants, 3)
        tree = rec.to_json()["tree"]
        self.assertEqual(_find_label(tree, "kb")["kind"], "constant")
        self.assertIsNone(_find_label(tree, "kB"))
        self.assertEqual(_find_label(tree, "h")["kind"], "constant")
        self.assertEqual(_find_label(tree, "R")["kind"], "constant")
        tex = rec.to_json()["latex"]
        self.assertIn(r"k_{b}", tex)
        self.assertNotIn("kB", tex)
        self.assertNotIn(r"k_{B}", tex)
        self.assertEqual(rec.to_json()["hbar1"]["formula"],
                         "kb*T/(2*pi)*exp(dG/(R*T))")

    def test_Wiki11_k_is_boltzmann_kb(self):
        rec = self.by_id["Wiki_11"]
        self.assertEqual(rec.original_formula, "-k*T*log(Zn)")
        self.assertEqual(rec.formula, "-kb*T*log(Z_n)")
        self.assertIn("kb", rec.constant_names)
        self.assertEqual(rec.n_variables, 2)
        self.assertEqual(rec.n_constants, 1)
        tree = rec.to_json()["tree"]
        self.assertEqual(_find_label(tree, "kb")["kind"], "constant")
        self.assertIsNone(_find_label(tree, "k"))
        self.assertEqual(_find_label(tree, "T")["kind"], "variable")
        self.assertEqual(_find_label(tree, "Z_n")["kind"], "variable")
        self.assertIsNone(_find_label(tree, "Zn"))
        self.assertIn(r"k_{b}", rec.to_json()["latex"])
        self.assertIn(r"Z_{n}", rec.to_json()["latex"])

    def test_Wiki12_k_is_a_rate_not_boltzmann(self):
        rec = self.by_id["Wiki_12"]
        self.assertEqual(rec.formula, "1-exp(-k*t**n)")
        self.assertNotIn("kb", rec.constant_names)
        self.assertEqual(_find_label(rec.to_json()["tree"], "k")["kind"],
                         "variable")

    def test_Wiki13_R_is_gas_constant_c_is_not_light(self):
        rec = self.by_id["Wiki_13"]
        self.assertIn("R", rec.constant_names)
        self.assertNotIn("c", rec.constant_names)
        tree = rec.to_json()["tree"]
        self.assertEqual(_find_label(tree, "R")["kind"], "constant")
        self.assertEqual(_count_label_kind(tree, "R", "constant"), 3)
        self.assertEqual(_find_label(tree, "c")["kind"], "variable")
        self.assertEqual(rec.n_variables, 22)
        self.assertEqual(rec.n_constants, 3)
        self.assertEqual(rec.original_n_variables, 22)
        self.assertEqual(rec.original_n_constants, 3)
        self.assertIsNone(rec.to_json().get("c1"))
        nat = C.rewrite_units(rec, "natural")
        self.assertIn("c", C._name_occurrences(nat))
        self.assertNotIn("R", C._name_occurrences(nat))

    def test_Wiki16_h_is_height_not_planck(self):
        rec = self.by_id["Wiki_16"]
        self.assertEqual(rec.formula, "q+rho*g*h")
        self.assertIn("g", rec.constant_names)
        self.assertNotIn("h", rec.constant_names)
        tree = rec.to_json()["tree"]
        self.assertEqual(_find_label(tree, "g")["kind"], "constant")
        self.assertEqual(rec.original_n_variables, 3)
        self.assertEqual(rec.original_n_constants, 1)
        self.assertNotIn("hbar1", rec.to_json())
        nat = C.rewrite_units(rec, "natural")
        self.assertIn("h", C._name_occurrences(nat))
        self.assertNotIn("g", C._name_occurrences(nat))

    def test_Wiki31_R0_is_not_the_gas_constant(self):
        rec = self.by_id["Wiki_31"]
        self.assertNotIn("R", rec.constant_names)
        self.assertIsNone(_find_label(rec.to_json()["tree"], "R"))
        self.assertEqual(_find_label(rec.to_json()["tree"], "R0")["kind"],
                         "variable")

    def test_Wiki40_epsilon0_and_kb_are_constants(self):
        rec = self.by_id["Wiki_40"]
        self.assertIn("kb", rec.formula)
        self.assertNotIn("kB", rec.formula)
        self.assertIn("epsilon0", rec.constant_names)
        self.assertIn("kb", rec.constant_names)
        tree = rec.to_json()["tree"]
        self.assertEqual(_find_label(tree, "epsilon0")["kind"], "constant")
        self.assertEqual(_find_label(tree, "kb")["kind"], "constant")
        self.assertEqual(_find_label(tree, "epsilonr")["kind"], "variable")
        self.assertIn(r"k_{b}", rec.to_json()["latex"])

    def test_Wiki26_k_in_kT_is_boltzmann_kb(self):
        rec = self.by_id["Wiki_26"]
        self.assertEqual(rec.original_formula, "A*j**(-n)*exp(Q/(k*T))")
        self.assertEqual(rec.formula, "A*j**(-n)*exp(Q/(kb*T))")
        self.assertIn("kb", rec.constant_names)
        self.assertEqual(_find_label(rec.to_json()["tree"], "kb")["kind"],
                         "constant")
        self.assertIsNone(_find_label(rec.to_json()["tree"], "k"))

    def test_Wiki39_f_D_is_darcy_friction(self):
        rec = self.by_id["Wiki_39"]
        self.assertEqual(rec.original_formula, "fD*L*rho*V**2/D/2")
        self.assertEqual(rec.formula, "f_D*L*rho*V**2/D/2")
        tree = rec.to_json()["tree"]
        self.assertEqual(_find_label(tree, "f_D")["kind"], "variable")
        self.assertIsNone(_find_label(tree, "fD"))
        self.assertIn(r"f_{D}", rec.to_json()["latex"])

    def test_Wiki41_D_v_is_a_velocity(self):
        rec = self.by_id["Wiki_41"]
        self.assertEqual(rec.original_formula, "(1+Dv/c)*f0")
        self.assertEqual(rec.formula, "(1+D_v/c)*f0")
        self.assertIn("c", rec.constant_names)
        tree = rec.to_json()["tree"]
        self.assertEqual(_find_label(tree, "D_v")["kind"], "variable")
        self.assertIsNone(_find_label(tree, "Dv"))
        tex = rec.to_json()["latex"]
        self.assertIn(r"D_{v}", tex)
        self.assertIn(r"f_{0}", tex)
        c1 = rec.to_json().get("c1")
        self.assertIsNotNone(c1)
        self.assertEqual(c1["formula"], "(1+D_v)*f0")

    def test_Wiki3_vp_vs_are_subscripts(self):
        rec = self.by_id["Wiki_3"]
        self.assertEqual(rec.original_formula, "vp**2 - 4/3*vs**2")
        self.assertEqual(rec.formula, "v_p**2 - 4/3*v_s**2")
        tree = rec.to_json()["tree"]
        self.assertEqual(_find_label(tree, "v_p")["kind"], "variable")
        self.assertEqual(_find_label(tree, "v_s")["kind"], "variable")
        tex = rec.to_json()["latex"]
        self.assertIn(r"v_{p}", tex)
        self.assertIn(r"v_{s}", tex)

    def test_Wiki27_T_m_is_a_subscript(self):
        rec = self.by_id["Wiki_27"]
        self.assertEqual(rec.original_formula, "p*(0.457*Tm+8.128)")
        self.assertEqual(rec.formula, "p*(0.457*T_m+8.128)")
        self.assertEqual(_find_label(rec.to_json()["tree"], "T_m")["kind"],
                         "variable")
        self.assertIn(r"T_{m}", rec.to_json()["latex"])

    def test_Wiki41_c_is_a_constant(self):
        rec = self.by_id["Wiki_41"]
        self.assertIn("c", rec.constant_names)
        self.assertIsNotNone(rec.to_json().get("c1"))


class InflationCatalogTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.recs = C.seed_inflation()
        cls.by_id = {r.id: r for r in cls.recs}

    def test_released_keeps_all_seventy(self):
        self.assertEqual(len(self.recs), 70)
        self.assertEqual(len(C.select(self.recs, "released")), 70)
        self.assertEqual(len(C.select(self.recs, "primary")), 70)
        self.assertTrue(all(r.status == "original" for r in self.recs))
        self.assertTrue(all(r.corpus == "inflation" for r in self.recs))
        self.assertEqual(set(C._EI_NAMES), {r.id for r in self.recs})
        self.assertEqual(len(set(r.id for r in self.recs)), 70)

    def test_every_formula_parses(self):
        for rec in self.recs:
            tree = rec.to_json()["tree"]
            self.assertEqual(rec.complexity, M.tree_size(tree), rec.id)
            self.assertNotIn("_sym", rec.to_json()["latex"], rec.id)
            self.assertTrue(rec.name, rec.id)
            self.assertEqual(rec.source_url, C.INFLATION_PDF_URL, rec.id)
            self.assertEqual(rec.source_kind, "pdf", rec.id)
            self.assertEqual(rec.output, "V", rec.id)
            self.assertNotIn("c", C._name_occurrences(rec.formula), rec.id)
            if "M_Pl" in C._name_occurrences(rec.formula):
                self.assertIn("M_Pl", rec.constant_names, rec.id)
                self.assertEqual(
                    _find_label(tree, "M_Pl")["kind"], "constant", rec.id)

    def test_vfmi_recovers_blank_filename(self):
        rec = self.by_id["VFMI"]
        self.assertEqual(rec.original_id, "")
        self.assertEqual(rec.original_formula, C._VFMI_CSV)
        self.assertEqual(rec.formula, C.rename_ei_c(C._EI_FORMULA_PATCHES["VFMI"]))
        self.assertNotEqual(rec.formula, C.rename_ei_c(C._VFMI_CSV))
        self.assertEqual(rec.formula.count("(2-a)/(2*sqrt(3)*b)"), 2)
        self.assertEqual(rec.formula.count("(x/M_Pl)"), 2)
        self.assertNotIn("(2-a)/2*sqrt(3)*b", rec.formula)
        self.assertNotIn("))*x/M_Pl", rec.formula)
        tree = rec.to_json()["tree"]
        ones = _add_nodes_with_one(tree)
        self.assertEqual(len(ones), 2, rec.formula)
        for node in ones:
            self.assertIsNotNone(_find_label(node, "x"),
                                 "x/M_Pl must sit inside 1+(…), not multiply it")
            self.assertIsNotNone(_find_label(node, "M_Pl"))
        tex = rec.to_json()["latex"]
        self.assertIn(r"\frac{2 - a}{2 \sqrt{3} b}", tex)
        self.assertEqual(rec.name, "Mukhanov Inflation")
        self.assertIn("Filename", rec.notes)
        self.assertIn("M_Pl", rec.constant_names)
        self.assertNotIn("c", rec.constant_names)
        self.assertNotIn("M", rec.constant_names)
        self.assertEqual(_find_label(rec.to_json()["tree"], "M_Pl")["kind"],
                         "constant")
        self.assertEqual(_find_label(rec.to_json()["tree"], "M")["kind"],
                         "variable")

    def test_pai_has_no_csv_number(self):
        rec = self.by_id["PAI"]
        self.assertEqual(rec.original_id, "PAI")
        self.assertEqual(rec.name, "Pure Arctan Inflation")
        self.assertIn("Number blank", rec.notes)
        self.assertEqual(rec.source_locator, "PAI")
        self.assertNotIn("c", rec.constant_names)

    def test_six_bracket_repairs(self):
        self.assertEqual(set(M.EI_BRACKET_REPAIRS),
                         {"CNAI", "CNBI", "TMI", "CNCI", "HBI", "SATI"})
        for acr, repaired in M.EI_BRACKET_REPAIRS.items():
            rec = self.by_id[acr]
            self.assertEqual(rec.formula, C.rename_ei_c(repaired), acr)
            self.assertNotEqual(rec.original_formula, rec.formula, acr)
            self.assertIn("EI_BRACKET_REPAIRS", rec.notes, acr)
            M.expression_tree(rec.formula)
            with self.assertRaises(SyntaxError):
                M.expression_tree(rec.original_formula)

    def test_si_mpl_is_constant_m_is_amplitude(self):
        rec = self.by_id["SI"]
        self.assertEqual(rec.name, "Starobinsky Inflation")
        self.assertIn("M_Pl", rec.formula)
        self.assertNotIn("c", C._name_occurrences(rec.formula))
        self.assertIn("c", C._name_occurrences(rec.original_formula))
        self.assertIn("M_Pl", rec.constant_names)
        self.assertNotIn("c", rec.constant_names)
        self.assertNotIn("M", rec.constant_names)
        tree = rec.to_json()["tree"]
        self.assertEqual(_find_label(tree, "M_Pl")["kind"], "constant")
        self.assertIsNone(_find_label(tree, "c"))
        self.assertEqual(_find_label(tree, "M")["kind"], "variable")
        self.assertEqual(_find_label(tree, "x")["kind"], "variable")
        self.assertIn(r"M_{\mathrm{Pl}}", rec.to_json()["latex"])
        self.assertNotIn("c1", rec.to_json())

    def test_ni_has_no_planck_mass(self):
        rec = self.by_id["NI"]
        self.assertNotIn("c", C._name_occurrences(rec.formula))
        self.assertNotIn("M_Pl", C._name_occurrences(rec.formula))
        self.assertNotIn("c", rec.constant_names)
        self.assertNotIn("M_Pl", rec.constant_names)
        self.assertIsNone(rec.to_json().get("c1"))
        self.assertEqual(_find_label(rec.to_json()["tree"], "M")["kind"],
                         "variable")

    def test_rchi_follows_martin_16_and_sqrt6(self):
        rec = self.by_id["RCHI"]
        csv = ("(M**4)*(1-2*exp(-sqrt(2/3)*x/c)"
               "+x/c*(A/pi**2*16*sqrt(6)))")
        self.assertEqual(rec.original_formula, csv)
        self.assertEqual(
            rec.formula, C.rename_ei_c(C._EI_FORMULA_PATCHES["RCHI"]))
        self.assertIn("16*pi**2*sqrt(6)", rec.formula)
        self.assertNotIn("A/pi**2*16", rec.formula)
        tree = rec.to_json()["tree"]
        self.assertEqual(_find_label(tree, "A")["kind"], "constant")
        self.assertEqual(_find_label(tree, "M_Pl")["kind"], "constant")
        self.assertEqual(_find_label(tree, "pi")["kind"], "constant")
        self.assertEqual(_find_label(tree, "M")["kind"], "variable")
        self.assertIn("A", rec.constant_names)
        self.assertEqual(rec.n_variables, 3)
        self.assertEqual(rec.n_constants, 4)
        self.assertEqual(rec.original_n_variables, 3)
        self.assertEqual(rec.original_n_constants, 4)
        self.assertEqual(rec.complexity, rec.original_complexity)
        tex = rec.to_json()["latex"]
        self.assertIn(r"16", tex)
        self.assertIn(r"\sqrt{6}", tex)
        self.assertIn(r"\pi^{2}", tex)
        # 16 and √6 sit in the denominator of A / (16 π² √6).
        self.assertIn(r"16 \pi^{2} \sqrt{6}", tex)
        self.assertFalse(C.csv_disagreements(rec))
        self.assertIn("A is a constant", rec.notes)

    def test_rchi_A_is_a_constant(self):
        rec = self.by_id["RCHI"]
        tree = rec.to_json()["tree"]
        self.assertEqual(_find_label(tree, "A")["kind"], "constant")
        self.assertEqual(rec.n_variables, 3)
        self.assertEqual(rec.n_constants, 4)

    def test_csv_count_mismatches_are_the_eleven(self):
        flagged = {r.id: C.csv_disagreements(r) for r in self.recs
                   if C.csv_disagreements(r)}
        self.assertEqual(set(flagged), set(C._EI_COUNT_NOTES))
        self.assertEqual(set(flagged), {
            "HF1I", "CNAI", "SFI", "II", "SAII",
            "VFMI", "FI", "DEI", "SABI", "VHI", "DSI",
        })
        self.assertNotIn("RCHI", flagged)
        self.assertFalse(C.csv_disagreements(self.by_id["RCHI"]))
        self.assertEqual(set(flagged["HF1I"]), {"complexity"})
        self.assertEqual(set(flagged["SFI"]), {"variables"})
        self.assertEqual(set(flagged["VHI"]), {"variables"})
        self.assertEqual(set(flagged["DSI"]), {"variables"})
        self.assertEqual(set(flagged["SABI"]), {"numbers"})
        self.assertEqual(set(flagged["VFMI"]),
                         {"variables", "numbers", "complexity"})
        self.assertFalse(C.csv_disagreements(self.by_id["SI"]))
        self.assertFalse(C.csv_disagreements(self.by_id["LFI"]))
        for rec in self.recs:
            blob = rec.to_json()
            self.assertEqual(blob["csv_disagreements"],
                             C.csv_disagreements(rec), rec.id)
            if rec.id in C._EI_COUNT_NOTES:
                self.assertIn("CSV", rec.notes, rec.id)


class ExpressionScatterTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.payload = C.expression_scatter()

    def test_one_point_per_fitted_expression(self):
        by = {c["id"]: c for c in self.payload["corpora"]}
        expected = ["feynman", "wiki", "inflation", "leemis", "fr"]
        self.assertEqual(list(by), expected)
        for c in by.values():
            self.assertEqual(c["corpus_beta_estimator"], "multinomial MLE")
            self.assertLess(c["corpus_beta_p16"], c["corpus_beta"])
        self.assertEqual(by["feynman"]["n_fit"], 100)
        self.assertEqual(by["wiki"]["n"], 41 + len(C._WIKI_EXTENSION))
        self.assertEqual(by["wiki"]["n_fit"], by["wiki"]["n"])
        self.assertEqual(by["inflation"]["n_fit"], 70)
        self.assertEqual(by["leemis"]["n"], 61)
        self.assertEqual(by["leemis"]["n_fit"], 60)
        self.assertEqual(by["fr"]["n"], 10)
        self.assertEqual(by["fr"]["n_fit"], 10)
        self.assertEqual(
            [row["id"] for row in by["leemis"]["omitted"]],
            ["leemis.standard_uniform"],
        )
        for corpus in by.values():
            self.assertGreater(corpus["corpus_beta"], 0.1)
            self.assertLess(corpus["corpus_beta"], 0.6)
            self.assertNotIn("null", corpus)
            for point in corpus["points"]:
                self.assertGreaterEqual(point["k"], 1)
                self.assertGreaterEqual(point["beta"], 0)
                self.assertLess(point["beta"], 2)
        self.assertNotIn("population", self.payload)


def _add_nodes_with_one(node):
    """Add nodes whose one child is the number 1 (the VFMI 1+(…) factors)."""
    found = []
    if node.get("label") == "Add":
        kids = node.get("children") or []
        if any(c.get("kind") == "number" and c.get("label") == "1" for c in kids):
            found.append(node)
    for child in node.get("children") or []:
        found.extend(_add_nodes_with_one(child))
    return found


def _find_label(node, label):
    if node.get("label") == label:
        return node
    for child in node.get("children") or []:
        found = _find_label(child, label)
        if found:
            return found
    return None


def _count_kind(node, kind):
    n = 1 if node.get("kind") == kind else 0
    return n + sum(_count_kind(c, kind) for c in node.get("children") or [])


def _count_label_kind(node, label, kind):
    n = 1 if node.get("label") == label and node.get("kind") == kind else 0
    return n + sum(_count_label_kind(c, label, kind)
                   for c in node.get("children") or [])


if __name__ == "__main__":
    unittest.main()


class FeynmanLecturesCatalogTests(unittest.TestCase):
    """Feynman Lectures catalogs: public with a source link per expression."""

    def test_public_and_loadable(self):
        found = False
        for key in C.FEYNMAN_LECTURES:
            path = C.CATALOG_DIR / f"{key}.json"
            if not path.exists():
                continue
            found = True
            payload = __import__("json").loads(path.read_text())
            self.assertTrue(payload["public"])
            self.assertTrue(payload["release_formulas"])
            self.assertTrue(all(r["public"] is True for r in payload["records"]))
            self.assertTrue(all(r["source_url"].startswith("https://www.feynmanlectures.caltech.edu/")
                                for r in payload["records"]))
            recs = C.scatter_records(key)
            self.assertEqual(len(recs), payload["n_records"])
            ids = [r.id for r in recs]
            self.assertEqual(len(ids), len(set(ids)))
            for r in recs:
                self.assertGreaterEqual(r.complexity, 2, r.id)
        if not found:
            self.skipTest("no Feynman Lectures catalog present")

    def test_not_in_expression_scatter(self):
        self.assertFalse(set(C.FEYNMAN_LECTURES) & {k for k, _ in C.SCATTER_CORPORA})
