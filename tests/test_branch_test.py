"""The Feynman corpus of scripts/branch_test.py: deduplication, branch sizes, β (no bootstrap)."""
import importlib.util
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
_spec = importlib.util.spec_from_file_location("branch_test", ROOT / "scripts" / "branch_test.py")
BT = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(BT)


class FeynmanCorpusTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.rows, cls.dedupe, cls.added = BT.feynman_corpus()
        cls.groups = BT.groups_of(cls.rows)

    def test_sizes(self):
        self.assertEqual(len(self.rows), 1016)
        self.assertEqual(self.added, {"v1": 278, "v2": 443, "v3": 195})
        self.assertEqual({k: len(v) for k, v in self.groups.items()},
                         {"SR": 89, "QM": 221, "EM": 312, "Waves+optics": 126, "Classical": 105,
                          "Thermo": 73, "Fluid+continuum": 42})

    def test_fixes_applied(self):
        import csv, json
        from physics_metalaw.catalog import CATALOG_DIR
        with open(BT.ROOT / "data" / "feynman_lectures_fixes.csv", newline="") as f:
            fixes = {r["id"]: r for r in csv.DictReader(f)}
        recs = {d["id"]: d for v in (1, 2, 3)
                for d in json.loads((CATALOG_DIR / f"feynman_lectures_v{v}.json").read_text())["records"]}
        for i, fx in fixes.items():
            if fx["action"] == "reject":
                self.assertNotIn(i, recs)
            elif fx["action"] == "replace":
                self.assertEqual(recs[i]["formula"], fx["new_formula"])
        rel = [d for d in recs.values() if (d["volume"], d["chapter"]) in BT.BF.BRANCH_OVERRIDE]
        self.assertTrue(rel)
        self.assertTrue(all(d["branch"] == "Special relativity" for d in rel))

    def test_no_benchmark_double_count(self):
        ids = [r.id for r, _, _ in self.rows]
        self.assertEqual(len(ids), len(set(ids)))

    def test_branch_beta(self):
        C = BT.Counts([r for r, _, _ in self.rows])
        want = {"SR": 0.308, "QM": 0.291, "EM": 0.352, "Waves+optics": 0.343, "Classical": 0.359,
                "Thermo": 0.408, "Fluid+continuum": 0.433}
        for k, b in want.items():
            self.assertAlmostEqual(BT.mle([C.obs(self.groups[k])])[0], b, places=3, msg=k)
        stat = BT.lrt(C, [self.groups[k] for k in BT.ORDER])[0]
        self.assertAlmostEqual(stat, 102.05, places=1)


if __name__ == "__main__":
    unittest.main()
