import unittest

import numpy as np

from physics_metalaw import metalaw as M
from physics_metalaw.beta import expression_beta


class MultinomialMLETest(unittest.TestCase):
    def test_recovers_beta_from_exact_counts(self):
        r = np.arange(1, 21)
        p = np.exp(-0.33 * r); p /= p.sum()
        self.assertAlmostEqual(M.fit_beta_mle(1e6 * p), 0.33, places=5)

    def test_two_ranks_closed_form(self):
        self.assertAlmostEqual(M.fit_beta_mle([2, 1]), np.log(2), places=6)

    def test_order_independent(self):
        self.assertAlmostEqual(M.fit_beta_mle([1, 5, 3]), M.fit_beta_mle([5, 3, 1]), places=10)

    def test_fisher_matches_simulation(self):
        rng = np.random.default_rng(1)
        r = np.arange(1, 16)
        p = np.exp(-0.3 * r); p /= p.sum()
        bs = []
        for _ in range(400):
            c = rng.multinomial(2000, p)
            c = c[c > 0]
            bs.append(M.fit_beta_mle(c))
        self.assertAlmostEqual(np.std(bs) / M.fisher_sigma(0.3, 15, 2000), 1.0, delta=0.2)

    def test_bootstrap_deterministic_and_bracketing(self):
        from physics_metalaw.catalog import record_as_entry, scatter_records
        e = [record_as_entry(r) for r in scatter_records("feynman")]
        a = M.beta_with_errors(e, reps=300)
        b = M.beta_with_errors(e, reps=300)
        self.assertEqual(a, b)
        self.assertLess(a["boot_p16"], a["beta"])
        self.assertLess(a["beta"], a["boot_p84"])
        self.assertGreater(a["boot_sigma"], 0)



if __name__ == "__main__":
    unittest.main()
