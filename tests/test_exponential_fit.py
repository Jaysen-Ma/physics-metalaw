"""Focused checks for exponential fitting, tail exclusion and simultaneous intervals.

These tests check the calculations. They do not estimate confidence-interval coverage.
"""
import unittest

import numpy as np

from physics_metalaw.exponential_fit import (
    BOOTSTRAP_SEED,
    QUANTILE_METHOD,
    STATUS_DEGENERATE,
    STATUS_FAILED,
    STATUS_OK,
    bootstrap_type_matrix,
    fit_exponential,
    loss,
    pairwise_differences,
    simultaneous_intervals,
    tail_exclusion,
)


def _exponential_counts(beta: float, n_ranks: int, total: float) -> np.ndarray:
    rank = np.arange(1, n_ranks + 1, dtype=float)
    p = np.exp(-beta * rank)
    p /= p.sum()
    return total * p


class ExponentialFitTests(unittest.TestCase):
    def test_recovers_known_exponential_frequencies(self):
        counts = _exponential_counts(0.33, 12, 1e6)
        fit = fit_exponential(counts, threshold=None)
        self.assertEqual(fit.status, STATUS_OK)
        self.assertAlmostEqual(fit.beta, 0.33, places=5)

    def test_exclusion_keeps_the_same_beta_when_frequencies_match(self):
        counts = _exponential_counts(0.33, 12, 1e6)
        fit = fit_exponential(counts, threshold=0.01)
        self.assertGreaterEqual(fit.n_retained, 2)
        self.assertAlmostEqual(fit.beta, 0.33, places=4)

    def test_order_of_input_counts_does_not_matter(self):
        counts = np.array([1, 4, 9, 2, 7], dtype=float)
        self.assertAlmostEqual(
            fit_exponential(counts, None).beta,
            fit_exponential(counts[::-1], None).beta,
            places=10,
        )


class TailExclusionTests(unittest.TestCase):
    def test_no_exclusion_retains_every_positive_rank(self):
        counts = np.array([10, 3, 1], dtype=float)
        exclusion = tail_exclusion(counts, None)
        self.assertTrue(exclusion.retained.all())
        self.assertEqual(exclusion.n_excluded, 0)
        self.assertEqual(exclusion.excluded_mass, 0.0)

    def test_first_tail_group_above_one_percent_is_excluded_together(self):
        # Three tied counts of 10. Their mass is 0.03, so the whole tie is excluded.
        counts = np.array([970, 10, 10, 10], dtype=float)
        exclusion = tail_exclusion(counts, 0.01)
        self.assertEqual(exclusion.retained_ranks, [1])
        self.assertEqual(exclusion.n_excluded, 3)
        self.assertAlmostEqual(exclusion.excluded_mass, 0.03)

    def test_exact_one_percent_does_not_stop_the_accumulation(self):
        # The two 5s are exactly 1% of 1000, which is not strictly above 1%.
        # The next group, 40, is included in the excluded tail.
        counts = np.array([800, 150, 40, 5, 5], dtype=float)
        exclusion = tail_exclusion(counts, 0.01)
        self.assertEqual(exclusion.retained_ranks, [1, 2])
        self.assertEqual(exclusion.n_excluded, 3)
        self.assertAlmostEqual(exclusion.excluded_mass, 0.05)

    def test_threshold_crossing_stops_before_the_head(self):
        counts = np.array([900, 80, 20], dtype=float)
        exclusion = tail_exclusion(counts, 0.01)
        self.assertEqual(exclusion.retained_ranks, [1, 2])
        self.assertAlmostEqual(exclusion.excluded_mass, 0.02)

    def test_normalisation_uses_original_ranks_and_total(self):
        counts = np.array([900, 80, 20], dtype=float)
        exclusion = tail_exclusion(counts, 0.01)
        beta = 0.2
        rank = np.arange(1, 4, dtype=float)
        logp = -beta * rank
        logp -= np.logaddexp.reduce(logp)
        logf = np.log(counts / counts.sum())
        expected = float(np.sum((logf[exclusion.retained] - logp[exclusion.retained]) ** 2))
        self.assertAlmostEqual(loss(beta, counts, exclusion.retained), expected)
        # The excluded rank is absent from the loss, and f_2 still uses N = 1000.
        self.assertNotAlmostEqual(np.log(80 / 980), np.log(80 / 1000))


class BootstrapRuleTests(unittest.TestCase):
    def test_resample_reranks_and_reapplies_exclusion(self):
        matrix = np.array([
            [50, 0, 0],
            [0, 1, 1],
        ], dtype=float)
        full = fit_exponential(matrix.sum(axis=0), 0.01)
        self.assertGreater(full.n_excluded, 0)

        class OnlyFirst:
            def integers(self, low, high, size):
                return np.zeros(size, dtype=int)

        _fit, betas, statuses = bootstrap_type_matrix(
            matrix, n_boot=1, threshold=0.01, rng=OnlyFirst())
        doubled = fit_exponential(matrix[0] * 2, 0.01)
        self.assertEqual(statuses[0], doubled.status)
        if doubled.beta is None:
            self.assertTrue(np.isnan(betas[0]))
        else:
            self.assertAlmostEqual(betas[0], doubled.beta)
        self.assertNotEqual(full.n_excluded, doubled.n_excluded)

    def test_reproducible_with_the_recorded_seed(self):
        matrix = np.array([[4, 1, 0], [3, 2, 1], [8, 0, 0], [1, 1, 1]], dtype=float)
        first = bootstrap_type_matrix(matrix, n_boot=20, seed=BOOTSTRAP_SEED, threshold=0.01)
        second = bootstrap_type_matrix(matrix, n_boot=20, seed=BOOTSTRAP_SEED, threshold=0.01)
        np.testing.assert_array_equal(first[1], second[1])
        self.assertEqual(list(first[2]), list(second[2]))

    def test_failed_draw_is_kept_as_missing(self):
        # One token in total: exclusion of any positive tail leaves too few ranks,
        # and a resample cannot invent a second rank.
        matrix = np.array([[1.0, 0.0]], dtype=float)
        fit, betas, statuses = bootstrap_type_matrix(matrix, n_boot=5, seed=1, threshold=None)
        self.assertEqual(fit.status, STATUS_DEGENERATE)
        self.assertTrue(np.isnan(betas).all())
        self.assertTrue(np.all(statuses == STATUS_DEGENERATE))


class PairwiseAndSimultaneousTests(unittest.TestCase):
    def test_identical_inputs_have_zero_difference(self):
        rows = pairwise_differences({"A": 0.34, "B": 0.34}, ["A", "B"])
        self.assertEqual(rows[0]["delta"], 0.0)

    def test_direction_follows_the_input_rates(self):
        rows = pairwise_differences({"A": 0.40, "B": 0.30}, ["A", "B"])
        self.assertGreater(rows[0]["delta"], 0.0)
        flipped = pairwise_differences({"A": 0.30, "B": 0.40}, ["A", "B"])
        self.assertLess(flipped[0]["delta"], 0.0)

    def test_intervals_are_centred_at_the_original_estimate(self):
        delta = np.array([0.10, -0.05])
        draws = np.array([
            [0.10, -0.05],
            [0.12, -0.04],
            [0.08, -0.06],
            [0.11, -0.05],
        ])
        result = simultaneous_intervals(delta, draws)
        self.assertEqual(result["centring"], "original_estimate")
        self.assertEqual(result["quantile_method"], QUANTILE_METHOD)
        for interval, hat in zip(result["intervals"], delta):
            lo, hi = interval["simultaneous_ci95"]
            self.assertAlmostEqual(0.5 * (lo + hi), hat, places=10)
            self.assertFalse(interval["zero_se"])

    def test_zero_standard_error_collapses_to_the_point(self):
        delta = np.array([0.2, 0.0])
        draws = np.column_stack([
            np.full(10, 0.2),
            np.linspace(-0.01, 0.01, 10),
        ])
        result = simultaneous_intervals(delta, draws)
        point = result["intervals"][0]
        self.assertTrue(point["zero_se"])
        self.assertEqual(point["simultaneous_ci95"], (0.2, 0.2))
        self.assertTrue(point["excludes_zero"])
        self.assertEqual(result["interval_validity"], "approximate_simultaneous")

    def test_zero_standard_error_that_is_not_the_original_estimate_is_unresolved(self):
        delta = np.array([0.2])
        draws = np.full((8, 1), 0.4)
        result = simultaneous_intervals(delta, draws)
        self.assertIsNone(result["critical_value"])
        self.assertIsNone(result["intervals"][0]["simultaneous_ci95"][0])
        self.assertEqual(result["interval_validity"], "unresolved")

    def test_incomplete_replicate_is_counted(self):
        delta = np.array([0.1, 0.2])
        draws = np.array([
            [0.1, 0.2],
            [0.12, np.nan],
            [0.08, 0.22],
        ])
        result = simultaneous_intervals(delta, draws)
        self.assertEqual(result["n_incomplete"], 1)
        self.assertEqual(result["interval_validity"], "unresolved")



class UntrimmedRegressionTests(unittest.TestCase):
    def test_default_retains_all_positive_ranks(self):
        fit = fit_exponential([339,153,123,107,103,72,70,49,29,9,9,7,2,1])
        self.assertEqual(fit.n_excluded, 0)
        self.assertEqual(fit.n_retained, 14)
        self.assertAlmostEqual(fit.beta, 0.382469242, places=6)

    def test_default_is_stable_to_breaking_the_nine_count_tie(self):
        a = fit_exponential([339,153,123,107,103,72,70,49,29,9,9,7,2,1])
        b = fit_exponential([339,153,123,107,103,72,70,49,29,10,8,7,2,1])
        self.assertLess(abs(a.beta-b.beta), .01)

    def test_dictionary_uses_counts_not_operator_names(self):
        a = fit_exponential({'x':50,'plus':20,'sin':5})
        b = fit_exponential([50,20,5])
        self.assertEqual(a.beta,b.beta)

    def test_cluster_multiplicity_is_shared_across_rows(self):
        from physics_metalaw.exponential_fit import stratified_cluster_weights
        class PickFirst:
            def integers(self, low, high, size):
                return np.zeros(size,dtype=int)
        # Rows 0 and 2 may belong to different physics branches, but share
        # a source chapter and therefore must get the same multiplicity.
        clusters = [np.array([0,2]), np.array([1]), np.array([3,4])]
        weights = stratified_cluster_weights(clusters,[np.array([0,1]),np.array([2])],5,PickFirst())
        np.testing.assert_array_equal(weights,[2,0,2,1,1])

    def test_coverage_generator_has_known_exponential_marginals(self):
        from scripts.exponential_coverage_check import exponential_probabilities
        p=exponential_probabilities(12)
        self.assertAlmostEqual(p.sum(),1.0)
        np.testing.assert_allclose(p[1:]/p[:-1],np.exp(-.35))

    def test_monte_carlo_summary_counts_unresolved_datasets(self):
        from scripts.exponential_coverage_check import binomial_summary
        result=binomial_summary([True,False,None])
        self.assertEqual(result['n'],2)
        self.assertEqual(result['n_missing'],1)
        self.assertEqual(result['rate'],.5)

if __name__ == "__main__":
    unittest.main()
