"""β of one expression, on the ranks that expression actually uses."""
from __future__ import annotations

import warnings

import numpy as np
from scipy.optimize import OptimizeWarning

from . import metalaw as M


def expression_beta(entry: M.EIEntry):
    """Multinomial MLE of β on the positive ranks of one expression.

    Same estimator as the corpus β: the expression's operator instances are
    treated as draws over its own ranks. (Log least squares, the former
    choice, is still M.fit_beta_logls.)

    A type the expression does not use is absent, not a zero count.
    Returns (beta, n_ranks). beta is None when fewer than two ranks appear.
    """
    vals = [v for v in M.rank_counts([entry]).values() if v > 0]
    if len(vals) < 2:
        return None, len(vals)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", OptimizeWarning)
        warnings.simplefilter("ignore", RuntimeWarning)
        try:
            beta = M.fit_beta_mle(vals)
        except (RuntimeError, ValueError):
            return None, len(vals)
    if not np.isfinite(beta):
        return None, len(vals)
    return float(beta), len(vals)
