"""The one random draw on the contributions page.

Every curve_fit on that page is deterministic: the loss, the start
(β = 0.5, or α = 2 and b = 2), and the catalog are fixed. The sample of
10 Encyclopaedia Inflationaris expressions is the exception. It uses
SAMPLE_SEED so the median, the band, and the histogram come back the same.
"""
from __future__ import annotations

import numpy as np
from scipy.optimize import curve_fit

from . import metalaw as M
from .catalog import record_as_entry

SAMPLE_SEED = 20260930
N_DRAWS = 2000
SAMPLE_SIZE = 10


def inflation_sample_betas(entries) -> np.ndarray:
    """β for each draw of SAMPLE_SIZE expressions, in draw order."""

    def one(subset) -> float:
        counts = list(M.rank_counts(subset).values())
        obs = np.sort(np.asarray(counts, dtype=float))[::-1]
        obs = obs[obs > 0]
        rank = np.arange(1, obs.size + 1, dtype=float)
        logf = np.log(obs / obs.sum())

        def model(r, beta):
            logy = -beta * r
            return logy - np.logaddexp.reduce(logy)

        popt, _ = curve_fit(model, rank, logf, p0=[0.5])
        return float(popt[0])

    rng = np.random.default_rng(SAMPLE_SEED)
    idx = np.arange(len(entries))
    betas = np.empty(N_DRAWS)
    for i in range(N_DRAWS):
        take = rng.choice(idx, size=SAMPLE_SIZE, replace=False)
        betas[i] = one([entries[j] for j in take])
    return betas


def inflation_entries():
    from .catalog import scatter_records
    return [record_as_entry(r) for r in scatter_records("inflation")]
