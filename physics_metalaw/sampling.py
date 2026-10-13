"""The one random draw behind the Encyclopaedia Inflationaris figure in the writeup.

Every other β is the multinomial MLE (M.fit_beta_mle); its bootstrap
uses M.BOOTSTRAP_SEED. The sample of
10 Encyclopaedia Inflationaris expressions is the exception. It uses
SAMPLE_SEED so the median, the band, and the histogram come back the same.
"""
from __future__ import annotations

import numpy as np

from . import metalaw as M
from .catalog import record_as_entry

SAMPLE_SEED = 20260930
N_DRAWS = 2000
SAMPLE_SIZE = 10


def inflation_sample_betas(entries) -> np.ndarray:
    """β for each draw of SAMPLE_SIZE expressions, in draw order."""

    def one(subset) -> float:
        return M.fit_beta_mle(list(M.rank_counts(subset).values()))

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
