"""Unweighted log-frequency least squares for a normalised exponential.

The primary analysis minimises

    sum_{r retained} (log f_r - log p_r(beta))^2

with p_r(beta) = exp(-beta r) / sum_{s=1}^R exp(-beta s). Ranks, N and R are
fixed before rare-operator exclusion. Exclusion drops residual terms only.
Frequencies stay f_r = n_r / N and the exponential stays normalised over the
original R ranks.

This module does not count operators. Callers pass positive rank counts, or a
matrix of per-equation operator counts.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

BOUNDS = (0.0, 5.0)
BOUND_TOL = 1e-4
PRIMARY_THRESHOLD = None  # All positive-count ranks enter the primary fit.
MIN_RETAINED = 2
QUANTILE_METHOD = "linear"  # numpy.quantile, Hyndman-Fan type 7
BOOTSTRAP_SEED = 20261013
BOOTSTRAP_REPS = 10000

STATUS_OK = "ok"
STATUS_BOUNDARY = "boundary"
STATUS_DEGENERATE = "degenerate"
STATUS_FAILED = "failed"


@dataclass
class Exclusion:
    """Which original ranks (1-based) remain in the fitting loss."""

    retained: np.ndarray  # bool, length R, True = kept in the loss
    excluded_mass: float
    n_excluded: int
    n_retained: int
    threshold: float | None

    @property
    def retained_ranks(self) -> list[int]:
        return (np.flatnonzero(self.retained) + 1).tolist()


@dataclass
class ExponentialFit:
    beta: float | None
    status: str
    n_instances: int
    n_ranks: int
    n_retained: int
    n_excluded: int
    excluded_mass: float
    retained_ranks: list[int] = field(default_factory=list)
    near_bound: bool = False
    threshold: float | None = PRIMARY_THRESHOLD
    counts: np.ndarray | None = None  # positive counts, descending

    def as_dict(self) -> dict:
        return {
            "beta": None if self.beta is None else float(self.beta),
            "status": self.status,
            "n_instances": int(self.n_instances),
            "n_ranks": int(self.n_ranks),
            "n_retained": int(self.n_retained),
            "n_excluded": int(self.n_excluded),
            "excluded_mass": float(self.excluded_mass),
            "retained_ranks": list(self.retained_ranks),
            "near_bound": bool(self.near_bound),
            "threshold": self.threshold,
        }


def positive_descending(counts) -> np.ndarray:
    """Positive counts sorted from most frequent to least frequent."""
    arr = np.asarray(list(counts.values()) if isinstance(counts, dict) else counts, dtype=float)
    arr = arr[np.isfinite(arr) & (arr > 0)]
    if arr.size == 0:
        return arr
    return np.sort(arr)[::-1]


def tail_exclusion(counts_desc: np.ndarray, threshold: float | None) -> Exclusion:
    """Drop a tail whose cumulative frequency first exceeds ``threshold``.

    Walk from the least frequent type toward the most frequent. Types that
    share a count are one group and are kept or dropped together. The group
    that first pushes the accumulated tail frequency strictly above
    ``threshold`` is excluded, along with every rarer group. ``threshold is
    None`` retains every positive-count rank.

    Comparison uses the integer count total: a group is past the threshold
    when its cumulative count is strictly greater than ``threshold * N``.
    """
    counts_desc = np.asarray(counts_desc, dtype=float)
    r = int(counts_desc.size)
    retained = np.ones(r, dtype=bool)
    if threshold is None or r == 0:
        return Exclusion(retained, 0.0, 0, r, threshold)
    total = float(counts_desc.sum())
    if total <= 0:
        return Exclusion(retained, 0.0, 0, 0, threshold)
    limit = float(threshold) * total
    excluded = np.zeros(r, dtype=bool)
    cum = 0.0
    i = r - 1
    while i >= 0:
        count = counts_desc[i]
        j = i
        while j >= 0 and counts_desc[j] == count:
            j -= 1
        cum += float(counts_desc[j + 1:i + 1].sum())
        excluded[j + 1:i + 1] = True
        i = j
        if cum > limit:
            break
    retained = ~excluded
    mass = float(counts_desc[excluded].sum() / total) if excluded.any() else 0.0
    return Exclusion(retained, mass, int(excluded.sum()), int(retained.sum()), float(threshold))


def _log_model(beta: float, n_ranks: int) -> np.ndarray:
    rank = np.arange(1, n_ranks + 1, dtype=float)
    logp = -beta * rank
    logp -= np.logaddexp.reduce(logp)
    return logp


def loss(beta: float, counts_desc: np.ndarray, retained: np.ndarray) -> float:
    """Unweighted squared log-frequency error on the retained original ranks.

    ``p_r`` is normalised over all ``R`` original ranks. ``f_r = n_r / N``
    uses the original total. Ranks that ``retained`` marks False contribute
    nothing.
    """
    counts_desc = np.asarray(counts_desc, dtype=float)
    if int(retained.sum()) == 0:
        return float("nan")
    logf = np.log(counts_desc / counts_desc.sum())
    resid = logf[retained] - _log_model(beta, counts_desc.size)[retained]
    return float(resid @ resid)


def _near_bound(beta: float) -> bool:
    lo, hi = BOUNDS
    return beta <= lo + BOUND_TOL or beta >= hi - BOUND_TOL


def fit_exponential(counts, threshold: float | None = PRIMARY_THRESHOLD) -> ExponentialFit:
    """Fit beta by unweighted log-frequency least squares."""
    from scipy.optimize import minimize_scalar

    counts_desc = positive_descending(counts)
    if counts_desc.size == 0:
        return ExponentialFit(None, STATUS_FAILED, 0, 0, 0, 0, 0.0, threshold=threshold)
    exclusion = tail_exclusion(counts_desc, threshold)
    total = int(round(float(counts_desc.sum())))
    base = dict(
        n_instances=total,
        n_ranks=int(counts_desc.size),
        n_retained=exclusion.n_retained,
        n_excluded=exclusion.n_excluded,
        excluded_mass=exclusion.excluded_mass,
        retained_ranks=exclusion.retained_ranks,
        threshold=threshold,
        counts=counts_desc,
    )
    if exclusion.n_retained < MIN_RETAINED:
        return ExponentialFit(None, STATUS_DEGENERATE, **base)
    retained = exclusion.retained

    def objective(beta: float) -> float:
        return loss(beta, counts_desc, retained)

    try:
        res = minimize_scalar(objective, bounds=BOUNDS, method="bounded", options={"xatol": 1e-10})
    except (ValueError, RuntimeError):
        return ExponentialFit(None, STATUS_FAILED, **base)
    if not res.success or not np.isfinite(res.x):
        return ExponentialFit(None, STATUS_FAILED, **base)
    beta = float(res.x)
    status = STATUS_BOUNDARY if _near_bound(beta) else STATUS_OK
    return ExponentialFit(beta, status, near_bound=status == STATUS_BOUNDARY, **base)


def fit_from_type_counts(type_counts: np.ndarray, threshold: float | None = PRIMARY_THRESHOLD) -> ExponentialFit:
    """Fit from one vector of per-type counts (zeros allowed)."""
    return fit_exponential(type_counts, threshold)


def _fit_positive(counts_desc: np.ndarray, threshold: float | None) -> tuple[float, str]:
    """Lean fit used inside the bootstrap. Returns (beta or nan, status)."""
    from scipy.optimize import minimize_scalar

    exclusion = tail_exclusion(counts_desc, threshold)
    if exclusion.n_retained < MIN_RETAINED:
        return float("nan"), STATUS_DEGENERATE
    retained = exclusion.retained
    try:
        res = minimize_scalar(
            lambda beta: loss(beta, counts_desc, retained),
            bounds=BOUNDS, method="bounded", options={"xatol": 1e-8})
    except (ValueError, RuntimeError):
        return float("nan"), STATUS_FAILED
    if not res.success or not np.isfinite(res.x):
        return float("nan"), STATUS_FAILED
    beta = float(res.x)
    if _near_bound(beta):
        return beta, STATUS_BOUNDARY
    return beta, STATUS_OK


def _fit_type_vector(type_counts: np.ndarray, threshold: float | None) -> tuple[float, str]:
    counts = np.asarray(type_counts, dtype=float)
    counts = counts[counts > 0]
    if counts.size == 0:
        return float("nan"), STATUS_FAILED
    return _fit_positive(np.sort(counts)[::-1], threshold)


def _weights(index: np.ndarray, n_rows: int) -> np.ndarray:
    return np.bincount(index, minlength=n_rows).astype(float)


def summarise_bootstrap(original: ExponentialFit, betas: np.ndarray, statuses: np.ndarray) -> dict:
    """Summarise bootstrap draws. Failed draws stay in ``betas`` as NaN.

    Percentile intervals and the standard error use the finite draws. If any
    draw failed or was degenerate, ``interval_validity`` is ``unresolved``:
    the interval is then only a summary of the successful draws.
    """
    betas = np.asarray(betas, dtype=float)
    statuses = np.asarray(statuses)
    finite = betas[np.isfinite(betas)]
    n_failed = int(np.sum(statuses == STATUS_FAILED))
    n_degenerate = int(np.sum(statuses == STATUS_DEGENERATE))
    n_boundary = int(np.sum(statuses == STATUS_BOUNDARY))
    n_ok = int(np.sum(statuses == STATUS_OK))
    unresolved = n_failed + n_degenerate > 0
    out = {
        "beta": original.beta,
        "status": original.status,
        "near_bound": original.near_bound,
        "n_instances": original.n_instances,
        "n_ranks": original.n_ranks,
        "n_retained": original.n_retained,
        "n_excluded": original.n_excluded,
        "excluded_mass": original.excluded_mass,
        "retained_ranks": original.retained_ranks,
        "threshold": original.threshold,
        "n_boot": int(betas.size),
        "n_finite": int(finite.size),
        "n_ok": n_ok,
        "n_boundary": n_boundary,
        "n_degenerate": n_degenerate,
        "n_failed": n_failed,
        "interval_validity": "unresolved" if unresolved else "approximate_percentile",
        "quantile_method": QUANTILE_METHOD,
    }
    if finite.size == 0 or original.beta is None:
        out.update(boot_se=None, boot_bias=None, ci95=(None, None))
        return out
    se = float(finite.std(ddof=1)) if finite.size > 1 else 0.0
    lo, hi = np.quantile(finite, [0.025, 0.975], method=QUANTILE_METHOD)
    out.update(
        boot_se=se,
        boot_bias=float(finite.mean() - original.beta),
        ci95=(float(lo), float(hi)),
    )
    return out


def bootstrap_type_matrix(
    matrix: np.ndarray,
    *,
    n_boot: int = BOOTSTRAP_REPS,
    seed: int = BOOTSTRAP_SEED,
    threshold: float | None = PRIMARY_THRESHOLD,
    clusters: list[np.ndarray] | None = None,
    rng: np.random.Generator | None = None,
) -> tuple[ExponentialFit, np.ndarray, np.ndarray]:
    """Resample rows of a per-equation type-count matrix and refit.

    Each row is one equation. A bootstrap draw resamples rows with
    replacement. If ``clusters`` is given, draws resample those row-groups
    instead, keeping every equation in a sampled cluster together. The
    original fit is computed on all rows. Returns ``(fit, betas, statuses)``.
    Status entries are the strings ``ok``, ``boundary``, ``degenerate`` and
    ``failed``. Non-numeric fits are NaN and are not deleted.
    """
    matrix = np.asarray(matrix, dtype=float)
    if matrix.ndim != 2:
        raise ValueError("matrix must have one row per equation")
    n_rows = matrix.shape[0]
    original = fit_from_type_counts(matrix.sum(axis=0), threshold)
    betas = np.full(n_boot, np.nan)
    statuses = np.empty(n_boot, dtype=object)
    if n_rows == 0:
        statuses[:] = STATUS_FAILED
        return original, betas, statuses
    generator = rng if rng is not None else np.random.default_rng(seed)
    if clusters is None:
        for i in range(n_boot):
            drawn = generator.integers(0, n_rows, n_rows)
            beta, status = _fit_type_vector(_weights(drawn, n_rows) @ matrix, threshold)
            betas[i] = beta
            statuses[i] = status
        return original, betas, statuses
    n_clusters = len(clusters)
    if n_clusters == 0:
        statuses[:] = STATUS_FAILED
        return original, betas, statuses
    for i in range(n_boot):
        picked = generator.integers(0, n_clusters, n_clusters)
        weight = np.zeros(n_rows)
        for cluster_id in picked:
            weight[clusters[cluster_id]] += 1.0
        beta, status = _fit_type_vector(weight @ matrix, threshold)
        betas[i] = beta
        statuses[i] = status
    return original, betas, statuses


def stratified_cluster_weights(clusters, strata, n_rows, rng):
    """Resample clusters within strata; one weight applies to all their rows."""
    weights = np.zeros(n_rows)
    for cluster_ids in strata:
        for picked in rng.integers(0, len(cluster_ids), len(cluster_ids)):
            weights[clusters[cluster_ids[picked]]] += 1
    return weights


def pairwise_differences(betas: dict[str, float], order: list[str]) -> list[dict]:
    """All ordered pairs i < j in ``order``, delta = beta_i - beta_j."""
    rows = []
    for i, left in enumerate(order):
        for right in order[i + 1:]:
            rows.append({
                "a": left,
                "b": right,
                "delta": float(betas[left] - betas[right]),
            })
    return rows


def simultaneous_intervals(
    delta_hat: np.ndarray,
    draws: np.ndarray,
    *,
    level: float = 0.95,
) -> dict:
    """Centred max-statistic simultaneous intervals.

    ``draws`` has shape (n_boot, n_comparisons). A draw is centred by
    subtracting the original estimate, not the bootstrap mean. ``s_j`` is the
    sample standard deviation of the finite draws of comparison ``j``.

    Replicates that are missing any comparison are counted and left out of
    the maximum statistic. If any comparison has a zero standard deviation,
    that comparison is reported explicitly. A zero standard deviation with
    every finite draw equal to the original estimate contributes nothing to
    the maximum, and the interval collapses to the point estimate. Any other
    zero-standard-deviation case is marked unresolved.
    """
    delta_hat = np.asarray(delta_hat, dtype=float)
    draws = np.asarray(draws, dtype=float)
    if draws.ndim != 2 or draws.shape[1] != delta_hat.size:
        raise ValueError("draws must be (n_boot, n_comparisons)")
    n_boot, n_comp = draws.shape
    finite = np.isfinite(draws)
    complete = finite.all(axis=1)
    n_incomplete = int((~complete).sum())
    s = np.empty(n_comp)
    zero_se = np.zeros(n_comp, dtype=bool)
    unresolved_se = np.zeros(n_comp, dtype=bool)
    for j in range(n_comp):
        column = draws[finite[:, j], j]
        if column.size <= 1:
            s[j] = 0.0 if column.size == 1 else float("nan")
        else:
            s[j] = float(column.std(ddof=1))
        if column.size == 0 or not np.isfinite(s[j]):
            unresolved_se[j] = True
            zero_se[j] = True
            continue
        if s[j] == 0.0:
            zero_se[j] = True
            if not np.allclose(column, delta_hat[j]):
                unresolved_se[j] = True
    validity = "approximate_simultaneous"
    if n_incomplete or unresolved_se.any() or not np.isfinite(delta_hat).all():
        validity = "unresolved"
    m_values = np.array([])
    if complete.any() and np.isfinite(s).all() and not unresolved_se.any():
        scale = s.copy()
        # A comparison with no bootstrap variation adds a zero term.
        varying = scale > 0
        centred = draws[complete] - delta_hat
        if varying.any():
            m_values = np.max(np.abs(centred[:, varying] / scale[varying]), axis=1)
        else:
            m_values = np.zeros(int(complete.sum()))
    if m_values.size == 0 or validity == "unresolved" and unresolved_se.any():
        c = None
    else:
        c = float(np.quantile(m_values, level, method=QUANTILE_METHOD))
    intervals = []
    for j in range(n_comp):
        if c is None or not np.isfinite(s[j]) or unresolved_se[j]:
            lo = hi = None
        elif s[j] == 0.0:
            lo = hi = float(delta_hat[j])
        else:
            half = c * s[j]
            lo = float(delta_hat[j] - half)
            hi = float(delta_hat[j] + half)
        excludes = None if lo is None else bool(lo > 0 or hi < 0)
        intervals.append({
            "delta": float(delta_hat[j]),
            "boot_se": None if not np.isfinite(s[j]) else float(s[j]),
            "zero_se": bool(zero_se[j]),
            "simultaneous_ci95": (lo, hi),
            "excludes_zero": excludes,
        })
    individual = []
    for j in range(n_comp):
        column = draws[finite[:, j], j]
        if column.size == 0:
            individual.append((None, None))
        else:
            lo, hi = np.quantile(column, [(1 - level) / 2, 1 - (1 - level) / 2], method=QUANTILE_METHOD)
            individual.append((float(lo), float(hi)))
    return {
        "level": level,
        "quantile_method": QUANTILE_METHOD,
        "critical_value": c,
        "n_boot": int(n_boot),
        "n_incomplete": n_incomplete,
        "n_complete": int(complete.sum()),
        "interval_validity": validity if c is not None or validity == "unresolved" else "unresolved",
        "intervals": intervals,
        "individual_percentile_ci95": individual,
        "centring": "original_estimate",
    }
