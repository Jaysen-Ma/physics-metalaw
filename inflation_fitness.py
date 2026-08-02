"""Fitness of single-field slow-roll inflationary potentials.

Reimplementation of the fitness value defined in Section II C of

    S. Abel, A. Constantin, T. R. Harvey, A. Lukas,
    "Cosmic Inflation and Genetic Algorithms", arXiv:2208.13804.

Everything is in reduced Planck units, M_p = 1/sqrt(8 pi G_N) = 1.

Pipeline for a potential V(phi)
-------------------------------
1.  Scan I_phi = [-30, 30] for local minima, take the lowest one, phi_m.
    Require phi_m in I_phi_m = [-5, 5], otherwise apply the no-vacuum penalty.
    Shift V by a constant so that V(phi_m) = 0.
2.  On each side of phi_m, the region of interest runs up to the nearest local
    maximum, or to the edge of I_phi if there is none.
3.  Walking outwards from phi_m, phi_e is the first field value where the
    slow-roll conditions eps <= 0.5 and |eta| <= 0.5 hold, with
        eps = (V'/V)^2 / 2,   eta = V''/V.
4.  Continuing outwards from phi_e, phi_i is the last field value where the
    slow-roll conditions *and* the no-eternal-inflation condition
        delta_phi_q = sqrt(V/3)/(2 pi) < 0.5 * delta_phi_c,
        delta_phi_c = 2 pi V' / (3 V)
    all still hold.
5.  Observables: N_total = |int_{phi_e}^{phi_i} V/V' dphi|,
    N_star = 58 + ln(V(phi_e))/6, phi_star fixed by
    |int_{phi_e}^{phi_star} V/V' dphi| = N_star, then
    n_s = 1 - 6 eps* + 2 eta*, r = 16 eps*, m = (V*/eps*)^(1/4).
6.  Fitness contributions (theta = Heaviside):
        f_N  = -(N_star - N_total)/10 * theta(N_star - N_total)
        f_ns = fh_ns * theta(f_term - fh_ns),
               fh_ns = -log10(1 + |n_s - 0.9649|/0.0042)
        f_r  = -log10(r/0.061) * theta(r - 0.061)
        f_m  = fh_m * theta(f_term - fh_m),
               fh_m  = -log10(1 + |m - 0.027|/0.0027)
    with f = f_N + f_ns + f_r + f_m and f_term = -0.3. A state is viable
    ("terminal") when f = 0.
7.  Both branches are scored and the larger fitness is returned.

Deviations from the paper, all deliberate:
  * The paper writes the slow-roll cut as "eta <= 0.5"; the absolute value
    |eta| <= 0.5 is used here (see ``FitnessConfig.use_abs_eta``).
  * The paper puts phi_i at the first field value where a condition *fails*;
    the last value where everything still holds is used instead, which avoids
    an integrable V/V' blow-up when the failure happens at a maximum. The two
    differ by one grid step.
  * The paper does not quote numerical values for its penalties, so the ones
    in ``FitnessConfig`` are a choice: they only need to rank a broken
    potential below any potential that at least produces a trajectory.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, replace
from typing import Callable, Optional, Sequence

import numpy as np

__all__ = [
    "Potential",
    "polynomial",
    "polynomial_plus_cosine",
    "polynomial_plus_exp_polynomial",
    "from_callable",
    "from_sympy",
    "FitnessConfig",
    "BranchReport",
    "FitnessReport",
    "evaluate",
    "fitness",
]

TWO_PI = 2.0 * math.pi


# --------------------------------------------------------------------------
# Potentials
# --------------------------------------------------------------------------


@dataclass(frozen=True)
class Potential:
    """A potential together with its first two derivatives.

    All three callables must accept and return numpy arrays. ``offset`` is the
    additive constant applied on top of ``v_raw``; it is what implements the
    paper's "adjust the additive constant such that V(phi_m) = 0" and it is not
    part of the model's state space.
    """

    v_raw: Callable[[np.ndarray], np.ndarray]
    dv: Callable[[np.ndarray], np.ndarray]
    d2v: Callable[[np.ndarray], np.ndarray]
    offset: float = 0.0
    label: str = ""

    def v(self, phi):
        return self.v_raw(phi) + self.offset

    def __call__(self, phi):
        return self.v(phi)

    def with_offset(self, offset: float) -> "Potential":
        return replace(self, offset=offset)

    def __str__(self) -> str:
        return self.label or "V(phi)"


def polynomial(coeffs: Sequence[float], label: str = "") -> Potential:
    """V(phi) = sum_k coeffs[k] phi^k, coefficients in ascending order."""
    c = np.asarray(coeffs, dtype=float)
    p = np.polynomial.Polynomial(c)
    p1, p2 = p.deriv(1), p.deriv(2)
    return Potential(p, p1, p2, 0.0, label or f"poly{tuple(c)}")


def polynomial_plus_cosine(
    poly_coeffs: Sequence[float], amplitude: float, frequency: float, label: str = ""
) -> Potential:
    """V(phi) = sum_k poly_coeffs[k] phi^k + amplitude * cos(frequency * phi)."""
    p = np.polynomial.Polynomial(np.asarray(poly_coeffs, dtype=float))
    p1, p2 = p.deriv(1), p.deriv(2)
    a, w = float(amplitude), float(frequency)

    def v(phi):
        return p(phi) + a * np.cos(w * phi)

    def dv(phi):
        return p1(phi) - a * w * np.sin(w * phi)

    def d2v(phi):
        return p2(phi) - a * w * w * np.cos(w * phi)

    return Potential(v, dv, d2v, 0.0, label or "poly+cos")


def polynomial_plus_exp_polynomial(
    poly_coeffs: Sequence[float],
    exp_rate: float,
    exp_poly_coeffs: Sequence[float],
    label: str = "",
) -> Potential:
    """V(phi) = P(phi) + Q(phi) * exp(exp_rate * phi)."""
    p = np.polynomial.Polynomial(np.asarray(poly_coeffs, dtype=float))
    q = np.polynomial.Polynomial(np.asarray(exp_poly_coeffs, dtype=float))
    p1, p2 = p.deriv(1), p.deriv(2)
    q1, q2 = q.deriv(1), q.deriv(2)
    k = float(exp_rate)

    def v(phi):
        return p(phi) + q(phi) * np.exp(k * phi)

    def dv(phi):
        return p1(phi) + (q1(phi) + k * q(phi)) * np.exp(k * phi)

    def d2v(phi):
        return p2(phi) + (q2(phi) + 2 * k * q1(phi) + k * k * q(phi)) * np.exp(k * phi)

    return Potential(v, dv, d2v, 0.0, label or "poly+exp*poly")


def from_callable(
    v: Callable[[np.ndarray], np.ndarray],
    dv: Optional[Callable[[np.ndarray], np.ndarray]] = None,
    d2v: Optional[Callable[[np.ndarray], np.ndarray]] = None,
    step: float = 1e-4,
    label: str = "",
) -> Potential:
    """Wrap an arbitrary callable, filling in missing derivatives numerically."""
    if dv is None:
        def dv(phi, _h=step):
            return (v(phi + _h) - v(phi - _h)) / (2 * _h)
    if d2v is None:
        def d2v(phi, _h=step):
            return (v(phi + _h) - 2 * v(phi) + v(phi - _h)) / (_h * _h)
    return Potential(v, dv, d2v, 0.0, label or "callable")


def from_sympy(expr, symbol=None, label: str = "") -> Potential:
    """Wrap a sympy expression, differentiating it symbolically."""
    import sympy as sy

    if symbol is None:
        free = tuple(expr.free_symbols)
        if len(free) != 1:
            raise ValueError(f"expected exactly one free symbol, got {free}")
        symbol = free[0]
    d1 = sy.diff(expr, symbol)
    d2 = sy.diff(d1, symbol)
    mods = ["numpy"]
    return Potential(
        sy.lambdify(symbol, expr, mods),
        sy.lambdify(symbol, d1, mods),
        sy.lambdify(symbol, d2, mods),
        0.0,
        label or str(expr),
    )


# --------------------------------------------------------------------------
# Configuration
# --------------------------------------------------------------------------


@dataclass(frozen=True)
class FitnessConfig:
    # Field ranges (Section II C).
    phi_min: float = -30.0          # I_phi
    phi_max: float = 30.0
    vacuum_min: float = -5.0        # I_phi_m, where the global minimum must sit
    vacuum_max: float = 5.0
    d_phi: float = 0.005            # scan step

    # Slow-roll and eternal-inflation cuts.
    eps_max: float = 0.5
    eta_max: float = 0.5
    use_abs_eta: bool = True
    qf_factor: float = 0.5          # delta_phi_q < qf_factor * delta_phi_c

    # Observables.
    n_star_offset: float = 58.0
    ns_obs: float = 0.9649
    ns_err: float = 0.0042
    r_max: float = 0.061
    m_obs: float = 0.027
    m_rel_err: float = 0.1          # Delta m = 0.1 m

    f_term: float = -0.3

    # Penalties for potentials that never reach the observable stage.
    no_vacuum_penalty: float = -100.0
    no_slow_roll_penalty: float = -50.0
    invalid_penalty: float = -100.0

    @property
    def m_err(self) -> float:
        return self.m_rel_err * self.m_obs


DEFAULT_CONFIG = FitnessConfig()


# --------------------------------------------------------------------------
# Reports
# --------------------------------------------------------------------------


@dataclass
class BranchReport:
    """Everything computed for one side of the global minimum."""

    direction: int                      # +1 right of phi_m, -1 left of it
    status: str                         # "ok" or a failure tag
    fitness: float
    region_end: float = math.nan        # nearest maximum, or edge of I_phi
    region_end_is_maximum: bool = False
    phi_e: float = math.nan
    phi_i: float = math.nan
    n_total: float = math.nan
    n_star: float = math.nan
    phi_star: float = math.nan
    phi_star_clamped: bool = False      # trajectory too short to reach phi_star
    eps_star: float = math.nan
    eta_star: float = math.nan
    v_star: float = math.nan
    ns: float = math.nan
    r: float = math.nan
    m: float = math.nan
    f_n: float = 0.0
    f_ns: float = 0.0
    f_r: float = 0.0
    f_m: float = 0.0

    @property
    def is_viable(self) -> bool:
        return self.status == "ok" and self.fitness == 0.0


@dataclass
class FitnessReport:
    fitness: float
    status: str
    label: str = ""
    phi_m: float = math.nan
    v_offset: float = 0.0               # constant added so that V(phi_m) = 0
    branches: tuple = ()
    best: Optional[BranchReport] = None

    @property
    def is_viable(self) -> bool:
        return self.status == "ok" and self.fitness == 0.0

    def summary(self) -> str:
        head = f"{self.label or 'V(phi)'}: fitness = {self.fitness:.4f} [{self.status}]"
        if self.best is None:
            return head
        b = self.best
        side = "right" if b.direction > 0 else "left"
        return (
            f"{head}\n"
            f"  phi_m      = {self.phi_m:+.4f}   (branch: {side})\n"
            f"  phi_i      = {b.phi_i:+.4f}   phi_e = {b.phi_e:+.4f}   "
            f"phi_* = {b.phi_star:+.4f}{' (clamped)' if b.phi_star_clamped else ''}\n"
            f"  N_total    = {b.n_total:.4g}   N_* = {b.n_star:.4g}\n"
            f"  n_s        = {b.ns:.4f}   r = {b.r:.4g}   m = {b.m:.4g}\n"
            f"  f_N/f_ns/f_r/f_m = {b.f_n:.3f} / {b.f_ns:.3f} / "
            f"{b.f_r:.3f} / {b.f_m:.3f}"
        )


# --------------------------------------------------------------------------
# Numerical helpers
# --------------------------------------------------------------------------


def _golden_min(f, a: float, b: float, tol: float = 1e-10, max_iter: int = 200):
    """Bounded golden-section minimisation of a scalar function."""
    invphi = (math.sqrt(5.0) - 1.0) / 2.0
    c, d = b - invphi * (b - a), a + invphi * (b - a)
    fc, fd = float(f(c)), float(f(d))
    for _ in range(max_iter):
        if b - a < tol:
            break
        if fc < fd:
            b, d, fd = d, c, fc
            c = b - invphi * (b - a)
            fc = float(f(c))
        else:
            a, c, fc = c, d, fd
            d = a + invphi * (b - a)
            fd = float(f(d))
    x = 0.5 * (a + b)
    return x, float(f(x))


def _bisect_root(f, a: float, b: float, tol: float = 1e-12, max_iter: int = 200):
    """Bisection for a sign change of f on [a, b]; returns the midpoint on failure."""
    fa, fb = float(f(a)), float(f(b))
    if not (math.isfinite(fa) and math.isfinite(fb)) or fa * fb > 0.0:
        return 0.5 * (a + b)
    for _ in range(max_iter):
        if b - a < tol:
            break
        mid = 0.5 * (a + b)
        fm = float(f(mid))
        if fm == 0.0:
            return mid
        if fa * fm < 0.0:
            b, fb = mid, fm
        else:
            a, fa = mid, fm
    return 0.5 * (a + b)


def _cumtrapz(y: np.ndarray, dx: float) -> np.ndarray:
    """Cumulative trapezoidal integral, same length as y, starting at 0."""
    out = np.zeros_like(y)
    if y.size > 1:
        np.cumsum(0.5 * (y[1:] + y[:-1]) * dx, out=out[1:])
    return out


def _scalar(f, x: float) -> float:
    return float(np.asarray(f(np.asarray(float(x)))))


# --------------------------------------------------------------------------
# Step 1: the global minimum
# --------------------------------------------------------------------------


def _find_global_minimum(pot: Potential, cfg: FitnessConfig):
    """Lowest local minimum of V inside I_phi. Returns (phi_m, V(phi_m)) or None."""
    n = int(round((cfg.phi_max - cfg.phi_min) / cfg.d_phi)) + 1
    grid = np.linspace(cfg.phi_min, cfg.phi_max, n)
    with np.errstate(all="ignore"):
        v = np.asarray(pot.v(grid), dtype=float)
    if not np.all(np.isfinite(v)):
        return None

    interior = v[1:-1]
    is_min = (interior < v[:-2]) & (interior <= v[2:])
    idx = np.flatnonzero(is_min) + 1
    if idx.size == 0:
        return None

    # Refine only the most promising candidates; the grid value already ranks them.
    order = idx[np.argsort(v[idx])][:16]
    best_phi, best_v = math.nan, math.inf
    for i in order:
        phi_c, v_c = _golden_min(lambda x: float(pot.v(np.asarray(x))),
                                 grid[i - 1], grid[i + 1])
        if v_c < best_v:
            best_phi, best_v = phi_c, v_c
    if not math.isfinite(best_v):
        return None
    return best_phi, best_v


# --------------------------------------------------------------------------
# Step 2-6: one branch
# --------------------------------------------------------------------------


def _region_end(pot: Potential, phi_m: float, direction: int, cfg: FitnessConfig):
    """Nearest local maximum away from phi_m, else the edge of I_phi."""
    boundary = cfg.phi_max if direction > 0 else cfg.phi_min
    span = abs(boundary - phi_m)
    n = max(int(span / cfg.d_phi) + 1, 2)
    s = np.linspace(0.0, span, n)
    phi = phi_m + direction * s
    with np.errstate(all="ignore"):
        slope = direction * np.asarray(pot.dv(phi), dtype=float)

    # V must rise away from the minimum; index 0 sits at the minimum itself.
    rising = slope > 0.0
    if n < 3 or not rising[1]:
        return phi_m, False
    turn = np.flatnonzero(~rising[1:]) + 1
    if turn.size == 0:
        return boundary, False
    j = int(turn[0])
    phi_max_pos = _bisect_root(
        lambda x: direction * float(pot.dv(np.asarray(x))), *sorted((phi[j - 1], phi[j]))
    )
    return phi_max_pos, True


def _branch(pot: Potential, phi_m: float, direction: int, cfg: FitnessConfig) -> BranchReport:
    end, end_is_max = _region_end(pot, phi_m, direction, cfg)

    def fail(tag: str) -> BranchReport:
        return BranchReport(direction, tag, cfg.no_slow_roll_penalty,
                            region_end=end, region_end_is_maximum=end_is_max)

    span = abs(end - phi_m)
    if span < 2 * cfg.d_phi:
        return fail("empty-region")

    n = int(span / cfg.d_phi) + 1
    s = np.linspace(0.0, span, n)
    dx = float(s[1] - s[0])
    phi = phi_m + direction * s

    with np.errstate(all="ignore"):
        v = np.asarray(pot.v(phi), dtype=float)
        d1 = np.asarray(pot.dv(phi), dtype=float)
        d2 = np.asarray(pot.d2v(phi), dtype=float)
        eps = 0.5 * (d1 / v) ** 2
        eta = d2 / v
        # delta_phi_q = sqrt(V/3)/(2 pi) < qf_factor * delta_phi_c,
        # delta_phi_c = 2 pi |V'| / (3 V)          [Eq. (4)]
        dphi_q = np.sqrt(np.maximum(v, 0.0) / 3.0) / TWO_PI
        dphi_c = TWO_PI * np.abs(d1) / (3.0 * v)

    finite = np.isfinite(v) & np.isfinite(d1) & np.isfinite(d2)
    eta_test = np.abs(eta) if cfg.use_abs_eta else eta
    slow_roll = finite & (v > 0.0) & (eps <= cfg.eps_max) & (eta_test <= cfg.eta_max)
    quantum_ok = finite & (dphi_q < cfg.qf_factor * dphi_c)

    # phi_e: first slow-roll point walking away from the minimum.
    candidates = np.flatnonzero(slow_roll[1:]) + 1
    if candidates.size == 0:
        return fail("no-slow-roll")
    i_e = int(candidates[0])

    # phi_i: walk on while slow roll and the no-eternal-inflation bound both hold.
    ok = slow_roll & quantum_ok
    if i_e + 1 >= n or not ok[i_e]:
        i_i = i_e
    else:
        broken = np.flatnonzero(~ok[i_e:])
        i_i = n - 1 if broken.size == 0 else i_e + int(broken[0]) - 1
    i_i = max(i_i, i_e)

    phi_e, phi_i = float(phi[i_e]), float(phi[i_i])
    v_e = float(v[i_e])
    if not (v_e > 0.0 and math.isfinite(v_e)):
        return fail("bad-v-at-phi-e")

    # e-folds measured from phi_e outwards, N(s) = |int V/V' dphi|.
    with np.errstate(all="ignore"):
        integrand = np.abs(v[i_e:] / d1[i_e:])
    integrand = np.where(np.isfinite(integrand), integrand, 0.0)
    n_cum = _cumtrapz(integrand, dx)

    n_total = float(n_cum[i_i - i_e])
    n_star = cfg.n_star_offset + math.log(v_e) / 6.0
    if not math.isfinite(n_star):
        return fail("bad-n-star")

    # phi_star sits N_star e-folds before the end of inflation, and cannot be
    # placed beyond where the slow-roll trajectory actually exists.
    reachable = n_cum[: i_i - i_e + 1]
    clamped = n_star > n_total
    if clamped:
        phi_star = phi_i
    else:
        s_star = float(np.interp(n_star, reachable, s[i_e : i_i + 1] - s[i_e]))
        phi_star = phi_e + direction * s_star

    with np.errstate(all="ignore"):
        v_s = _scalar(pot.v, phi_star)
        d1_s = _scalar(pot.dv, phi_star)
        d2_s = _scalar(pot.d2v, phi_star)
    if not (math.isfinite(v_s) and v_s > 0.0 and math.isfinite(d1_s) and math.isfinite(d2_s)):
        return fail("bad-pivot")

    eps_s = 0.5 * (d1_s / v_s) ** 2
    eta_s = d2_s / v_s
    if not math.isfinite(eps_s) or eps_s <= 0.0:
        return fail("bad-pivot")

    ns = 1.0 - 6.0 * eps_s + 2.0 * eta_s
    r = 16.0 * eps_s
    m = (v_s / eps_s) ** 0.25

    f_n = -0.1 * (n_star - n_total) if n_total < n_star else 0.0

    fh_ns = -math.log10(1.0 + abs(ns - cfg.ns_obs) / cfg.ns_err)
    f_ns = fh_ns if fh_ns < cfg.f_term else 0.0

    f_r = -math.log10(r / cfg.r_max) if r > cfg.r_max else 0.0

    fh_m = -math.log10(1.0 + abs(m - cfg.m_obs) / cfg.m_err)
    f_m = fh_m if fh_m < cfg.f_term else 0.0

    total = f_n + f_ns + f_r + f_m

    return BranchReport(
        direction=direction,
        status="ok",
        fitness=total,
        region_end=end,
        region_end_is_maximum=end_is_max,
        phi_e=phi_e,
        phi_i=phi_i,
        n_total=n_total,
        n_star=n_star,
        phi_star=phi_star,
        phi_star_clamped=clamped,
        eps_star=eps_s,
        eta_star=eta_s,
        v_star=v_s,
        ns=ns,
        r=r,
        m=m,
        f_n=f_n,
        f_ns=f_ns,
        f_r=f_r,
        f_m=f_m,
    )


# --------------------------------------------------------------------------
# Public entry points
# --------------------------------------------------------------------------


def evaluate(pot: Potential, cfg: FitnessConfig = DEFAULT_CONFIG) -> FitnessReport:
    """Score a potential and return the full diagnostic report."""
    label = str(pot)
    try:
        found = _find_global_minimum(pot, cfg)
    except (ArithmeticError, ValueError, TypeError):
        return FitnessReport(cfg.invalid_penalty, "evaluation-error", label)

    if found is None:
        return FitnessReport(cfg.no_vacuum_penalty, "no-minimum", label)

    phi_m, v_m = found
    if not (cfg.vacuum_min <= phi_m <= cfg.vacuum_max):
        return FitnessReport(cfg.no_vacuum_penalty, "minimum-outside-window", label,
                             phi_m=phi_m)

    shifted = pot.with_offset(pot.offset - v_m)

    try:
        branches = tuple(_branch(shifted, phi_m, d, cfg) for d in (+1, -1))
    except (ArithmeticError, ValueError, TypeError):
        return FitnessReport(cfg.invalid_penalty, "evaluation-error", label, phi_m=phi_m)

    best = max(branches, key=lambda b: b.fitness)
    return FitnessReport(
        fitness=best.fitness,
        status=best.status,
        label=label,
        phi_m=phi_m,
        v_offset=-v_m,
        branches=branches,
        best=best if best.status == "ok" else None,
    )


def fitness(pot: Potential, cfg: FitnessConfig = DEFAULT_CONFIG) -> float:
    """Scalar fitness, for use as a genetic-algorithm objective."""
    return evaluate(pot, cfg).fitness


if __name__ == "__main__":
    # The five worked examples of Section III, with the values quoted in the
    # paper for comparison.
    examples = [
        (
            "sextic poly, large field (eq. 14)",
            polynomial([3.57328e-10, 1.77313e-10, 1.77313e-11, -1.13662e-12,
                        1.45e-16, 1.95409e-15, -1.0e-16]),
            dict(phi_m=-3.7, phi_i=14.8, phi_e=0.0, phi_star=9.4,
                 n_total=2900, n_star=53.4, ns=0.9620, r=0.063, m=0.029),
        ),
        (
            "sextic poly, small field (eq. 15)",
            polynomial([2.32562e-11, 1.44128e-13, -1.33159e-13, 5.47169e-15,
                        -2.44753e-13, 8.63806e-15, 1.25252e-14]),
            dict(phi_m=-4.0, phi_i=0.5, phi_e=-2.2, phi_star=-0.2,
                 n_total=330, n_star=53.9, ns=0.9661, r=0.001, m=0.028),
        ),
        (
            "poly + cosine (eq. 17)",
            polynomial_plus_cosine(
                [2.80939e-11, 5.92241e-14, 1.13662e-11, -3.79642e-16, -4.10847e-15],
                amplitude=1.3662e-10, frequency=1.05,
            ),
            dict(phi_m=-2.5, phi_i=float("nan"), phi_e=-5.2, phi_star=-13.2,
                 n_total=265, n_star=54.4, ns=0.9611, r=0.048, m=0.028),
        ),
        (
            "poly + exp*poly, small field (eq. 19)",
            polynomial_plus_exp_polynomial(
                [5.26402e-11, 1.56e-14, 1e-15, 3.79642e-16, -2.1025e-16],
                exp_rate=1.2,
                exp_poly_coeffs=[-6.9349e-12, 0.0, -1.33159e-13, 1e-16, 1.68795e-13],
            ),
            dict(phi_m=2.1, phi_i=-3.8, phi_e=0.5, phi_star=-2.1,
                 n_total=1860, n_star=53.9, ns=0.9645, r=0.002, m=0.025),
        ),
        (
            "poly + exp*poly, large field (eq. 20)",
            polynomial_plus_exp_polynomial(
                [4.25907e-10, -1.33159e-12, 1.77313e-11, -7.28605e-13],
                exp_rate=-1.2,
                exp_poly_coeffs=[1.1641e-13, -1.05012e-9, 2.76609e-11,
                                 4.67054e-13, -3.04863e-16],
            ),
            dict(phi_m=0.8, phi_i=float("nan"), phi_e=float("nan"), phi_star=11.0,
                 n_total=743, n_star=54.4, ns=0.9613, r=0.058, m=0.025),
        ),
    ]

    row = "  {:<10} {:>14} {:>14}"
    for name, pot, paper in examples:
        report = evaluate(pot)
        report.label = name
        print(report.summary())
        b = report.best
        if b is not None:
            print(row.format("", "this code", "paper"))
            mine = dict(phi_m=report.phi_m, phi_i=b.phi_i, phi_e=b.phi_e,
                        phi_star=b.phi_star, n_total=b.n_total, n_star=b.n_star,
                        ns=b.ns, r=b.r, m=b.m)
            for key, ref in paper.items():
                print(row.format(key, f"{mine[key]:.4g}", f"{ref:.4g}"))
        print()
