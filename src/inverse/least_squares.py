"""
Inverse Solver — Tikhonov-Regularised Least Squares
=====================================================
Recovers θ = [D_left, D_right] from M_obs = [T_AL, A_AL, T_AR, A_AR, T_TDR, A_TDR]
using gradient-free Nelder-Mead minimisation of:

  J(θ) = ||W · (F(θ) − M_obs)||²  +  λ ||θ − θ_prior||²

W is a diagonal weight matrix that normalises observations with different
units (times in s, amplitudes in Pa) to comparable scale.

Physical domain:
  D_left, D_right ∈ [1 mm, 25 mm]

We implement a simple bounded Nelder-Mead using only numpy to avoid scipy
dependency. For production use, scipy.optimize.minimize with 'L-BFGS-B'
is preferred — the interface is identical.

Usage
-----
  from src.inverse.least_squares import InverseSolver
  solver = InverseSolver(forward_fn, theta_prior, theta_bounds)
  result = solver.solve(M_obs)
  print(result.theta_estimated, result.n_iter, result.success)
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, List, Optional, Tuple

import numpy as np


# ---------------------------------------------------------------------------
# Minimal bounded Nelder-Mead (numpy only)
# ---------------------------------------------------------------------------

def _nelder_mead(
    func,
    x0:       np.ndarray,
    bounds:   List[Tuple[float, float]],
    max_iter: int   = 500,
    xatol:    float = 1e-9,
    fatol:    float = 1e-12,
) -> Tuple[np.ndarray, float, int, bool]:
    """
    Bounded Nelder-Mead simplex minimiser (pure numpy).

    Bounds are enforced by reflecting out-of-bound vertices back inside.
    Returns (x_best, f_best, n_eval, converged).
    """
    n = len(x0)
    lo = np.array([b[0] for b in bounds])
    hi = np.array([b[1] for b in bounds])

    def clip(x):
        return np.clip(x, lo, hi)

    # Initial simplex
    simplex = np.zeros((n + 1, n))
    simplex[0] = clip(x0)
    step = np.maximum((hi - lo) * 0.05, 1e-4)
    for i in range(n):
        v = simplex[0].copy()
        v[i] += step[i]
        simplex[i + 1] = clip(v)

    fvals = np.array([func(simplex[i]) for i in range(n + 1)])
    n_eval = n + 1

    alpha, gamma, rho, sigma = 1.0, 2.0, 0.5, 0.5

    for _ in range(max_iter):
        order = np.argsort(fvals)
        simplex = simplex[order]
        fvals   = fvals[order]

        # Convergence checks
        if (np.max(np.abs(simplex[1:] - simplex[0])) < xatol and
                np.max(np.abs(fvals[1:] - fvals[0])) < fatol):
            return simplex[0], fvals[0], n_eval, True

        centroid = simplex[:-1].mean(axis=0)

        # Reflection
        xr = clip(centroid + alpha * (centroid - simplex[-1]))
        fr = func(xr);  n_eval += 1

        if fr < fvals[0]:
            # Expansion
            xe = clip(centroid + gamma * (xr - centroid))
            fe = func(xe);  n_eval += 1
            if fe < fr:
                simplex[-1], fvals[-1] = xe, fe
            else:
                simplex[-1], fvals[-1] = xr, fr
        elif fr < fvals[-2]:
            simplex[-1], fvals[-1] = xr, fr
        else:
            # Contraction
            if fr < fvals[-1]:
                xc = clip(centroid + rho * (xr - centroid))
                fc = func(xc);  n_eval += 1
                if fc <= fr:
                    simplex[-1], fvals[-1] = xc, fc
                else:
                    # Shrink
                    for i in range(1, n + 1):
                        simplex[i] = clip(simplex[0] + sigma * (simplex[i] - simplex[0]))
                        fvals[i]   = func(simplex[i]);  n_eval += 1
            else:
                xc = clip(centroid + rho * (simplex[-1] - centroid))
                fc = func(xc);  n_eval += 1
                if fc < fvals[-1]:
                    simplex[-1], fvals[-1] = xc, fc
                else:
                    for i in range(1, n + 1):
                        simplex[i] = clip(simplex[0] + sigma * (simplex[i] - simplex[0]))
                        fvals[i]   = func(simplex[i]);  n_eval += 1

    return simplex[0], fvals[0], n_eval, False


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

@dataclass
class InverseResult:
    """Output of one inverse solve."""
    theta_estimated: np.ndarray
    theta_true:      Optional[np.ndarray]
    M_obs:           np.ndarray
    M_pred:          np.ndarray
    residual_norm:   float
    final_loss:      float
    n_iter:          int
    success:         bool
    message:         str

    def summary(self, param_names=None, obs_names=None) -> str:
        pnames = param_names or [f'θ_{i}' for i in range(len(self.theta_estimated))]
        lines = ['InverseResult:']
        for i, name in enumerate(pnames):
            truth_str = ''
            if self.theta_true is not None:
                truth_str = f'  (true: {self.theta_true[i]*1e3:.3f} mm)'
            lines.append(f'  {name} = {self.theta_estimated[i]*1e3:.3f} mm{truth_str}')
        lines.append(f'  Residual norm = {self.residual_norm:.4e}')
        lines.append(f'  Loss          = {self.final_loss:.4e}')
        lines.append(f'  Evaluations   = {self.n_iter}')
        lines.append(f'  Converged     = {self.success}  ({self.message})')
        return '\n'.join(lines)


class InverseSolver:
    """
    Tikhonov-regularised inverse solver using bounded Nelder-Mead.

    Parameters
    ----------
    forward_fn    : callable θ → M_pred
    theta_prior   : np.ndarray — prior / initial point [m]
    bounds        : list of (lo, hi) tuples [m]
    lambda_reg    : float — regularisation weight (default 1e-5)
    weights       : np.ndarray or None — diagonal of W; auto if None
    max_iter      : int
    """

    def __init__(
        self,
        forward_fn:   Callable[[np.ndarray], np.ndarray],
        theta_prior:  np.ndarray,
        bounds:       List[Tuple[float, float]],
        lambda_reg:   float = 1e-5,
        weights:      Optional[np.ndarray] = None,
        max_iter:     int = 500,
    ):
        self.forward_fn  = forward_fn
        self.theta_prior = np.asarray(theta_prior, dtype=float)
        self.bounds      = bounds
        self.lambda_reg  = lambda_reg
        self.weights     = weights
        self.max_iter    = max_iter

    def _auto_weights(self, M_obs: np.ndarray) -> np.ndarray:
        return np.where(np.abs(M_obs) > 1e-12, 1.0 / np.abs(M_obs), 1.0)

    def _make_starts(self, theta_init: np.ndarray) -> List[np.ndarray]:
        """
        Multi-start grid: try several initial points to avoid local minima.

        Generates 5 candidates: the user-supplied init plus 4 evenly-spaced
        points across each parameter's range (lo, 25%, 50%, 75%, hi).
        All candidates are clipped to bounds.
        """
        lo = np.array([b[0] for b in self.bounds])
        hi = np.array([b[1] for b in self.bounds])
        candidates = [theta_init]
        # For each param, try low/mid-low/mid/mid-high initializations
        for fraction in [0.1, 0.3, 0.5]:
            alt = lo + fraction * (hi - lo)
            candidates.append(np.clip(alt, lo, hi))
        return candidates

    def solve(
        self,
        M_obs:      np.ndarray,
        theta_init: Optional[np.ndarray] = None,
        theta_true: Optional[np.ndarray] = None,
    ) -> InverseResult:
        M_obs = np.asarray(M_obs, dtype=float)
        theta_init = (np.asarray(theta_init, dtype=float)
                      if theta_init is not None
                      else self.theta_prior.copy())
        W = (self.weights if self.weights is not None
             else self._auto_weights(M_obs))

        def objective(theta):
            M_pred = self.forward_fn(theta)
            resid  = W * (M_pred - M_obs)
            return float(np.dot(resid, resid)) + self.lambda_reg * float(
                np.dot(theta - self.theta_prior, theta - self.theta_prior)
            )

        # Multi-start: run from each candidate and keep best result
        starts = self._make_starts(theta_init)
        x_best, f_best, n_eval_total, converged_any = None, float('inf'), 0, False

        for x0 in starts:
            x, f, n, conv = _nelder_mead(
                objective, x0, self.bounds, max_iter=self.max_iter,
            )
            n_eval_total += n
            if f < f_best:
                x_best, f_best = x, f
                converged_any = conv

        M_pred = self.forward_fn(x_best)
        W2     = self._auto_weights(M_obs)
        resid  = W2 * (M_pred - M_obs)

        return InverseResult(
            theta_estimated = x_best,
            theta_true      = np.asarray(theta_true) if theta_true is not None else None,
            M_obs           = M_obs,
            M_pred          = M_pred,
            residual_norm   = float(np.linalg.norm(resid)),
            final_loss      = f_best,
            n_iter          = n_eval_total,
            success         = converged_any,
            message         = 'Converged' if converged_any else 'Max iterations reached',
        )
