"""
Identifiability Analysis — Numerical Jacobian ∂M/∂θ
=====================================================
Checks whether the observation vector M = [T_AL, A_AL, T_AR, A_AR, T_TDR, A_TDR]
is sufficiently sensitive to the parameters θ = [D_left, D_right] to allow
unique recovery of θ from M.

For each scenario we:
1. Build the forward map F(θ) using WaveSolver
2. Compute ∂M/∂θ numerically (finite differences)
3. Inspect rank and condition number of the Jacobian
4. Identify which observations drive sensitivity

Output:
  experiments/identifiability_report.csv
  experiments/figures/jacobian_heatmap.png

A rank-deficient or very ill-conditioned Jacobian signals non-identifiability.
"""

from __future__ import annotations

import csv
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from src.datasets.synthetic import make_y_bifurcation
from src.forward.wave_solver import WaveSolver

SOLVER_PARAMS = dict(
    source_amplitude=1.0,
    min_amplitude=1e-5,
    max_reflections=6,
    max_time=0.6,
    include_stenosis_reflections=True,
)

D_HEALTHY      = 0.010       # 10 mm
JUNCTION_TOF_S = 2.0 * 0.350 / 4.5       # ≈ 0.1556 s
TDR_LO = JUNCTION_TOF_S + 1e-3           # 156.56 ms
TDR_HI = JUNCTION_TOF_S + 20e-3          # 175.56 ms — daughter roundtrip limit
TDR_THRESHOLD  = TDR_LO                  # kept for compat

OBS_NAMES = ['T_AL', 'A_AL', 'T_AR', 'A_AR', 'T_TDR', 'A_TDR']
PARAM_NAMES = ['D_left', 'D_right']

FIGURES_DIR = os.path.join(os.path.dirname(__file__), 'figures')
os.makedirs(FIGURES_DIR, exist_ok=True)


def forward_map(theta: np.ndarray, stenosis_edge=None, stenosis_severity=0.0) -> np.ndarray:
    """
    Map θ = [D_left, D_right] → M = [T_AL, A_AL, T_AR, A_AR, T_TDR, A_TDR].

    D_left and D_right are the full diameters [m] of e_left and e_right.
    The stenosis severity is derived from D_left if stenosis_edge='e_left'.
    """
    d_left, d_right = float(theta[0]), float(theta[1])

    # Derive stenosis severity from diameter reduction
    if stenosis_edge == 'e_left':
        sev = max(0.0, min(0.999, 1.0 - d_left / D_HEALTHY))
        scen = make_y_bifurcation(stenosis_edge='e_left', stenosis_severity=sev)
    elif stenosis_edge == 'e_right':
        sev = max(0.0, min(0.999, 1.0 - d_right / D_HEALTHY))
        scen = make_y_bifurcation(stenosis_edge='e_right', stenosis_severity=sev)
    else:
        # healthy or bilateral — use direct severity from diameters
        sev_l = max(0.0, 1.0 - d_left  / D_HEALTHY)
        sev_r = max(0.0, 1.0 - d_right / D_HEALTHY)
        if sev_l > 0.01 and sev_r > 0.01:
            # bilateral: use left as primary (simplified)
            scen = make_y_bifurcation(stenosis_edge='e_left', stenosis_severity=sev_l)
        else:
            scen = make_y_bifurcation()

    solver = WaveSolver(scen.graph, scen.source_node, ['A', 'L', 'R'],
                        **SOLVER_PARAMS)
    records = solver.run()
    tof_data = solver.extract_tof_amplitude(records, threshold_fraction=1e-5)

    def first(sid):
        arr = tof_data.get(sid, [])
        return arr[0] if arr else (0.0, 0.0)

    t_al, a_al = first('L')
    t_ar, a_ar = first('R')

    arrivals_A = tof_data.get('A', [])
    sten = [(t0, A) for t0, A in arrivals_A if TDR_LO < t0 < TDR_HI]
    t_tdr, a_tdr = sten[0] if sten else (0.0, 0.0)

    return np.array([t_al, a_al, t_ar, a_ar, t_tdr, a_tdr])


def numerical_jacobian(theta_0: np.ndarray, forward_fn, eps: float = 2e-5) -> np.ndarray:
    """
    Compute Jacobian J[i,j] = ∂M_i / ∂θ_j via central differences.

    Returns shape (n_obs, n_params).
    """
    M_0 = forward_fn(theta_0)
    n_obs, n_params = len(M_0), len(theta_0)
    J = np.zeros((n_obs, n_params))

    for j in range(n_params):
        th_plus = theta_0.copy()
        th_minus = theta_0.copy()
        th_plus[j]  += eps
        th_minus[j] -= eps
        M_plus  = forward_fn(th_plus)
        M_minus = forward_fn(th_minus)
        J[:, j] = (M_plus - M_minus) / (2.0 * eps)

    return J


def analyze_scenario(name, stenosis_edge, stenosis_severity, theta_0):
    """Compute and return Jacobian diagnostics for one scenario."""
    import functools
    fwd = functools.partial(forward_map,
                            stenosis_edge=stenosis_edge,
                            stenosis_severity=stenosis_severity)
    J = numerical_jacobian(theta_0, fwd)

    rank = np.linalg.matrix_rank(J, tol=1e-10)
    # Condition number of J (using SVD)
    sv = np.linalg.svd(J, compute_uv=False)
    cond = sv[0] / (sv[-1] if sv[-1] > 1e-15 else 1e-15)

    return J, rank, cond, sv


def main():
    D_HEALTHY = 0.010

    scenarios = [
        ('healthy',           None,     0.00),
        ('left_mild_30',     'e_left',  0.30),
        ('left_moderate_50', 'e_left',  0.50),
        ('left_severe_70',   'e_left',  0.70),
        ('right_severe_70',  'e_right', 0.70),
    ]

    print("Identifiability Analysis — Numerical Jacobian ∂M/∂θ")
    print(f"  θ = [D_left, D_right],  θ_0 = [{D_HEALTHY*1000:.0f}mm, {D_HEALTHY*1000:.0f}mm] (healthy)")
    print(f"  M = {OBS_NAMES}")
    print()

    report_rows = []
    all_J = {}

    header = (f"  {'Scenario':<22} | {'rank':>5} | {'cond_J':>12} | "
              f"  ∂M/∂D_left [top 3 sensitivities]")
    print(header)
    print("  " + "-" * 80)

    for name, edge, sev in scenarios:
        theta_0 = np.array([D_HEALTHY, D_HEALTHY])
        J, rank, cond, sv = analyze_scenario(name, edge, sev, theta_0)
        all_J[name] = J

        # Identify which observations have largest sensitivity to D_left
        sens_left = np.abs(J[:, 0])
        top_obs_idx = np.argsort(sens_left)[::-1][:3]
        top_obs_str = ', '.join(
            f'{OBS_NAMES[i]}={J[i,0]:.4g}' for i in top_obs_idx
        )

        print(f"  {name:<22} | {rank:>5} | {cond:>12.3e} |  {top_obs_str}")

        report_rows.append({
            'scenario': name,
            'rank': rank,
            'cond_J': cond,
            'sv_max': sv[0],
            'sv_min': sv[-1],
            'J_dAL_dDleft':  J[1, 0],
            'J_dTDR_dDleft': J[4, 0],
            'J_dAL_dDright': J[1, 1],
        })

    # Save CSV
    csv_path = os.path.join(os.path.dirname(__file__), 'identifiability_report.csv')
    with open(csv_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(report_rows[0].keys()))
        writer.writeheader()
        writer.writerows(report_rows)
    print(f"\n  Saved → experiments/identifiability_report.csv")

    # Heatmap figure — Jacobian for each scenario
    fig, axes = plt.subplots(1, len(scenarios), figsize=(14, 5), sharey=True)
    fig.suptitle('Jacobian ∂M/∂θ — Normalized Sensitivity\n'
                 'Y Bifurcation (θ = [D_left, D_right] at healthy baseline)',
                 fontsize=12)

    for ax, (name, _, __) in zip(axes, scenarios):
        J = all_J[name]
        # Normalize by max abs value for visualization
        J_norm = J / (np.max(np.abs(J)) + 1e-15)
        im = ax.imshow(J_norm, aspect='auto', cmap='RdBu', vmin=-1, vmax=1)
        ax.set_xticks([0, 1])
        ax.set_xticklabels(['D_left', 'D_right'], fontsize=8)
        ax.set_title(name.replace('_', '\n'), fontsize=8)
        if ax == axes[0]:
            ax.set_yticks(range(len(OBS_NAMES)))
            ax.set_yticklabels(OBS_NAMES, fontsize=8)

    fig.colorbar(im, ax=axes, label='Normalized sensitivity', fraction=0.02)
    fig.savefig(os.path.join(FIGURES_DIR, 'jacobian_heatmap.png'), dpi=150,
                bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved → experiments/figures/jacobian_heatmap.png")

    print()
    print("  Interpretation:")
    for row in report_rows:
        identifiable = "IDENTIFIABLE" if row['rank'] == 2 else "DEGENERATE"
        print(f"    {row['scenario']:<22}: rank={row['rank']} → {identifiable}  "
              f"(cond={row['cond_J']:.2e})")

    return report_rows


if __name__ == '__main__':
    main()
