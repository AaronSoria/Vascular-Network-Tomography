"""
Experiment: Forward Simulation
================================
Runs the forward wave solver on three scenarios and plots:
  1. Single tube  — healthy vs. 70% stenosis
  2. Y bifurcation — stenosis on left daughter
  3. 5-generation tree — femoral stenosis
  4. Observability analysis — sensitivity of observables to stenosis severity

Output: plots saved to experiments/figures/
"""

from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from src.datasets.synthetic import (
    make_single_tube,
    make_y_bifurcation,
    make_arterial_tree_5,
)
from src.forward.wave_solver import WaveSolver, gaussian_pulse

FIGURES_DIR = os.path.join(os.path.dirname(__file__), 'figures')
os.makedirs(FIGURES_DIR, exist_ok=True)

T_MAX     = 0.6
N_SAMPLES = 6000
SIGMA     = 2e-3
SOLVER_PARAMS = dict(
    source_amplitude=1.0,
    min_amplitude=5e-4,
    max_reflections=6,
    max_time=T_MAX,
    include_stenosis_reflections=True,
)

t  = np.linspace(0, T_MAX, N_SAMPLES)
pf = gaussian_pulse(sigma=SIGMA, amplitude=1.0)


# ===========================================================================
# Scenario 1 — Single tube
# ===========================================================================
def run_single_tube():
    print("\n=== Scenario 1: Single Tube ===")

    fig, axes = plt.subplots(2, 2, figsize=(12, 7))
    fig.suptitle('Single Tube — Healthy vs. 70% Stenosis', fontsize=13)

    for col, sev in enumerate([0.0, 0.70]):
        scen = make_single_tube(stenosis_severity=sev)
        solver = WaveSolver(scen.graph, scen.source_node,
                            scen.sensor_nodes, **SOLVER_PARAMS)
        records = solver.run()
        wf = solver.waveforms(records, t, pf)

        for row, sid in enumerate(['A', 'B']):
            ax = axes[row][col]
            color = 'steelblue' if sid == 'A' else 'darkorange'
            ax.plot(t * 1e3, wf[sid], lw=1.2, color=color)
            if row == 0:
                ax.set_title(f'{"Healthy" if sev == 0 else f"{sev*100:.0f}% stenosis"}')
            ax.set_ylabel(f'Sensor {sid} [Pa]')
            ax.set_xlabel('Time [ms]')
            ax.set_xlim(0, T_MAX * 1e3)
            ax.grid(alpha=0.3)

        tof = solver.extract_tof_amplitude(records)
        for sensor_id, arrivals in tof.items():
            label = 'Healthy' if sev == 0 else 'Stenosed'
            print(f"  {label:8s} | {sensor_id}: "
                  + ', '.join(f't={a[0]*1e3:.2f}ms A={a[1]:.4f}Pa'
                              for a in arrivals[:3]))

    plt.tight_layout()
    fig.savefig(os.path.join(FIGURES_DIR, 'single_tube.png'), dpi=150)
    plt.close(fig)
    print(f"  -> figures/single_tube.png")


# ===========================================================================
# Scenario 2 — Y bifurcation
# ===========================================================================
def run_y_bifurcation():
    print("\n=== Scenario 2: Y Bifurcation ===")

    configs = [
        (None,      0.0,  'Healthy'),
        ('e_left',  0.60, '60% stenosis — left'),
        ('e_right', 0.60, '60% stenosis — right'),
    ]
    sensor_colors = {'A': 'steelblue', 'L': 'darkorange', 'R': 'forestgreen'}

    fig, axes = plt.subplots(3, 3, figsize=(14, 9))
    fig.suptitle('Y Bifurcation — Stenosis Location Comparison', fontsize=13)

    for col, (eid, sev, title) in enumerate(configs):
        scen = make_y_bifurcation(stenosis_edge=eid, stenosis_severity=sev)
        solver = WaveSolver(scen.graph, scen.source_node,
                            scen.sensor_nodes, **SOLVER_PARAMS)
        records = solver.run()
        wf = solver.waveforms(records, t, pf)

        for row, sid in enumerate(['A', 'L', 'R']):
            ax = axes[row][col]
            ax.plot(t * 1e3, wf[sid], lw=1.2, color=sensor_colors[sid])
            if row == 0:
                ax.set_title(title, fontsize=10)
            if col == 0:
                ax.set_ylabel(f'Sensor {sid}\n[Pa]')
            ax.set_xlabel('Time [ms]')
            ax.set_xlim(0, 350)
            ax.grid(alpha=0.3)

        tof = solver.extract_tof_amplitude(records)
        print(f"  {title}:")
        for sid, arrivals in tof.items():
            first = arrivals[0] if arrivals else (float('nan'), 0)
            print(f"    {sid}: first t={first[0]*1e3:.2f}ms, A={first[1]:.4f}Pa "
                  f"({len(arrivals)} arrivals)")

    plt.tight_layout()
    fig.savefig(os.path.join(FIGURES_DIR, 'y_bifurcation.png'), dpi=150)
    plt.close(fig)
    print(f"  -> figures/y_bifurcation.png")


# ===========================================================================
# Scenario 3 — 5-generation arterial tree
# ===========================================================================
def run_arterial_tree():
    print("\n=== Scenario 3: 5-Generation Arterial Tree ===")

    configs = [
        ([],                                 'Healthy'),
        ([('fem_L', 0.70)],                  'Left femoral 70%'),
        ([('fem_L', 0.70), ('pop_R', 0.50)], 'Bilateral stenosis'),
    ]
    sensor_ids = ['root', 'int_il_L', 'int_il_R',
                  'tib_L1', 'tib_L2', 'tib_R1', 'tib_R2']

    fig, axes = plt.subplots(len(configs), 7, figsize=(18, 8))
    fig.suptitle('5-Generation Arterial Tree — Forward Simulation', fontsize=13)

    for row, (sten_cfg, title) in enumerate(configs):
        scen = make_arterial_tree_5(stenosis_configs=sten_cfg)
        solver = WaveSolver(scen.graph, scen.source_node,
                            scen.sensor_nodes, **SOLVER_PARAMS)
        records = solver.run()
        wf = solver.waveforms(records, t, pf)
        tof = solver.extract_tof_amplitude(records)
        print(f"\n  {title}:")

        for col, sid in enumerate(sensor_ids):
            ax = axes[row][col]
            ax.plot(t * 1e3, wf.get(sid, np.zeros_like(t)),
                    lw=1.0, color='steelblue')
            if row == 0:
                ax.set_title(sid.replace('_', '\n'), fontsize=8)
            if col == 0:
                ax.set_ylabel(title[:18], fontsize=7)
            ax.set_xlim(0, 450)
            ax.set_xlabel('ms', fontsize=7)
            ax.tick_params(labelsize=6)
            ax.grid(alpha=0.3)

            arrivals = tof.get(sid, [])
            if arrivals:
                print(f"    {sid:12s}: t0={arrivals[0][0]*1e3:.1f}ms, "
                      f"A={arrivals[0][1]:.4f}Pa ({len(arrivals)} arrivals)")

    plt.tight_layout()
    fig.savefig(os.path.join(FIGURES_DIR, 'arterial_tree_5.png'), dpi=150)
    plt.close(fig)
    print(f"  -> figures/arterial_tree_5.png")


# ===========================================================================
# Observability analysis
# ===========================================================================
def run_observability_analysis():
    """
    Tracks three observables vs stenosis severity in the Y bifurcation.

    Observable 1 — TOF at distal sensor L (shifts slightly: higher wave speed
                   through stenotic throat shortens transit time).
    Observable 2 — Amplitude at distal sensor L (amplified by impedance mismatch;
                   strongest signal).
    Observable 3 — TDR echo amplitude at proximal sensor A (the stenosis-
                   reflected echo that bounces back from the throat; the key
                   network-tomography-style signal).
    """
    print("\n=== Observability Analysis ===")
    print("  Tracking 3 observables vs stenosis severity (e_left, Y bifurcation)")

    severities = np.linspace(0.0, 0.90, 19)
    tof_L_list       = []
    amp_L_list       = []
    echo_amp_A_list  = []

    # Healthy baseline
    scen_h = make_y_bifurcation()
    solver_h = WaveSolver(scen_h.graph, scen_h.source_node,
                          ['A', 'L', 'R'], **SOLVER_PARAMS)
    records_h = solver_h.run()
    tof_h = solver_h.extract_tof_amplitude(records_h, threshold_fraction=1e-4)
    amp_L_healthy = tof_h['L'][0][1] if tof_h.get('L') else 0.0
    tof_L_healthy = tof_h['L'][0][0] * 1e3 if tof_h.get('L') else 0.0

    # Expected arrival time of the junction echo at A: 2 * L_parent / c_parent
    junction_tof_ms = 2.0 * 0.350 / 4.5 * 1e3  # ≈ 155.6 ms

    for sev in severities:
        scen = make_y_bifurcation(stenosis_edge='e_left',
                                  stenosis_severity=float(sev))
        solver = WaveSolver(scen.graph, scen.source_node,
                            ['A', 'L', 'R'], **SOLVER_PARAMS)
        records = solver.run()
        tof_data = solver.extract_tof_amplitude(records, threshold_fraction=1e-4)

        arrivals_L = tof_data.get('L', [])
        tof_L_list.append(arrivals_L[0][0] * 1e3 if arrivals_L else float('nan'))
        amp_L_list.append(arrivals_L[0][1] if arrivals_L else 0.0)

        # Stenosis echo at A: arrives AFTER the junction echo
        # (junction echo at ~155.6ms; stenosis echo at ~155.6 + 2*pos*L_left/c_left)
        arrivals_A = tof_data.get('A', [])
        sten_echoes = [(t0, A) for t0, A in arrivals_A
                       if t0 * 1e3 > junction_tof_ms + 0.5]
        echo_amp_A_list.append(abs(sten_echoes[0][1]) if sten_echoes else 0.0)

    # --- Plot ---
    fig, axes = plt.subplots(1, 3, figsize=(14, 4))
    fig.suptitle(
        'Sensitivity of Observables to Stenosis Severity\n'
        'Y Bifurcation — Stenosis on left daughter (e_left)',
        fontsize=12,
    )

    ax = axes[0]
    ax.axhline(tof_L_healthy, color='gray', ls='--', lw=1, label='Healthy')
    ax.plot(severities * 100, tof_L_list, 'o-', color='steelblue',
            ms=4, label='Stenosed')
    ax.set_xlabel('Stenosis severity [%]')
    ax.set_ylabel('First arrival TOF at sensor L [ms]')
    ax.set_title('Obs 1: Distal TOF\n(minor signal, from c increase)')
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3)

    ax = axes[1]
    ax.axhline(amp_L_healthy, color='gray', ls='--', lw=1, label='Healthy')
    ax.plot(severities * 100, amp_L_list, 's-', color='darkorange',
            ms=4, label='Stenosed')
    ax.set_xlabel('Stenosis severity [%]')
    ax.set_ylabel('First arrival amplitude at sensor L [Pa]')
    ax.set_title('Obs 2: Distal Amplitude\n(dominant signal)')
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3)

    ax = axes[2]
    ax.plot(severities * 100, echo_amp_A_list, '^-', color='forestgreen', ms=4)
    ax.set_xlabel('Stenosis severity [%]')
    ax.set_ylabel('Stenosis echo amplitude at sensor A [Pa]')
    ax.set_title('Obs 3: TDR Echo at Source A\n(key network tomography signal)')
    ax.grid(alpha=0.3)

    plt.tight_layout()
    fig.savefig(os.path.join(FIGURES_DIR, 'observability_analysis.png'), dpi=150)
    plt.close(fig)
    print(f"  -> figures/observability_analysis.png")

    # Numerical summary
    print(f"\n  {'Severity':>8} | {'TOF-L [ms]':>12} | "
          f"{'Amp-L [Pa]':>12} | {'Echo-A [Pa]':>12}")
    print(f"  {'-'*52}")
    for sev, tof, amp, echo in zip(
        severities[::3], tof_L_list[::3],
        amp_L_list[::3], echo_amp_A_list[::3]
    ):
        flag = ' <-- detectable echo' if echo > 5e-4 else ''
        print(f"  {sev*100:>7.0f}% | {tof:>12.3f} | "
              f"{amp:>12.4f} | {echo:>12.5f}{flag}")


# ===========================================================================
# Entry point
# ===========================================================================
if __name__ == '__main__':
    print("Forward Wave Simulation — Vascular Network Tomography")
    print("=" * 54)
    run_single_tube()
    run_y_bifurcation()
    run_arterial_tree()
    run_observability_analysis()
    print("\nAll scenarios complete.")
