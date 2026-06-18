"""
Dataset Generation — Minimum Validation Experiment
====================================================
Generates M_obs for 6 scenarios on the Y-bifurcation topology.

For each scenario, the forward solver runs and we extract the following
observation vector:
  M = [T_AL, A_AL, T_AR, A_AR, T_TDR, A_TDR]

where:
  T_AL  = TOF of first arrival at sensor L (from source A via e_parent → J → e_left)
  A_AL  = Amplitude of first arrival at sensor L
  T_AR  = TOF of first arrival at sensor R
  A_AR  = Amplitude of first arrival at sensor R
  T_TDR = TOF of stenosis echo at sensor A (0.0 if not present)
  A_TDR = Amplitude of stenosis echo at sensor A (0.0 if not present)

Dataset saved to: experiments/dataset_y_bifurcation.json
"""

from __future__ import annotations

import json
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))

from src.datasets.synthetic import make_y_bifurcation
from src.forward.wave_solver import WaveSolver

SOLVER_PARAMS = dict(
    source_amplitude=1.0,
    min_amplitude=1e-5,
    max_reflections=6,
    max_time=0.6,
    include_stenosis_reflections=True,
)

# Expected junction echo TOF at A: 2 * L_parent / c_parent = 2 * 0.350 / 4.5
JUNCTION_TOF_MS = 2.0 * 0.350 / 4.5 * 1e3   # ≈ 155.56 ms
# TDR echo can only arrive in (junction_tof + 1ms, junction_tof + daughter_roundtrip)
# daughter roundtrip = 2 * L_daughter / c_daughter = 2 * 0.070 / 7.0 = 20ms
TDR_LO_MS = JUNCTION_TOF_MS + 1.0           # 156.56 ms
TDR_HI_MS = JUNCTION_TOF_MS + 20.0          # 175.56 ms


def extract_observations(records, solver) -> dict:
    """
    Extract the 6-element observation vector from solver records.

    Returns dict with keys: T_AL, A_AL, T_AR, A_AR, T_TDR, A_TDR
    (times in seconds, amplitudes in Pa)
    """
    tof_data = solver.extract_tof_amplitude(records, threshold_fraction=1e-5)

    def first_arrival(sensor_id):
        arrivals = tof_data.get(sensor_id, [])
        if arrivals:
            return arrivals[0]   # (time_s, amplitude_Pa)
        return (float('nan'), 0.0)

    t_al, a_al = first_arrival('L')
    t_ar, a_ar = first_arrival('R')

    # TDR echo at A: arrives in (junction_tof+1ms, junction_tof+20ms)
    # Beyond 20ms window = 2nd junction roundtrip, not a stenosis echo.
    arrivals_A = tof_data.get('A', [])
    sten_echoes = [(t0, A) for t0, A in arrivals_A
                   if TDR_LO_MS * 1e-3 < t0 < TDR_HI_MS * 1e-3]
    if sten_echoes:
        t_tdr, a_tdr = sten_echoes[0]
    else:
        t_tdr, a_tdr = 0.0, 0.0

    return {
        'T_AL':  t_al,
        'A_AL':  a_al,
        'T_AR':  t_ar,
        'A_AR':  a_ar,
        'T_TDR': t_tdr,
        'A_TDR': a_tdr,
    }


def run_scenario(name, stenosis_edge=None, stenosis_severity=0.0,
                 d_left_true=None, d_right_true=None, severity_true=0.0):
    """Run one scenario and return the dataset record."""
    scen = make_y_bifurcation(
        stenosis_edge=stenosis_edge,
        stenosis_severity=stenosis_severity,
    )
    solver = WaveSolver(
        scen.graph, scen.source_node,
        ['A', 'L', 'R'],
        **SOLVER_PARAMS,
    )
    records = solver.run()
    obs = extract_observations(records, solver)

    D_HEALTHY = 0.010   # 10 mm baseline

    record = {
        'scenario': name,
        'ground_truth': {
            'stenosis_edge':     stenosis_edge,
            'stenosis_severity': stenosis_severity,
            'D_left_m':  d_left_true  if d_left_true  else D_HEALTHY,
            'D_right_m': d_right_true if d_right_true else D_HEALTHY,
        },
        'observations': obs,
    }
    return record


def main():
    D_HEALTHY = 0.010

    scenarios = [
        # name,              edge,      sev,  d_left (true),       d_right (true)
        ('healthy',           None,     0.00, D_HEALTHY,           D_HEALTHY),
        ('left_mild_30',     'e_left',  0.30, D_HEALTHY * 0.70,   D_HEALTHY),
        ('left_moderate_50', 'e_left',  0.50, D_HEALTHY * 0.50,   D_HEALTHY),
        ('left_severe_70',   'e_left',  0.70, D_HEALTHY * 0.30,   D_HEALTHY),
        ('right_severe_70',  'e_right', 0.70, D_HEALTHY,          D_HEALTHY * 0.30),
        ('bilateral_50',     'e_left',  0.50, D_HEALTHY * 0.50,   D_HEALTHY * 0.50),
    ]

    dataset = []
    print("Generating Y-bifurcation dataset...")
    print(f"  {'Scenario':<22} | {'T_AL [ms]':>10} | {'A_AL [Pa]':>10} | "
          f"{'T_TDR [ms]':>11} | {'A_TDR [Pa]':>11}")
    print("  " + "-" * 73)

    for name, edge, sev, d_left, d_right in scenarios:
        rec = run_scenario(name, stenosis_edge=edge, stenosis_severity=sev,
                           d_left_true=d_left, d_right_true=d_right)
        obs = rec['observations']
        dataset.append(rec)

        t_tdr_ms = obs['T_TDR'] * 1e3 if obs['T_TDR'] else 0.0
        print(f"  {name:<22} | {obs['T_AL']*1e3:>10.3f} | {obs['A_AL']:>10.5f} | "
              f"{t_tdr_ms:>11.3f} | {obs['A_TDR']:>11.6f}")

    # Save
    out_path = os.path.join(os.path.dirname(__file__), 'dataset_y_bifurcation.json')
    with open(out_path, 'w') as f:
        json.dump(dataset, f, indent=2)
    print(f"\n  Saved {len(dataset)} scenarios -> experiments/dataset_y_bifurcation.json")
    return dataset


if __name__ == '__main__':
    main()

