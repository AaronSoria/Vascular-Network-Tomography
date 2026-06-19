# Vascular Network Tomography via Distributed Reflectometry

**Research project** — Engineering / Biomedical / Signal Processing

## Hypothesis

A distributed sensing network based on reflectometry over arterial graphs can identify and localize vascular stenosis through network tomography techniques with clinically useful accuracy under simulated conditions.

---

## Repository Structure

```
src/
├── graph/             # VascularGraph, VascularNode, VascularEdge, Stenosis
├── forward/           # 1D wave solver (ray-tracing) + boundary conditions
├── sensors/           # Sensor recording
├── datasets/          # Synthetic scenario generators
├── metrics/           # Evaluation metrics (sensitivity, F1, RMSE)
└── inverse/           # Tikhonov-regularised Nelder-Mead inverse solver

experiments/
├── run_forward.py           # Phase I: forward simulation across all scenarios
├── generate_dataset.py      # Phase III: generate M_obs for 6 Y-bifurcation scenarios
├── identifiability_analysis.py  # Numerical Jacobian — rank and condition number
├── run_inverse.py           # Phase III: recover θ from M_obs, check success criteria
└── figures/                 # Output plots

simulations/
└── pressure_tdr_simulation.html  # Interactive P-TDR browser simulation (open in browser)

docs/
├── research-proposal-v1.md
├── literature-review-v1.md
├── validation-experiment-v1.md
├── hardware-version-A-cuff.md
├── hardware-version-B-ultrasound.md
└── hardware-comparative-analysis.md
```

---

## Physical Model

The vascular network is modelled as a graph G = (V, E):
- **Nodes V** = arterial bifurcations / endpoints
- **Edges E** = arterial segments, each with parameters (L, D, c)

Wave propagation uses the linearised 1D model (transmission line analogy):

| Quantity | Formula |
|---|---|
| Area | A = π(D/2)² |
| Impedance | Z = ρc/A |
| Admittance | Y = A/(ρc) |
| Transit time | T = L/c |

**Junction conditions** (pressure continuity + flow conservation):

```
Γ = (Y_parent − ΣY_daughters) / (Y_parent + ΣY_daughters)
τ = 2·Y_parent / (Y_parent + ΣY_daughters)
```

**Stenosis model**: lumen diameter reduced by severity s → D' = D(1-s).
This creates an internal impedance discontinuity that generates:
1. A **TDR echo** back toward the source (key diagnostic signal)
2. **Amplitude amplification** at distal sensors
3. A slight **TOF reduction** due to higher wave speed in the stiffened throat

---

## Key Results (Phase I)

From `experiments/run_forward.py`:

| Observable | Change with 70% stenosis | Diagnostic value |
|---|---|---|
| Distal TOF | -0.45 ms (subtle) | Low alone |
| Distal amplitude | +90% (1.63 → 3.10 Pa) | **High** |
| TDR echo at source | +0.25 Pa at t=164ms | **High (location-specific)** |
| Contralateral sensor | No change | Confirms asymmetric lesion |

---

## Quick Start

```bash
pip install numpy matplotlib scipy
# Phase I — forward simulation
python experiments/run_forward.py
# Phase III — validation experiment (inverse problem)
python experiments/generate_dataset.py
python experiments/identifiability_analysis.py
python experiments/run_inverse.py
```

Figures are saved to `experiments/figures/`.

## Interactive Simulation

Open `simulations/pressure_tdr_simulation.html` directly in any browser (no server needed).

Controls: stenosis side, severity (10–90%), wave speed c, chirp frequency range, excitation amplitude.
Displays: animated wave propagation in the arterial network, pressure signals at all three sensors (A/L/R), TDR impulse response with detection window, and computed metrics (TOF, echo time, stenosis distance, SNR).

## Validation Results (Phase III)

Minimum experiment: Y-bifurcation, topology known, no noise, 5 scenarios.

| Metric | Result | Threshold | Status |
|---|---|---|---|
| Localization accuracy | 100% | 100% | ✓ |
| Diameter error pass rate | 100% | ≥ 80% | ✓ |
| Severity error pass rate | 100% | ≥ 80% | ✓ |

**Identifiability finding:** Single-sided stenosis problems are rank-1 (one effective degree of freedom — the affected diameter). Full rank-2 identifiability (bilateral stenosis) is deferred to Phase IV.

**Inverse method:** Tikhonov-regularised Nelder-Mead with 4-start grid initialisation. Converges in 474–663 forward evaluations per scenario (~0.5s total).

---

## Phase Roadmap

| Phase | Status | Description |
|---|---|---|
| I — Forward Model | ✅ Complete | Ray-tracing solver, 3 scenario levels |
| II — Inverse Problem | ✅ Complete | Tikhonov-regularised Nelder-Mead (multi-start) |
| III — Validation | ✅ Complete | Synthetic dataset, identifiability, 100% criteria met |
| IV — Hardware Concept | ✅ Complete | Version A (cuff P-TDR) + Version B (guided wave), comparative analysis |
