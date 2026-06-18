# Vascular Network Tomography via Distributed Reflectometry

## Research Proposal v1

### Motivation

Current vascular diagnostic techniques such as catheterization, angiography, IVUS and OCT provide high-quality information but can be invasive, costly, or operationally complex.

This project investigates whether concepts from:

* Graph Theory
* Network Tomography
* Radar and Sonar Systems
* Time Domain Reflectometry (TDR)
* Distributed Systems
* Signal Processing
* Inverse Problems

can be combined to develop a new framework for detecting vascular stenosis through distributed sensing and reconstruction algorithms.

---

# 1. Mathematical Formulation

## Objective

Reconstruct the state of a vascular network from signals emitted and received by distributed sensing points.

The vascular system is modeled as a graph:

[
G=(V,E)
]

where:

* (V) = arterial bifurcations
* (E) = arterial segments

Each edge contains physical properties:

[
e=(L_e,D_e,Z_e)
]

where:

* (L_e) = length
* (D_e) = diameter
* (Z_e) = effective impedance

---

## Vascular State

Define:

[
\theta =
{
L_e,
D_e,
Z_e
}_{e\in E}
]

A stenosis modifies:

[
D_e \rightarrow D'_e
]

and therefore:

[
Z_e \rightarrow Z'_e
]

---

## Signal Propagation

A signal traverses a path:

[
P=(e_1,e_2,\ldots,e_n)
]

Total travel time:

[
T(P)=
\sum_{i=1}^{n}
\frac{L_i}{c_i}
]

where:

[
c_i
]

is the local propagation speed.

---

## Observation Matrix

Given a set of sensors:

[
S=
{
s_1,\ldots,s_k
}
]

we obtain:

[
M=
{
T_{ij},
A_{ij}
}
]

where:

* (T_{ij}) = time-of-flight
* (A_{ij}) = received amplitude

---

## Inverse Problem

Define:

[
F(\theta)=M
]

Goal:

[
\hat{\theta}
============

\arg\min_{\theta}
|
F(\theta)-M_{obs}
|^2
]

Interpretation:

Find the vascular configuration that best explains the observed measurements.

---

# 2. First Simulation Scope

The initial simulator intentionally ignores physiological complexity.

## Assumptions

* Static graph
* No cardiac pulse
* No vessel elasticity
* Constant propagation speed
* Simplified reflections
* Noise-free measurements

Goal:

Validate whether the inverse problem can be solved under ideal conditions.

---

# 3. Repository Architecture

```text
vascular-network-simulator/
│
├── cmd/
│   └── simulator/
│       └── main.go
│
├── internal/
│   ├── graph/
│   │   ├── node.go
│   │   ├── edge.go
│   │   └── graph.go
│   │
│   ├── signal/
│   │   ├── pulse.go
│   │   ├── propagation.go
│   │   └── reflection.go
│   │
│   ├── sensors/
│   │   ├── sensor.go
│   │   └── measurements.go
│   │
│   ├── stenosis/
│   │   └── stenosis.go
│   │
│   ├── inverse/
│   │   ├── solver.go
│   │   └── optimizer.go
│   │
│   └── visualization/
│       └── export.go
│
├── configs/
│   ├── healthy.json
│   ├── mild_stenosis.json
│   └── severe_stenosis.json
│
├── docs/
│   ├── research-proposal-v1.md
│   ├── math-model.md
│   ├── hardware-concept.md
│   └── roadmap.md
│
├── Dockerfile
├── docker-compose.yml
├── Makefile
└── README.md
```

---

# 4. Phase I — Forward Simulation

## Input

```json
{
  "nodes": [...],
  "edges": [...],
  "sensors": [...]
}
```

Process:

1. Build graph
2. Inject signal
3. Compute propagation paths
4. Compute reflections
5. Generate measurements

Output:

```json
{
  "sensorA": {
    "arrival_times": [...],
    "amplitudes": [...]
  }
}
```

---

# 5. Phase II — Reconstruction

The reconstruction algorithm receives only measurements.

Input:

```json
{
  "measurements": [...]
}
```

Expected output:

```json
{
  "suspected_stenosis": [
    {
      "edge": "e12",
      "confidence": 0.91
    }
  ]
}
```

---

# 6. Hardware Concepts

## Option A — Distributed Ultrasound

Components:

* Ultrasound transducer
* High-speed ADC
* Microcontroller
* Precision clock synchronization

Advantages:

* Existing technology
* Better understood physics
* Lower implementation risk

---

## Option B — Electromagnetic Reflectometry

Components:

* Miniaturized antenna
* UWB pulse generator
* Sensitive receiver
* DSP processing unit

Advantages:

* Closer to radar-inspired concept

Disadvantages:

* Strong attenuation in biological tissue
* More difficult signal interpretation

---

# 7. Scientific Metrics

## Localization Error

[
E_L
===

|x_{real}-x_{estimated}|
]

---

## Severity Error

[
E_S
===

|s_{real}-s_{estimated}|
]

---

## Coverage

[
Coverage
========

\frac{N_{detected}}
{N_{total}}
]

---

## Sensitivity

[
Sensitivity
===========

\frac{TP}
{TP+FN}
]

---

# 8. Research Hypothesis

A distributed sensing network based on reflectometry over arterial graphs can identify and localize vascular stenosis through network tomography techniques with clinically useful accuracy under simulated conditions.

---

# Next Steps

1. Literature review.
2. State-of-the-art analysis.
3. Graph generation engine.
4. Signal propagation engine.
5. Synthetic dataset generation.
6. Inverse problem solver.
7. Performance benchmarking.
8. Hardware feasibility study.
9. Comparison with IVUS, OCT and angiography.
10. Publication-quality experimental design.
