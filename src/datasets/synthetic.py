"""
Synthetic Dataset Generator
============================
Builds vascular network scenarios for simulation and inverse problem validation.

Three canonical topologies are provided:
  1. single_tube      — one segment, two sensors (ground truth trivial)
  2. y_bifurcation    — parent + two daughters, 3 sensors
  3. arterial_tree_5  — a 5-generation tree (~31 edges) mimicking
                        a simplified femoral / coronary topology

Each scenario returns:
  - graph           : VascularGraph
  - source_node     : str
  - sensor_nodes    : list[str]
  - ground_truth    : dict mapping edge_id → true parameters
  - stenosis_edges  : list[str]  (empty if healthy scenario)

Physiological parameter ranges (Westerhof et al. 2010; Olufsen 1999):
  Aorta          : D = 25 mm,  c = 4.5 m/s,  L = 350 mm
  Iliac common   : D = 10 mm,  c = 7.0 m/s,  L = 70  mm
  Femoral        : D = 8  mm,  c = 8.0 m/s,  L = 250 mm
  Popliteal      : D = 6  mm,  c = 9.0 m/s,  L = 200 mm
  Tibial         : D = 3  mm,  c = 12 m/s,   L = 300 mm
  Carotid common : D = 7  mm,  c = 8.0 m/s,  L = 120 mm
  Coronary LAD   : D = 3  mm,  c = 12 m/s,   L = 80  mm
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

from src.graph.graph import Stenosis, VascularEdge, VascularGraph, VascularNode


@dataclass
class Scenario:
    """Container for one simulation scenario."""
    name:           str
    graph:          VascularGraph
    source_node:    str
    sensor_nodes:   List[str]
    ground_truth:   Dict[str, Dict]   # edge_id → {length_m, diameter_m, wave_speed}
    stenosis_edges: List[str]
    description:    str = ""


# ---------------------------------------------------------------------------
# Level 1 — Single tube
# ---------------------------------------------------------------------------

def make_single_tube(
    length_m:   float = 0.30,
    diameter_m: float = 0.010,
    wave_speed: float = 8.0,
    stenosis_severity: float = 0.0,
    stenosis_position: float = 0.5,
) -> Scenario:
    """
    One arterial segment with a sensor at each end.

      [sensor_A] ——— e0 ——— [sensor_B]

    Validation target: can the inverse recover D_e0 and L_e0?
    """
    g = VascularGraph()
    g.add_node(VascularNode('A', label='Proximal end', terminal=True))
    g.add_node(VascularNode('B', label='Distal end',   terminal=True))

    stenosis = None
    if stenosis_severity > 0:
        stenosis = Stenosis(
            severity=stenosis_severity,
            position=stenosis_position,
        )

    g.add_edge(VascularEdge(
        edge_id    = 'e0',
        node_a     = 'A',
        node_b     = 'B',
        length_m   = length_m,
        diameter_m = diameter_m,
        wave_speed = wave_speed,
        stenosis   = stenosis,
    ))

    return Scenario(
        name        = 'single_tube',
        graph       = g,
        source_node = 'A',
        sensor_nodes = ['A', 'B'],
        ground_truth = {
            'e0': dict(length_m=length_m, diameter_m=diameter_m, wave_speed=wave_speed)
        },
        stenosis_edges = ['e0'] if stenosis_severity > 0 else [],
        description = f'Single tube, D={diameter_m*1e3:.1f}mm, L={length_m*1e2:.1f}cm, '
                      f'stenosis={stenosis_severity*100:.0f}%',
    )


# ---------------------------------------------------------------------------
# Level 2 — Y bifurcation
# ---------------------------------------------------------------------------

def make_y_bifurcation(
    stenosis_edge:    Optional[str] = None,
    stenosis_severity: float = 0.0,
) -> Scenario:
    """
    Simple Y-shaped network.

      [source=A] ——— e_parent ——— [junction=J]
                                        ├—— e_left ——— [sensor_L]
                                        └—— e_right —— [sensor_R]

    Sensor at A and both terminals.
    Validation target: identify which daughter (or parent) has the stenosis.

    Parameters
    ----------
    stenosis_edge     : 'e_parent', 'e_left', 'e_right', or None
    stenosis_severity : 0–1
    """
    g = VascularGraph()
    g.add_node(VascularNode('A', label='Aortic root',     terminal=True,
                             terminal_reflection=0.3))  # partial cardiac reflection
    g.add_node(VascularNode('J', label='Aortic bifurcation'))
    g.add_node(VascularNode('L', label='Left iliac terminal',  terminal=True))
    g.add_node(VascularNode('R', label='Right iliac terminal', terminal=True))

    def _sten(eid):
        if eid == stenosis_edge and stenosis_severity > 0:
            return Stenosis(severity=stenosis_severity)
        return None

    g.add_edge(VascularEdge('e_parent', 'A', 'J',
                             length_m=0.350, diameter_m=0.025, wave_speed=4.5,
                             stenosis=_sten('e_parent')))
    g.add_edge(VascularEdge('e_left',   'J', 'L',
                             length_m=0.070, diameter_m=0.010, wave_speed=7.0,
                             stenosis=_sten('e_left')))
    g.add_edge(VascularEdge('e_right',  'J', 'R',
                             length_m=0.070, diameter_m=0.010, wave_speed=7.0,
                             stenosis=_sten('e_right')))

    gt = {
        'e_parent': dict(length_m=0.350, diameter_m=0.025, wave_speed=4.5),
        'e_left':   dict(length_m=0.070, diameter_m=0.010, wave_speed=7.0),
        'e_right':  dict(length_m=0.070, diameter_m=0.010, wave_speed=7.0),
    }

    return Scenario(
        name         = 'y_bifurcation',
        graph        = g,
        source_node  = 'A',
        sensor_nodes = ['A', 'L', 'R'],
        ground_truth = gt,
        stenosis_edges = [stenosis_edge] if stenosis_edge and stenosis_severity > 0 else [],
        description  = f'Y bifurcation — stenosis on {stenosis_edge} '
                       f'({stenosis_severity*100:.0f}%)',
    )


# ---------------------------------------------------------------------------
# Level 3 — Arterial tree (5 levels, simplified lower limb)
# ---------------------------------------------------------------------------

def make_arterial_tree_5(
    stenosis_configs: Optional[List[Tuple[str, float]]] = None,
) -> Scenario:
    """
    5-generation tree approximating a simplified lower-limb arterial network.

    Topology (proximal → distal):
                         [root]
                           |
                       [aorta_abdom]
                           |
                     [aortic_bifurc]
                      /           \\
              [iliac_L]         [iliac_R]
               /    \\             /    \\
          [fem_L] [int_il_L]  [fem_R] [int_il_R]  ← generation 3
            |                    |
         [pop_L]              [pop_R]              ← generation 4
          / \\                   / \\
    [tib_L1][tib_L2]     [tib_R1][tib_R2]         ← generation 5 (terminals)

    int_il = internal iliac (also terminal)

    Sensors at: root + all 6 terminal nodes (tib_L1, tib_L2, tib_R1, tib_R2,
                int_il_L, int_il_R)

    Parameters
    ----------
    stenosis_configs : list of (edge_id, severity), e.g. [('fem_L', 0.70)]
    """
    if stenosis_configs is None:
        stenosis_configs = []

    sten_map = {eid: sev for eid, sev in stenosis_configs}

    def _sten(eid):
        sev = sten_map.get(eid, 0.0)
        return Stenosis(severity=sev) if sev > 0 else None

    g = VascularGraph()

    # --- Nodes ---
    nodes = [
        ('root',          'Aortic root',                  True,  False),
        ('aortic_bifurc', 'Aortic bifurcation',            False, False),
        ('iliac_L_node',  'Left iliac bifurcation',         False, False),
        ('iliac_R_node',  'Right iliac bifurcation',        False, False),
        ('fem_L_dist',    'Left femoral distal',            False, False),
        ('fem_R_dist',    'Right femoral distal',           False, False),
        ('int_il_L',      'Left internal iliac terminal',   True,  False),
        ('int_il_R',      'Right internal iliac terminal',  True,  False),
        ('pop_L_dist',    'Left popliteal distal',          False, False),
        ('pop_R_dist',    'Right popliteal distal',         False, False),
        ('tib_L1',        'Left tibial 1 terminal',         True,  False),
        ('tib_L2',        'Left tibial 2 terminal',         True,  False),
        ('tib_R1',        'Right tibial 1 terminal',        True,  False),
        ('tib_R2',        'Right tibial 2 terminal',        True,  False),
    ]
    for nid, label, terminal, _ in nodes:
        g.add_node(VascularNode(nid, label=label, terminal=terminal,
                                terminal_reflection=0.0))

    # --- Edges ---
    # (edge_id, node_a, node_b, L[m], D[m], c[m/s])
    edges_def = [
        # Generation 0–1
        ('aorta_abdom', 'root',          'aortic_bifurc', 0.20,  0.020, 5.0),
        # Generation 1–2
        ('iliac_L',     'aortic_bifurc', 'iliac_L_node',  0.070, 0.010, 7.0),
        ('iliac_R',     'aortic_bifurc', 'iliac_R_node',  0.070, 0.010, 7.0),
        # Generation 2–3
        ('fem_L',       'iliac_L_node',  'fem_L_dist',    0.250, 0.008, 8.0),
        ('int_il_L',    'iliac_L_node',  'int_il_L',      0.050, 0.006, 9.0),
        ('fem_R',       'iliac_R_node',  'fem_R_dist',    0.250, 0.008, 8.0),
        ('int_il_R',    'iliac_R_node',  'int_il_R',      0.050, 0.006, 9.0),
        # Generation 3–4
        ('pop_L',       'fem_L_dist',    'pop_L_dist',    0.200, 0.006, 9.0),
        ('pop_R',       'fem_R_dist',    'pop_R_dist',    0.200, 0.006, 9.0),
        # Generation 4–5
        ('tib_L_ant',   'pop_L_dist',    'tib_L1',        0.300, 0.003, 12.0),
        ('tib_L_post',  'pop_L_dist',    'tib_L2',        0.300, 0.003, 12.0),
        ('tib_R_ant',   'pop_R_dist',    'tib_R1',        0.300, 0.003, 12.0),
        ('tib_R_post',  'pop_R_dist',    'tib_R2',        0.300, 0.003, 12.0),
    ]
    ground_truth = {}
    for eid, na, nb, L, D, c in edges_def:
        g.add_edge(VascularEdge(eid, na, nb,
                                length_m=L, diameter_m=D, wave_speed=c,
                                stenosis=_sten(eid)))
        ground_truth[eid] = dict(length_m=L, diameter_m=D, wave_speed=c)

    sensor_nodes = [
        'root', 'int_il_L', 'int_il_R',
        'tib_L1', 'tib_L2', 'tib_R1', 'tib_R2',
    ]

    return Scenario(
        name          = 'arterial_tree_5',
        graph         = g,
        source_node   = 'root',
        sensor_nodes  = sensor_nodes,
        ground_truth  = ground_truth,
        stenosis_edges = [eid for eid, _ in stenosis_configs],
        description   = (
            f'5-generation lower-limb tree, '
            f'{len(g.edges)} edges, {len(g.nodes)} nodes, '
            f'stenoses: {stenosis_configs}'
        ),
    )
