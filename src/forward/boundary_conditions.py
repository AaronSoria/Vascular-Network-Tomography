"""
Junction Boundary Conditions
==============================
Computes reflection and transmission coefficients at arterial bifurcations
using linearised transmission-line theory.

Physical basis
--------------
At a junction with one incoming edge (e_in) and N outgoing edges {e_out_i},
the linearised wave equations give:

  Pressure continuity:   p_in + p_r = p_t  (same junction pressure for all branches)
  Flow conservation:     Y_in*(p_in - p_r) = Σ_i Y_out_i * p_t

Solving:
  Γ       = (Y_in - Σ Y_out_i) / (Y_in + Σ Y_out_i)     [reflection, for pressure]
  τ_total = 2 * Y_in / (Y_in + Σ Y_out_i)               [transmission, same for all branches]

Note: τ is the transmitted *pressure* amplitude.  All outgoing branches
receive the same junction pressure.  Flow is split proportionally to their
admittances.

For a *stenotic* edge, we model the stenosis as a discontinuity in cross-
section within the edge.  The reflection at the stenosis throat is computed
with the same formula, treating the throat as a junction between two segments
of the same edge with different impedances (one branch each side).

References
----------
- Parker & Jones (2009). Forward and backward waves in the arterial system.
  Med. Biol. Eng. Comput. 47(2):107-110.
- Olufsen MS et al. (2000). Numerical simulation and experimental validation of
  blood flow in arteries with structured-tree outflow conditions.
  Ann. Biomed. Eng. 28:1281-1299.
- Impedance matching at arterial bifurcations (1993). J. Biomech. 26:599-606.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import List

from src.graph.graph import VascularEdge, VascularGraph, VascularNode


@dataclass
class JunctionCoefficients:
    """
    Reflection and transmission coefficients at a vascular junction.

    Attributes
    ----------
    gamma : float
        Pressure reflection coefficient [-1, 1].
        Positive → forward reflection (same sign as incident wave).
        Negative → inverted reflection.
        0       → perfectly matched junction (no reflection).
    tau   : float
        Pressure transmission coefficient [0, 2].
        For a matched junction: tau = 1.
        Note: tau = 1 + gamma (conservation identity).
    incoming_edge_id : str
        The edge carrying the incident wave.
    node_id : str
        The junction node.
    outgoing_edge_ids : list[str]
        Edges that receive the transmitted wave.
    """
    gamma:               float
    tau:                 float
    incoming_edge_id:    str
    node_id:             str
    outgoing_edge_ids:   List[str]

    def __post_init__(self):
        # Sanity: tau = 1 + gamma (within floating-point tolerance)
        assert abs(self.tau - (1.0 + self.gamma)) < 1e-9, (
            f"tau={self.tau}, gamma={self.gamma}: conservation violated."
        )

    def __repr__(self) -> str:
        return (
            f"Junction({self.node_id}, "
            f"in={self.incoming_edge_id}, "
            f"Γ={self.gamma:+.4f}, τ={self.tau:.4f})"
        )


def compute_junction_coefficients(
    graph: VascularGraph,
    node_id: str,
    incoming_edge_id: str,
) -> JunctionCoefficients:
    """
    Compute Γ and τ at a junction node for a wave arriving along `incoming_edge_id`.

    Parameters
    ----------
    graph : VascularGraph
    node_id : str
        The junction node where the wave arrives.
    incoming_edge_id : str
        Edge id of the arriving wave (determines which edge is "source").

    Returns
    -------
    JunctionCoefficients
    """
    node = graph.node(node_id)

    # --- Terminal node: use terminal_reflection directly ----------------
    if node.terminal:
        gamma = node.terminal_reflection
        return JunctionCoefficients(
            gamma=gamma,
            tau=1.0 + gamma,
            incoming_edge_id=incoming_edge_id,
            node_id=node_id,
            outgoing_edge_ids=[],
        )

    # --- General junction -----------------------------------------------
    incoming_edge = graph.edge(incoming_edge_id)
    Y_in = incoming_edge.admittance()

    # All other edges at this node are "outgoing" from the wave's perspective
    outgoing = [
        e for e, _ in graph.neighbours(node_id)
        if e.edge_id != incoming_edge_id
    ]
    Y_out_total = sum(e.admittance() for e in outgoing)

    gamma = (Y_in - Y_out_total) / (Y_in + Y_out_total)
    tau   = 2.0 * Y_in / (Y_in + Y_out_total)

    return JunctionCoefficients(
        gamma=gamma,
        tau=tau,
        incoming_edge_id=incoming_edge_id,
        node_id=node_id,
        outgoing_edge_ids=[e.edge_id for e in outgoing],
    )


def compute_stenosis_coefficients(edge: VascularEdge) -> tuple[float, float]:
    """
    Compute the reflection and transmission coefficients at the proximal
    face of a stenotic throat within a single edge.

    The stenosis is modelled as an abrupt change from healthy impedance Z_h
    to stenosed impedance Z_s (one-way junction, single outgoing branch).

    Returns
    -------
    gamma_stenosis : float
        Reflection coefficient at the proximal stenosis face.
    tau_stenosis   : float
        Transmission coefficient into the stenotic segment.
    """
    if edge.stenosis is None:
        return 0.0, 1.0

    Y_h = edge.admittance()
    Y_s = 1.0 / edge.stenosis_impedance()

    # Single outgoing branch: Y_out = Y_s
    gamma = (Y_h - Y_s) / (Y_h + Y_s)
    tau   = 2.0 * Y_h / (Y_h + Y_s)
    return gamma, tau
