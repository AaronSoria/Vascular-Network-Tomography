"""
Vascular Network Graph Model
=============================
Models the arterial system as a directed graph G = (V, E) where:
  - V = arterial bifurcations / endpoints
  - E = arterial segments (edges)

Each edge carries physical parameters that determine wave propagation.
A stenosis is represented as a local reduction in lumen diameter.

Physical basis:
  - Characteristic impedance: Z_e = rho * c_e / A_e  [Pa·s/m³]
  - Admittance:               Y_e = A_e / (rho * c_e) [m³/(Pa·s)]
  - Wave speed (Moens-Korteweg): c_e = sqrt(E_e * h_e / (2 * rho_blood * R_e))
    Simplified here as a direct parameter c_e [m/s].

References:
  - Westerhof et al. (2009). The arterial Windkessel.
  - Parker et al. (2009). Forward and backward waves in the arterial system.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

# ---------------------------------------------------------------------------
# Physical constants
# ---------------------------------------------------------------------------
RHO_BLOOD = 1060.0   # kg/m³ — blood density
MU_BLOOD  = 0.003    # Pa·s  — dynamic viscosity (Newtonian approximation)


@dataclass
class VascularNode:
    """
    Represents a point in the arterial network: a bifurcation, anastomosis,
    or terminal endpoint.

    Attributes
    ----------
    node_id : str
        Unique identifier (e.g. 'aorta_root', 'iliac_bifurcation').
    label : str
        Human-readable anatomical label.
    terminal : bool
        True if this is a leaf node (no outgoing branches toward periphery).
    terminal_reflection : float
        Reflection coefficient at a terminal node.
        0.0  → matched load (absorbed, no echo)  — default for Phase I
        1.0  → closed end (fully reflected, same sign)
        -1.0 → open end  (fully reflected, inverted sign)
    """
    node_id:              str
    label:                str = ""
    terminal:             bool = False
    terminal_reflection:  float = 0.0  # matched load

    def __repr__(self) -> str:
        t = " [TERMINAL]" if self.terminal else ""
        return f"VascularNode({self.node_id}{t})"


@dataclass
class Stenosis:
    """
    Describes a stenotic lesion on an edge.

    Attributes
    ----------
    severity : float
        Fractional reduction in diameter: 0.0 = healthy, 0.7 = 70% stenosis.
        Diameter after stenosis: D_stenosed = D_healthy * (1 - severity).
    position : float
        Normalised position along the edge [0, 1].  0 = proximal, 1 = distal.
    length_fraction : float
        Fractional length of the stenotic segment relative to total edge length.
    """
    severity:         float   # 0–1
    position:         float = 0.5
    length_fraction:  float = 0.1


@dataclass
class VascularEdge:
    """
    Represents one arterial segment.

    All edges are undirected at the graph level; wave packets carry their own
    direction flag (proximal→distal or distal→proximal).

    Key derived quantities
    ----------------------
    area()        : cross-sectional area  [m²]
    impedance()   : characteristic impedance Z = rho*c/A  [Pa·s/m³]
    admittance()  : Y = 1/Z                              [m³/(Pa·s)]
    transit_time(): T = L/c                               [s]

    Stenosis effect
    ---------------
    When a stenosis is present the edge is modelled as three sections in series:
      proximal healthy | stenotic throat | distal healthy
    Each section has its own impedance.  The effective transfer matrix of the
    composite edge is computed by boundary_conditions.py.
    """
    edge_id:    str
    node_a:     str          # proximal node id
    node_b:     str          # distal  node id

    # --- Geometric and elastic parameters ---
    length_m:   float        # L_e  [m]
    diameter_m: float        # D_e  [m]  healthy lumen diameter
    wave_speed: float        # c_e  [m/s] — treated as independent parameter

    # --- Optional stenosis ---
    stenosis: Optional[Stenosis] = None

    # --- Viscous attenuation (per metre) ---
    # alpha = 8*pi*mu / A  (linearised Womersley/Poiseuille friction)
    # Stored for amplitude attenuation: A_out = A_in * exp(-alpha * L)
    viscous_attenuation: Optional[float] = None  # [Np/m]; None → auto-computed

    def __post_init__(self):
        if self.viscous_attenuation is None:
            # Linearised Poiseuille friction coefficient [1/m]
            # Derivation: 8*pi*mu*c / (rho * c² * A) ≈ 8*pi*mu / (rho*c*A)
            A = self.area()
            self.viscous_attenuation = (8.0 * math.pi * MU_BLOOD) / (
                RHO_BLOOD * self.wave_speed * A
            )

    # ------------------------------------------------------------------
    # Derived physical quantities
    # ------------------------------------------------------------------
    def area(self) -> float:
        """Cross-sectional area of the healthy lumen [m²]."""
        return math.pi * (self.diameter_m / 2.0) ** 2

    def stenosed_diameter(self) -> float:
        """Effective minimum diameter considering stenosis [m]."""
        if self.stenosis is None:
            return self.diameter_m
        return self.diameter_m * (1.0 - self.stenosis.severity)

    def stenosed_area(self) -> float:
        """Cross-sectional area at stenosis throat [m²]."""
        D = self.stenosed_diameter()
        return math.pi * (D / 2.0) ** 2

    def impedance(self) -> float:
        """
        Characteristic impedance of the healthy segment [Pa·s/m³].
        Z = rho * c / A
        """
        return RHO_BLOOD * self.wave_speed / self.area()

    def admittance(self) -> float:
        """Y = 1/Z [m³/(Pa·s)]."""
        return 1.0 / self.impedance()

    def stenosis_impedance(self) -> float:
        """Characteristic impedance at the stenosis throat [Pa·s/m³]."""
        if self.stenosis is None:
            return self.impedance()
        # Wave speed increases in stiffened/narrowed segment.
        # Simplified: c_stenosed = c_healthy * sqrt(D_healthy / D_stenosed)
        # (from Moens-Korteweg with constant E*h, R ∝ D/2)
        ratio = self.diameter_m / self.stenosed_diameter()
        c_s = self.wave_speed * math.sqrt(ratio)
        A_s = self.stenosed_area()
        return RHO_BLOOD * c_s / A_s

    def transit_time(self) -> float:
        """Time for wave to traverse the full edge [s]. T = L/c."""
        return self.length_m / self.wave_speed

    def amplitude_factor(self) -> float:
        """
        Amplitude multiplicative factor for viscous attenuation along the edge.
        F = exp(-alpha * L).  Dimensionless.
        """
        return math.exp(-self.viscous_attenuation * self.length_m)

    def __repr__(self) -> str:
        s = f" [STENOSIS {self.stenosis.severity*100:.0f}%]" if self.stenosis else ""
        return (
            f"VascularEdge({self.edge_id}: {self.node_a}→{self.node_b}, "
            f"L={self.length_m*100:.1f}cm, D={self.diameter_m*1000:.1f}mm, "
            f"c={self.wave_speed:.1f}m/s{s})"
        )


class VascularGraph:
    """
    The complete vascular network as an undirected graph.

    Internally stores adjacency as: node_id → list of (edge, neighbour_id).
    Wave direction is handled at the solver level.

    Usage
    -----
    >>> g = VascularGraph()
    >>> g.add_node(VascularNode('root', label='Aortic root'))
    >>> g.add_node(VascularNode('iliac_L', label='Left iliac', terminal=True))
    >>> g.add_edge(VascularEdge('aorta', 'root', 'iliac_L', 0.35, 0.025, 5.0))
    """

    def __init__(self):
        self._nodes: Dict[str, VascularNode] = {}
        self._edges: Dict[str, VascularEdge] = {}
        # adjacency: node_id → list of (VascularEdge, neighbour_id)
        self._adj: Dict[str, List[Tuple[VascularEdge, str]]] = {}

    # ------------------------------------------------------------------
    # Graph construction
    # ------------------------------------------------------------------
    def add_node(self, node: VascularNode) -> None:
        if node.node_id in self._nodes:
            raise ValueError(f"Node '{node.node_id}' already exists.")
        self._nodes[node.node_id] = node
        self._adj[node.node_id] = []

    def add_edge(self, edge: VascularEdge) -> None:
        for nid in (edge.node_a, edge.node_b):
            if nid not in self._nodes:
                raise ValueError(f"Node '{nid}' not found. Add nodes before edges.")
        if edge.edge_id in self._edges:
            raise ValueError(f"Edge '{edge.edge_id}' already exists.")
        self._edges[edge.edge_id] = edge
        self._adj[edge.node_a].append((edge, edge.node_b))
        self._adj[edge.node_b].append((edge, edge.node_a))

    # ------------------------------------------------------------------
    # Accessors
    # ------------------------------------------------------------------
    def node(self, node_id: str) -> VascularNode:
        return self._nodes[node_id]

    def edge(self, edge_id: str) -> VascularEdge:
        return self._edges[edge_id]

    @property
    def nodes(self) -> Dict[str, VascularNode]:
        return dict(self._nodes)

    @property
    def edges(self) -> Dict[str, VascularEdge]:
        return dict(self._edges)

    def neighbours(self, node_id: str) -> List[Tuple[VascularEdge, str]]:
        """Returns list of (edge, neighbour_node_id) for a given node."""
        return list(self._adj[node_id])

    def edges_at_node(self, node_id: str) -> List[VascularEdge]:
        """All edges incident on a node."""
        return [e for e, _ in self._adj[node_id]]

    def degree(self, node_id: str) -> int:
        return len(self._adj[node_id])

    def terminal_nodes(self) -> List[VascularNode]:
        return [n for n in self._nodes.values() if n.terminal]

    # ------------------------------------------------------------------
    # Validation
    # ------------------------------------------------------------------
    def validate(self) -> None:
        """Basic sanity checks. Raises ValueError on failure."""
        for nid, adj in self._adj.items():
            node = self._nodes[nid]
            if node.terminal and len(adj) != 1:
                raise ValueError(
                    f"Terminal node '{nid}' must have exactly 1 incident edge, "
                    f"found {len(adj)}."
                )
            if not node.terminal and len(adj) < 2:
                raise ValueError(
                    f"Non-terminal node '{nid}' must have at least 2 incident "
                    f"edges (is it a dangling node?), found {len(adj)}."
                )
        for eid, e in self._edges.items():
            if e.length_m <= 0:
                raise ValueError(f"Edge '{eid}' has non-positive length.")
            if e.diameter_m <= 0:
                raise ValueError(f"Edge '{eid}' has non-positive diameter.")
            if e.wave_speed <= 0:
                raise ValueError(f"Edge '{eid}' has non-positive wave speed.")
            if e.stenosis is not None:
                s = e.stenosis.severity
                if not (0.0 <= s < 1.0):
                    raise ValueError(
                        f"Edge '{eid}' stenosis severity {s} out of range [0, 1)."
                    )

    # ------------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------------
    def summary(self) -> str:
        lines = [
            f"VascularGraph: {len(self._nodes)} nodes, {len(self._edges)} edges",
        ]
        for nid, node in self._nodes.items():
            degree = self.degree(nid)
            t = " [terminal]" if node.terminal else ""
            lines.append(f"  {nid} (deg={degree}){t}")
        for eid, edge in self._edges.items():
            lines.append(f"  {edge}")
        return "\n".join(lines)
