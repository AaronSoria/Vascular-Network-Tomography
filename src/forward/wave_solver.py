"""
Forward Wave Solver — Linearised 1D Ray-Tracing
=================================================
Simulates the propagation of a pressure wave pulse through a vascular network
using a ray-tracing (method of characteristics) approach.

Physical model
--------------
The linearised 1D wave equation for a compliant tube admits solutions as
superpositions of forward and backward traveling wave packets.  Each packet
carries a pressure amplitude; flow is derived from amplitude / impedance.

At every junction the boundary conditions (pressure continuity + flow
conservation) yield reflection coefficient Γ and transmission coefficient τ.
At a stenotic throat inside an edge, the same conditions produce an internal
echo — the key TDR-style diagnostic signal.

Algorithm
---------
1. Inject one packet per edge incident on the source node.
2. Use a min-heap (priority queue) ordered by arrival time.
3. At each node: record at sensors, spawn reflected + transmitted packets.
4. If an outgoing edge has a stenosis: also spawn the stenosis echo packet.
5. Prune packets below min_amplitude, exceeding max_reflections, or max_time.

References
----------
- Sherwin SJ et al. (2003). Computational modelling of 1D blood flow.
- Parker KH (2009). A brief history of wave intensity analysis.
- Impedance matching at arterial bifurcations (1993). J. Biomech.
"""

from __future__ import annotations

import heapq
import math
from dataclasses import dataclass, field
from typing import Callable, Dict, List, Optional, Tuple

import numpy as np

from src.forward.boundary_conditions import (
    compute_junction_coefficients,
    compute_stenosis_coefficients,
)
from src.graph.graph import VascularEdge, VascularGraph


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

@dataclass(order=True)
class WavePacket:
    """
    A single propagating wave packet.

    Ordered by arrival_time for min-heap scheduling.

    Attributes
    ----------
    arrival_time    : absolute time of arrival at next_node_id [s]
    amplitude       : pressure amplitude [Pa]
    next_node_id    : node where this packet will be processed
    current_edge_id : edge just traversed
    path            : history of (edge_id, tag) tuples
    n_reflections   : total reflection events accumulated
    """
    arrival_time:    float
    amplitude:       float               = field(compare=False)
    next_node_id:    str                 = field(compare=False)
    current_edge_id: str                 = field(compare=False)
    path:            List[Tuple[str, str]] = field(
        default_factory=list, compare=False, repr=False
    )
    n_reflections:   int                 = field(default=0, compare=False)


@dataclass
class SensorRecord:
    """Accumulates impulse arrivals at a sensor node."""
    node_id:  str
    arrivals: List[Tuple[float, float]] = field(default_factory=list)

    def record(self, time: float, amplitude: float) -> None:
        self.arrivals.append((time, amplitude))

    def waveform(
        self,
        t: np.ndarray,
        pulse_fn: Callable[[np.ndarray], np.ndarray],
    ) -> np.ndarray:
        """
        Convolve impulse arrivals with pulse_fn to produce the sensor waveform.

        Parameters
        ----------
        t        : time axis [s]
        pulse_fn : callable, called as pulse_fn(t - t0)

        Returns
        -------
        np.ndarray — pressure waveform [Pa]
        """
        p = np.zeros_like(t, dtype=float)
        for t0, A in self.arrivals:
            p += A * pulse_fn(t - t0)
        return p


# ---------------------------------------------------------------------------
# Source pulse shapes
# ---------------------------------------------------------------------------

def gaussian_pulse(sigma: float = 2e-3, amplitude: float = 1.0) -> Callable:
    """Gaussian pressure pulse centred at t=0. sigma in seconds."""
    def _fn(t: np.ndarray) -> np.ndarray:
        return amplitude * np.exp(-0.5 * (t / sigma) ** 2)
    return _fn


def raised_cosine_pulse(duration: float = 5e-3, amplitude: float = 1.0) -> Callable:
    """Raised cosine pulse of finite duration [0, duration]."""
    def _fn(t: np.ndarray) -> np.ndarray:
        out = np.zeros_like(t, dtype=float)
        mask = (t >= 0) & (t <= duration)
        out[mask] = amplitude * 0.5 * (1.0 - np.cos(2.0 * math.pi * t[mask] / duration))
        return out
    return _fn


# ---------------------------------------------------------------------------
# Main solver
# ---------------------------------------------------------------------------

class WaveSolver:
    """
    Linearised 1D ray-tracing wave solver for a VascularGraph.

    Parameters
    ----------
    graph            : VascularGraph
    source_node_id   : str
    sensor_node_ids  : list[str]
    source_amplitude : float  — initial pressure [Pa]
    min_amplitude    : float  — prune packets below this [Pa]
    max_reflections  : int    — prune packets exceeding this reflection count
    max_time         : float  — simulation window [s]
    include_stenosis_reflections : bool
                       Generate internal echo packets at stenosis throats.
    """

    def __init__(
        self,
        graph:               VascularGraph,
        source_node_id:      str,
        sensor_node_ids:     List[str],
        source_amplitude:    float = 1.0,
        min_amplitude:       float = 1e-4,
        max_reflections:     int   = 5,
        max_time:            float = 2.0,
        include_stenosis_reflections: bool = True,
    ):
        self.graph               = graph
        self.source_node_id      = source_node_id
        self.sensor_node_ids     = set(sensor_node_ids)
        self.source_amplitude    = source_amplitude
        self.min_amplitude       = min_amplitude
        self.max_reflections     = max_reflections
        self.max_time            = max_time
        self.include_stenosis_reflections = include_stenosis_reflections
        graph.validate()

    # ------------------------------------------------------------------
    # Core simulation
    # ------------------------------------------------------------------

    def run(self) -> Dict[str, SensorRecord]:
        """
        Execute the ray-tracing simulation.

        Returns
        -------
        dict mapping sensor_node_id -> SensorRecord
        """
        sensors = {nid: SensorRecord(nid) for nid in self.sensor_node_ids}
        heap: List[WavePacket] = []

        # Inject initial packets from source node outward
        for edge, neighbour_id in self.graph.neighbours(self.source_node_id):
            dt    = self._edge_transit_time(edge)
            atten = self._edge_amplitude_factor(edge)
            pkt = WavePacket(
                arrival_time    = dt,
                amplitude       = self.source_amplitude * atten,
                next_node_id    = neighbour_id,
                current_edge_id = edge.edge_id,
                path            = [(edge.edge_id, 'fwd')],
                n_reflections   = 0,
            )
            heapq.heappush(heap, pkt)

            # Stenosis on the source-to-neighbour edge spawns an echo
            # back to the source immediately
            if edge.stenosis and self.include_stenosis_reflections:
                self._spawn_stenosis_echo(
                    heap,
                    WavePacket(0.0, self.source_amplitude, self.source_node_id,
                               edge.edge_id, [], 0),
                    edge,
                    self.source_node_id,
                    1.0,   # tau=1 at source (pure injection, no junction)
                )

        # Record source node at t=0
        if self.source_node_id in self.sensor_node_ids:
            sensors[self.source_node_id].record(0.0, self.source_amplitude)

        # ----- Event loop -----
        while heap:
            packet = heapq.heappop(heap)

            if packet.arrival_time > self.max_time:
                continue
            if abs(packet.amplitude) < self.min_amplitude:
                continue
            if packet.n_reflections > self.max_reflections:
                continue

            node_id = packet.next_node_id

            # Record at sensor
            if node_id in self.sensor_node_ids:
                sensors[node_id].record(packet.arrival_time, packet.amplitude)

            # Junction coefficients
            jc = compute_junction_coefficients(
                self.graph, node_id, packet.current_edge_id
            )

            # --- Reflected packet (back along incoming edge) ---
            if abs(jc.gamma) >= self.min_amplitude:
                inc = self.graph.edge(packet.current_edge_id)
                back_node = (
                    inc.node_a if inc.node_b == node_id else inc.node_b
                )
                refl = WavePacket(
                    arrival_time    = packet.arrival_time + self._edge_transit_time(inc),
                    amplitude       = packet.amplitude * jc.gamma * self._edge_amplitude_factor(inc),
                    next_node_id    = back_node,
                    current_edge_id = inc.edge_id,
                    path            = packet.path + [(inc.edge_id, 'bwd')],
                    n_reflections   = packet.n_reflections + 1,
                )
                heapq.heappush(heap, refl)

            # --- Transmitted packets ---
            for eid in jc.outgoing_edge_ids:
                out = self.graph.edge(eid)
                far_node = out.node_b if out.node_a == node_id else out.node_a
                trans_amp = packet.amplitude * jc.tau * self._edge_amplitude_factor(out)
                trans = WavePacket(
                    arrival_time    = packet.arrival_time + self._edge_transit_time(out),
                    amplitude       = trans_amp,
                    next_node_id    = far_node,
                    current_edge_id = eid,
                    path            = packet.path + [(eid, 'fwd')],
                    n_reflections   = packet.n_reflections,
                )
                heapq.heappush(heap, trans)

                # Internal stenosis echo: bounces back from the throat to node_id
                if out.stenosis and self.include_stenosis_reflections:
                    self._spawn_stenosis_echo(
                        heap, packet, out, node_id, jc.tau
                    )

        return sensors

    # ------------------------------------------------------------------
    # Stenosis echo generation
    # ------------------------------------------------------------------

    def _spawn_stenosis_echo(
        self,
        heap:           list,
        incoming_pkt:   WavePacket,
        edge:           VascularEdge,
        return_node_id: str,
        junction_tau:   float,
    ) -> None:
        """
        Spawn an echo packet from the proximal stenosis face back to return_node_id.

        This models the TDR-style reflection: wave enters the edge, travels the
        proximal healthy section, partially reflects at the stenosis throat, and
        returns to the entry node.

        Echo amplitude:
          A_echo = A_incident * tau_junction * atten_prox * gamma_sten * atten_prox

        Echo arrival time:
          t_echo = t_entry + 2 * t_proximal_section
        """
        s   = edge.stenosis
        pos = s.position
        lf  = s.length_fraction
        L   = edge.length_m
        c   = edge.wave_speed

        L_prox = max(pos - lf / 2.0, 0.0) * L
        t_prox = L_prox / c
        atten_prox = math.exp(-edge.viscous_attenuation * L_prox)

        gamma_sten, _ = compute_stenosis_coefficients(edge)

        A_echo = (
            incoming_pkt.amplitude
            * junction_tau
            * atten_prox
            * gamma_sten
            * atten_prox
        )

        if abs(A_echo) < self.min_amplitude:
            return

        echo = WavePacket(
            arrival_time    = incoming_pkt.arrival_time + 2.0 * t_prox,
            amplitude       = A_echo,
            next_node_id    = return_node_id,
            current_edge_id = edge.edge_id,
            path            = incoming_pkt.path + [(edge.edge_id, 'sten_echo')],
            n_reflections   = incoming_pkt.n_reflections + 1,
        )
        heapq.heappush(heap, echo)

    # ------------------------------------------------------------------
    # Edge helpers
    # ------------------------------------------------------------------

    def _edge_transit_time(self, edge: VascularEdge) -> float:
        """
        Full transit time across an edge [s].

        For stenosed edges uses the three-section model:
          t = t_proximal_healthy + t_stenosis_throat + t_distal_healthy
        """
        if edge.stenosis is None or not self.include_stenosis_reflections:
            return edge.transit_time()

        s   = edge.stenosis
        pos = s.position
        lf  = s.length_fraction
        L   = edge.length_m
        c   = edge.wave_speed

        t_prox = max(pos - lf / 2.0, 0.0) * L / c

        ratio = edge.diameter_m / edge.stenosed_diameter()
        c_s   = c * math.sqrt(ratio)
        t_sten = lf * L / c_s

        t_dist = max(1.0 - pos - lf / 2.0, 0.0) * L / c

        return t_prox + t_sten + t_dist

    def _edge_amplitude_factor(self, edge: VascularEdge) -> float:
        """
        Net amplitude factor for a wave traversing the full edge.

        Combines viscous attenuation and transmission through the stenosis throat.
        """
        base = edge.amplitude_factor()
        if edge.stenosis is None or not self.include_stenosis_reflections:
            return base
        _, tau_sten = compute_stenosis_coefficients(edge)
        return base * tau_sten

    # ------------------------------------------------------------------
    # Convenience wrappers
    # ------------------------------------------------------------------

    def waveforms(
        self,
        records:  Dict[str, SensorRecord],
        t:        np.ndarray,
        pulse_fn: Callable,
    ) -> Dict[str, np.ndarray]:
        """Build pressure waveform arrays for all sensors."""
        return {nid: rec.waveform(t, pulse_fn) for nid, rec in records.items()}

    def extract_tof_amplitude(
        self,
        records:             Dict[str, SensorRecord],
        threshold_fraction:  float = 0.1,
    ) -> Dict[str, List[Tuple[float, float]]]:
        """
        Extract (TOF, amplitude) pairs above threshold from impulse records.

        Returns
        -------
        dict mapping sensor_node_id -> sorted list of (t [s], A [Pa])
        """
        cutoff = threshold_fraction * self.source_amplitude
        result = {}
        for nid, rec in records.items():
            filtered = [(t0, A) for t0, A in rec.arrivals if abs(A) >= cutoff]
            filtered.sort(key=lambda x: x[0])
            result[nid] = filtered
        return result
