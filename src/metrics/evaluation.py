"""
Evaluation Metrics
==================
Quantitative metrics for assessing inverse problem reconstruction quality.

Notation follows the research proposal:
  E_L   = localization error  [m]
  E_S   = severity error      [dimensionless, 0-1]
  Coverage = N_detected / N_total
  Sensitivity = TP / (TP + FN)
  Specificity = TN / (TN + FP)
  F1    = 2*TP / (2*TP + FP + FN)
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple


@dataclass
class StenosisDetection:
    """One detection event from the reconstruction algorithm."""
    edge_id:    str
    severity:   float   # estimated stenosis severity [0, 1]
    confidence: float   # algorithm confidence score [0, 1]


@dataclass
class ReconstructionResult:
    """Output of the inverse solver for one scenario."""
    estimated_params: Dict[str, Dict]   # edge_id → {diameter_m, wave_speed, ...}
    detections:       List[StenosisDetection] = field(default_factory=list)


@dataclass
class EvaluationReport:
    """Aggregated metrics for one reconstruction."""
    # --- Continuous errors ---
    diameter_rmse_m:    float = 0.0   # root mean squared error on diameter [m]
    wavespeed_rmse:     float = 0.0   # RMSE on wave speed [m/s]
    severity_mae:       float = 0.0   # mean absolute error on stenosis severity
    localization_error_m: float = 0.0 # mean edge-localisation error [m]

    # --- Classification metrics ---
    tp: int = 0
    fp: int = 0
    fn: int = 0
    tn: int = 0

    @property
    def sensitivity(self) -> float:
        return self.tp / (self.tp + self.fn) if (self.tp + self.fn) > 0 else 0.0

    @property
    def specificity(self) -> float:
        return self.tn / (self.tn + self.fp) if (self.tn + self.fp) > 0 else 0.0

    @property
    def precision(self) -> float:
        return self.tp / (self.tp + self.fp) if (self.tp + self.fp) > 0 else 0.0

    @property
    def f1(self) -> float:
        p, r = self.precision, self.sensitivity
        return 2 * p * r / (p + r) if (p + r) > 0 else 0.0

    @property
    def coverage(self) -> float:
        total = self.tp + self.fn
        return self.tp / total if total > 0 else 0.0

    def summary(self) -> str:
        return (
            f"EvaluationReport:\n"
            f"  Diameter RMSE : {self.diameter_rmse_m*1e3:.3f} mm\n"
            f"  Wave speed RMSE: {self.wavespeed_rmse:.3f} m/s\n"
            f"  Severity MAE   : {self.severity_mae:.4f}\n"
            f"  Localization Err: {self.localization_error_m*1e2:.2f} cm\n"
            f"  Sensitivity    : {self.sensitivity:.3f}\n"
            f"  Specificity    : {self.specificity:.3f}\n"
            f"  Precision      : {self.precision:.3f}\n"
            f"  F1 score       : {self.f1:.3f}\n"
            f"  Coverage       : {self.coverage:.3f}\n"
            f"  TP={self.tp} FP={self.fp} FN={self.fn} TN={self.tn}"
        )


def evaluate(
    ground_truth:  Dict[str, Dict],         # edge_id → true params
    stenosis_edges: List[str],              # truly stenosed edges
    result:        ReconstructionResult,
    detection_threshold: float = 0.10,     # min severity to count as detection
) -> EvaluationReport:
    """
    Compare reconstruction against ground truth.

    Parameters
    ----------
    ground_truth      : true edge parameters (from Scenario.ground_truth)
    stenosis_edges    : list of edge_ids that are truly stenosed
    result            : output of inverse solver
    detection_threshold : minimum estimated severity to call a detection

    Returns
    -------
    EvaluationReport
    """
    report = EvaluationReport()

    # --- Continuous parameter errors ---
    d_sq_errors, c_sq_errors = [], []
    for eid, true_p in ground_truth.items():
        est_p = result.estimated_params.get(eid)
        if est_p is None:
            continue
        d_sq_errors.append((true_p['diameter_m'] - est_p.get('diameter_m', true_p['diameter_m'])) ** 2)
        c_sq_errors.append((true_p['wave_speed'] - est_p.get('wave_speed', true_p['wave_speed'])) ** 2)

    if d_sq_errors:
        report.diameter_rmse_m = math.sqrt(sum(d_sq_errors) / len(d_sq_errors))
    if c_sq_errors:
        report.wavespeed_rmse = math.sqrt(sum(c_sq_errors) / len(c_sq_errors))

    # --- Classification: TP/FP/FN/TN ---
    detected_edges = {
        det.edge_id for det in result.detections
        if det.severity >= detection_threshold
    }
    stenosis_set = set(stenosis_edges)
    all_edges    = set(ground_truth.keys())

    report.tp = len(detected_edges & stenosis_set)
    report.fp = len(detected_edges - stenosis_set)
    report.fn = len(stenosis_set - detected_edges)
    report.tn = len((all_edges - detected_edges) - stenosis_set)

    # --- Severity error (for TP detections only) ---
    sev_errors = []
    for det in result.detections:
        if det.edge_id in stenosis_set:
            # Ground truth severity: look for stenosis object
            # (passed here as a float in ground_truth for convenience)
            true_sev = ground_truth[det.edge_id].get('stenosis_severity', 0.0)
            sev_errors.append(abs(true_sev - det.severity))
    if sev_errors:
        report.severity_mae = sum(sev_errors) / len(sev_errors)

    return report
