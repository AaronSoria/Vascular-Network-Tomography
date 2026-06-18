# Experimento de Validación Mínimo
## Vascular Network Tomography — Phase I Validation

**Versión:** 1.0  
**Fecha:** 2026-06-18  
**Propósito:** Definir el experimento más pequeño posible que sea suficiente para demostrar (o falsificar) la hipótesis central.

---

## 1. Hipótesis a Demostrar

> "Las observaciones M = {T_ij, A_ij} de un sistema distribuido de K≥2 sensores son suficientes para localizar correctamente una estenosis en un grafo arterial conocido, bajo condiciones ideales (sin ruido, topología conocida)."

Este es el mínimo científicamente defensible: si la hipótesis falla bajo condiciones ideales, no tiene sentido avanzar a condiciones realistas.

---

## 2. Topología Elegida: Y Bifurcation (Nivel 2)

```
[source=A] ——— e_parent (35cm, D=25mm, c=4.5m/s) ——— [junction=J]
                                                         |
                              ┌──────────────────────────┤
                              |                          |
                       e_left (7cm, D=10mm, c=7m/s)  e_right (7cm, D=10mm, c=7m/s)
                              |                          |
                         [sensor=L]                [sensor=R]

Sensores: A (source), L, R  →  K=3
```

**Justificación de la elección:**
- Topología mínima no trivial (tiene una bifurcación real)
- 3 edges → 9 parámetros a reconstruir (L, D, c por arista)
- 3 sensores × múltiples llegadas → O(10-30) observaciones independientes
- Existe solución analítica para validar el inverso
- Simetría izquierda/derecha: un buen inverso debe romper la simetría cuando la estenosis está solo en un lado

---

## 3. Protocolo del Experimento

### 3.1 Generación de Ground Truth

```
Para cada escenario s ∈ {Healthy, Stenosis_L_30%, Stenosis_L_50%, Stenosis_L_70%, Stenosis_R_70%, Bilateral_50%}:
  1. Construir el grafo Y con parámetros fisiológicos conocidos (ground_truth)
  2. Ejecutar WaveSolver (Forward Model) → obtener M_obs = {T_ij, A_ij}
  3. Guardar (ground_truth, M_obs, escenario) como dataset
```

**Parámetros fisiológicos (ground truth):**

| Edge | L [cm] | D [mm] | c [m/s] |
|---|---|---|---|
| e_parent | 35.0 | 25.0 | 4.5 |
| e_left | 7.0 | 10.0 | 7.0 |
| e_right | 7.0 | 10.0 | 7.0 |

**Estenosis target:** `D_stenosed = D * (1 - severity)` solo en e_left o e_right.

### 3.2 Formulación del Problema Inverso

**Input al inverso:** Solo `M_obs` (no se provee ground_truth).

**Variables a estimar:**
```
θ = {D_left, D_right}   # solo diámetros — L y c se asumen conocidos en esta validación
```

**Reducción justificada:** En la validación mínima, asumimos que L y c son conocidos (obtenidos de un atlas anatómico o medición previa). Solo D es incógnita porque la estenosis lo modifica. Esto reduce el problema de 9 a 2 variables y permite análisis de identifiabilidad exacto.

**Función objetivo:**
```
θ̂ = argmin_{θ} ||F(θ) - M_obs||²  +  λ ||θ - θ_prior||²

donde:
  F(θ) = WaveSolver(graph(θ))        # Forward model
  M_obs = vector de (T_ij, A_ij)     # observaciones
  θ_prior = D_healthy = 10 mm        # regularización hacia diámetro sano
  λ = parámetro de regularización (buscar por validación cruzada)
```

### 3.3 Observaciones Utilizadas

Para el experimento mínimo, usar solo las **primeras llegadas** en cada sensor:

| Observación | Descripción | Sensibilidad esperada |
|---|---|---|
| T_ij (A→L) | TOF primera llegada A→L | Baja (0.4ms por 70% stenosis) |
| A_ij (A→L) | Amplitud primera llegada A→L | Alta (90% cambio por 70% stenosis) |
| T_ij (A→R) | TOF primera llegada A→R | Cero (e_right intacto) |
| A_ij (A→R) | Amplitud primera llegada A→R | Baja (~15% cambio por e_left stenosis) |
| T_TDR (A) | TOF del eco TDR en A | Alta (varía con posición de estenosis) |
| A_TDR (A) | Amplitud del eco TDR en A | Alta (varía con severidad) |

**Vector de observaciones:** `M_obs = [T_AL, A_AL, T_AR, A_AR, T_TDR, A_TDR]` → 6 escalares.

Con 2 incógnitas y 6 observaciones el sistema está **sobredeterminado** (redundante) — buena señal para el inverso.

### 3.4 Criterios de Éxito

El experimento se declara exitoso si:

| Criterio | Umbral |
|---|---|
| Localización correcta | La arista estimada como stenótica coincide con la real (accuracy = 100% en validación sin ruido) |
| Error de diámetro | |D_estimated - D_true| < 0.5 mm para ≥ 80% de escenarios |
| Error de severidad | |s_estimated - s_true| < 0.05 para ≥ 80% de escenarios |
| Convergencia | Optimizador converge en < 100 iteraciones |
| No falsos positivos | Con escenario healthy, D_estimated ≈ D_true (error < 5%) |

### 3.5 Criterios de Falsificación

El experimento refuta la hipótesis si:

- El mínimo del objetivo es plano (gradiente ≈ 0) respecto a D → problema no identificable
- Múltiples θ distintos producen el mismo M_obs (no unicidad sin regularización)
- El error de localización > 50% incluso sin ruido → los observables no contienen suficiente información

---

## 4. Análisis de Identifiabilidad (antes del inverso)

Antes de implementar el optimizador, verificar que las observaciones son sensibles a los parámetros.

**Test de sensibilidad:** Calcular la Jacobiana numérica ∂M/∂θ mediante diferencias finitas:

```python
def jacobian(theta_0, forward_fn, eps=1e-4):
    M_0 = forward_fn(theta_0)
    J = []
    for i, th in enumerate(theta_0):
        theta_plus = theta_0.copy()
        theta_plus[i] += eps
        M_plus = forward_fn(theta_plus)
        J.append((M_plus - M_0) / eps)
    return np.array(J).T   # shape: (n_obs, n_params)
```

**Criterio de identifiabilidad:** Si `rank(J) == len(theta_0)` → el sistema es (localmente) identificable.

Esto debe verificarse para cada escenario y para el vector completo de observaciones.

---

## 5. Plan de Implementación

### Fase A: Generación del Dataset (semana 1)

```python
# experiments/generate_dataset.py
scenarios = [
    make_y_bifurcation(),                              # healthy
    make_y_bifurcation('e_left', 0.30),                # mild
    make_y_bifurcation('e_left', 0.50),                # moderate
    make_y_bifurcation('e_left', 0.70),                # severe
    make_y_bifurcation('e_right', 0.70),               # right side
    make_y_bifurcation('e_left', 0.50, e_right=0.50),  # bilateral (Phase II)
]
# → guarda JSON: {scenario_name, ground_truth, M_obs}
```

### Fase B: Análisis de Identifiabilidad (semana 1)

```python
# experiments/identifiability_analysis.py
for scenario in scenarios:
    J = jacobian(theta_healthy, forward_fn)
    rank, cond_number = np.linalg.matrix_rank(J), np.linalg.cond(J)
    print(f'{scenario}: rank={rank}, cond={cond_number:.1e}')
```

### Fase C: Inverso Simple (semana 2)

```python
# src/inverse/least_squares.py
from scipy.optimize import minimize

def inverse_solve(M_obs, forward_fn, theta_init, lambda_reg=1e-3):
    theta_prior = theta_init.copy()
    def objective(theta):
        M_pred = forward_fn(theta)
        residual = np.sum((M_pred - M_obs)**2)
        regularizer = lambda_reg * np.sum((theta - theta_prior)**2)
        return residual + regularizer
    result = minimize(objective, theta_init, method='L-BFGS-B',
                      bounds=[(1e-3, 25e-3)] * len(theta_init))
    return result.x, result.fun
```

### Fase D: Evaluación (semana 2)

```python
# experiments/run_inverse.py
for scenario in test_scenarios:
    theta_estimated, loss = inverse_solve(scenario.M_obs, forward_fn, theta_init=theta_healthy)
    report = evaluate(scenario.ground_truth, scenario.stenosis_edges,
                      ReconstructionResult(theta_estimated))
    print(report.summary())
```

---

## 6. Extensiones Inmediatas (si el mínimo tiene éxito)

Una vez que el experimento mínimo (sin ruido, topología conocida) es exitoso:

1. **Ruido gaussiano:** Agregar σ = 1%, 5%, 10% a M_obs y repetir.
2. **Incertidumbre en topología:** Perturbar L y c en ±5% y observar degradación.
3. **Estenosis bilateral:** Dos incógnitas independientes — probar identifiabilidad.
4. **Árbol de 5 generaciones:** Escalar el inverso a 13 aristas y 7 sensores.
5. **Múltiples fuentes activas:** Inyectar desde 3 nodos distintos — verificar que aumenta el rango de la Jacobiana.

---

## 7. Métricas de Reporte (para publicación)

| Métrica | Definición | Unidad |
|---|---|---|
| E_D | RMSE de diámetro estimado vs real | mm |
| E_s | MAE de severidad estimada vs real | — |
| E_loc | % de escenarios con arista correcta identificada | % |
| Sensitivity | TP / (TP + FN) | — |
| Specificity | TN / (TN + FP) | — |
| F1 | 2·P·R / (P + R) | — |
| N_iter | Iteraciones del optimizador hasta convergencia | — |
| Δt_TDR | Diferencia de TOF entre healthy y stenosed | ms |

---

## 8. Ejemplo Numérico Esperado (predicción a priori)

Basado en los resultados del Forward Model (Phase I):

**Escenario: e_left 70% stenosis**

```
Observaciones generadas por WaveSolver:
  T_AL  = 87.33 ms   (vs 87.78 ms healthy → Δ = -0.45 ms)
  A_AL  = 3.10 Pa    (vs 1.63 Pa healthy → +90%)
  T_AR  = 87.78 ms   (sin cambio)
  A_AR  = 1.63 Pa    (sin cambio significativo)
  T_TDR = 164.56 ms  (solo aparece con stenosis → 9ms después de eco J)
  A_TDR = 0.249 Pa   (0 en healthy)

Estimación esperada del inverso:
  D_left_estimated  ≈ 3.0 mm (true: 3.0 mm = 10 mm × 0.30)
  D_right_estimated ≈ 10.0 mm (true: 10 mm)
  → Localización correcta (left): ✓
  → Error de diámetro: < 0.5 mm → ✓ (criterio de éxito)
```

---

## 9. Documentación de Resultados

Todos los resultados se guardarán en:
```
experiments/
├── dataset_y_bifurcation.json    # ground truth + M_obs por escenario
├── identifiability_report.csv    # rank y cond de Jacobiana por escenario
├── inverse_results.json          # θ_estimated vs θ_true
└── figures/
    ├── jacobian_heatmap.png
    ├── inverse_convergence.png
    └── reconstruction_comparison.png
```
