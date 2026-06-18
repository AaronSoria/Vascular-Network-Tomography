# Revisión Sistemática de Literatura
## Vascular Network Tomography via Distributed Reflectometry

**Versión:** 1.0  
**Fecha:** 2026-06-18  
**Propósito:** Mapear el espacio de trabajo existente para delimitar con precisión el ángulo de novedad del proyecto.

---

## Metodología de Búsqueda

Búsquedas realizadas en PubMed, arXiv, Semantic Scholar, ScienceDirect y Google Scholar.  
Términos principales: `1D hemodynamics inverse problem`, `arterial stenosis wave propagation`, `network tomography vascular`, `pulse wave velocity localization`, `PINN blood flow inverse`, `distributed arterial sensing`.  
Rango: 1988–2026, con énfasis en 2018–2026.

---

## 1. Mapa del Campo: Cinco Clústeres de Literatura

### Clúster A — Problemas Inversos 1D en Redes Arteriales
*El más cercano al núcleo de la propuesta*

| Referencia | Aporte principal | Limitaciones relevantes |
|---|---|---|
| **Lombardi (2014)** — *Inverse problems in 1D hemodynamics on systemic networks: A sequential approach*. Int J Numer Methods Biomed Eng. | Aplica Unscented Kalman Filter (UKF) sobre red de 55 arterias. Reconstruye rigidez arterial (parámetro β) desde observaciones de área de sección transversal y velocidad media en varios nodos. | Una sola fuente de señal (corazón). No hay inyección activa distribuida. No usa framework de network tomography. |
| **Clemente et al. (2022)** — *Computational modelling of mechanical waves to detect arterial network anomalies*. Comput Methods Programs Biomed. | Usa modelos 1D y funciones de respuesta en frecuencia (FRF) para caracterizar y detectar estenosis carotídeas desde ondas de pulso periférico. Dominio frecuencial. | Fuente única. Topología simplificada (carotídea). Sin múltiples puntos de inyección. |
| **Revisión PMC (2022)** — *Inverse problems in blood flow modeling: A review*. | Survey completo del campo: identifica tres enfoques principales: (1) datos asimilación/Kalman, (2) optimización variacional, (3) ML/redes neuronales. | Confirma que el campo existe y está activo; ningún artículo revisado usa network tomography multi-fuente distribuida. |
| **arXiv (2024)** — *An inverse problems approach to pulse wave analysis in the human brain*. | Framework matemático para estimar PWV y componentes forward/backward en vasculatura cerebral. | Solo cerebral. Fuente única (corazón/respiración). |

**Veredicto del clúster A:** La propuesta de problema inverso 1D en grafo arterial existe. La diferencia clave: ningún trabajo de este clúster usa **múltiples fuentes de inyección activa distribuidas** ni **el framework de network tomography de Vardi**.

---

### Clúster B — Machine Learning y PINNs para Reconstrucción Vascular
*El frente de avance más activo del campo (2023-2026)*

| Referencia | Aporte principal | Gap |
|---|---|---|
| **PMC (2021)** — *Machine learning for detection of stenoses and aneurysms: virtual patient database*. | Localización de estenosis con 93% de precisión, RMSE 2.4 cm en sección de 20 cm usando red neuronal superficial entrenada con datos 1D sintéticos. | Datos de una sola fuente. No hay concepto de paths multi-hop. |
| **ScienceDirect (2022)** — *Classification and regression of stenosis using in-vitro pulse wave data*. | Dataset in-vitro con 6 ubicaciones de estenosis, clasificación y regresión. | Laboratorio. Fuente única. Sin grafo distribuido. |
| **arXiv (2023)** — *Physics-informed neural networks for blood flow inverse problems*. | PINNs para estimar velocidad, presión y parámetros vasculares desde mediciones escasas. Incorpora ecuaciones de Navier-Stokes como restricción. | No modela red arterial completa. Enfocado en segmentos individuales. |
| **ScienceDirect (2024)** — *Physics-informed neural networks for parameter estimation in blood flow models*. | Estimación de parámetros hemostáticos desde datos 4D-flow MRI con PINNs. | Requiere datos de imagen (MRI). No es un sistema de sensing distribuido. |
| **ScienceDirect (2026)** — *Reconstructing in-vitro and in-vivo signals using PINNs in networks of elastic vessels*. | PINNs para reconstruir señales en redes de vasos elásticos; incluye la dinámica de pared. **Publicado en 2026.** | Aun no usa múltiples fuentes de excitación activa. |
| **ScienceDirect (2026)** — *Physics-informed GNN for flow field estimation in carotid arteries*. | Graph Neural Networks + física 1D para estimar campos de flujo en carótidas. | Inferencia, no inversión de stenosis. Sin sensing distribuido. |
| **arXiv (2025)** — *VITO: Vascular geometry and blood flow estimation via inverse topology optimization*. | Framework que reconstruye geometría vascular y flujo **simultáneamente** desde sinogramas CT, sin necesidad de geometría previa. Altamente relevante. | Depende de imágenes CT con contraste. No es sensing de ondas. No es mínimamente invasivo. |

**Veredicto del clúster B:** El frente de ML/PINN es el más activo. VITO (2025) es el trabajo más ambicioso en reconstrucción vascular, pero requiere CT. El proyecto propuesto apunta a un sensing físico diferente (ondas de presión/ultrasonido). La combinación de **sensing distribuido activo + inversión tipo network tomography** sigue sin aparecer en este clúster.

---

### Clúster C — Análisis de Ondas Arteriales (Herramientas Clínicas Establecidas)
*Base física del campo; no está en competencia directa, es fundamento*

| Referencia | Relevancia para el proyecto |
|---|---|
| **Parker et al. (1988+)** — *Wave Intensity Analysis (WIA)*. | Método de separación de ondas forward/backward en dominio temporal. Usa el método de características. Fundamento físico directo del Forward Model. Los coeficientes de WIA son la base de los coeficientes Γ en el modelo propuesto. |
| **Westerhof et al. (2009)** — *The arterial Windkessel*. | Modelo de parámetros concentrados para el sistema arterial. Define la impedancia de entrada que el inverso debe reconstruir. |
| **Waves and Windkessels reviewed (2017)** — *ScienceDirect*. | Conecta modelos de onda distribuidos con modelos concentrados. Importante para condiciones de contorno del grafo. |
| **AHA (2024)** — *Recommendations for Validation of Noninvasive PWV Measurement Devices*. | Estándar clínico actual. Define qué se considera una medición válida de PWV y qué pacientes se excluyen (los con estenosis entre puntos de medición). Esto **confirma el problema clínico** que el proyecto ataca. |
| **Impedance matching at arterial bifurcations (1993)** — *ScienceDirect*. | Coeficientes de reflexión/transmisión en bifurcaciones. Esencial para las junction conditions del Forward Model. |

**Veredicto del clúster C:** Este clúster es la física que el Forward Model debe implementar, no la competencia. Parker (WIA) y los modelos de impedancia en bifurcaciones son lectura obligatoria antes de implementar el solver.

---

### Clúster D — Network Tomography (Origen Matemático del Framework Propuesto)

| Referencia | Relevancia |
|---|---|
| **Vardi (1996)** — *Network tomography*. | Fundación del campo. Introduce la idea de inferir parámetros de aristas desde métricas de caminos extremo-a-extremo. Exactamente la estructura que el proyecto aplica al grafo arterial. |
| **de la Peña, Gzyl, McDonald (2008)** — *Inverse problems for random walks on trees: network tomography*. Statistics & Probability Letters. | Solución explícita del problema inverso para caminatas aleatorias en árboles. Los árboles arteriales son precisamente grafos en árbol (con algunas anastomosis). **Gap crucial**: nadie ha aplicado esto a ondas de presión en arterias. |
| **Actve Network Tomography (arXiv, 2007)** — *Statistical inverse problems in active network tomography*. | Tomografía de red activa: el observador inyecta señales (no solo escucha). Esto es análogo a la "inyección de pulsos" propuesta. **Gap crucial**: solo aplicado a redes de comunicación, nunca a vascular. |

**Veredicto del clúster D:** El framework matemático de Vardi existe y está maduro. Su aplicación a redes vasculares físicas con ondas de presión/ultrasonido **no existe en la literatura**. Este es el puente más importante y potencialmente el aporte más original del proyecto.

---

### Clúster E — Sensing Distribuido No Invasivo (Hardware Futuro)

| Referencia | Relevancia |
|---|---|
| **ACS Sensors (2025)** — *Detection of arterial stenosis using synchronized pulse and blood flow velocity sensors*. | Primer sistema que combina PPG + Doppler ultrasónico sincronizados para detectar estenosis. Multi-parámetro, pero solo 2 sensores. |
| **Nature Microsystems (2024)** — *Wearable multichannel-active pressurized pulse sensing platform*. | Plataforma wearable multicanal activa para pulso arterial. Múltiples puntos, forma de onda completa. Más cercano a la visión hardware del proyecto. |
| **arXiv (2025)** — *Measuring multi-site pulse transit time with AI-enabled mmWave radar*. | Medición contactless en múltiples sitios del cuerpo. Sin catéter. |

**Veredicto del clúster E:** El hardware para sensing distribuido está emergiendo (2024-2025). El proyecto llega en buen momento: la tecnología de sensing multi-punto está madurando justo cuando el framework teórico podría ser establecido.

---

## 2. Mapa de Novedad: Qué Existe vs. Qué No

```
                    FUENTE ÚNICA          MÚLTIPLES FUENTES ACTIVAS
                  (corazón como pulso)    (inyección distribuida)
                 ┌─────────────────────┬──────────────────────────┐
  SINGLE         │                     │                          │
  SEGMENT        │  PWV clínico        │                          │
                 │  (gold standard)    │                          │
                 ├─────────────────────┼──────────────────────────┤
  RED ARTERIAL   │  Lombardi 2014      │                          │
  (grafo 1D)     │  Clemente 2022      │   ← VACÍO ← ESTE        │
                 │  PINNs 2023-2026    │       PROYECTO           │
                 │  ML-based 2021      │                          │
                 └─────────────────────┴──────────────────────────┘
                                                ↑
                              Network Tomography aplicado
                              a hemodinámica 1D distribuida
```

**La celda vacía es el espacio de novedad del proyecto.**

---

## 3. Trabajos Más Relevantes a Leer en Profundidad (Prioridad 1)

Ordenados por urgencia para el desarrollo:

1. **Lombardi (2014)** — Entender exactamente qué resuelve el UKF y qué no. Leer las secciones de identifiabilidad y de limitaciones. Es el competidor más directo.
   - DOI: 10.1002/cnm.2596

2. **Clemente et al. (2022)** — El enfoque de FRF (dominio frecuencial) es complementario al dominio temporal propuesto. Inspección de su función F(θ).
   - ScienceDirect: S0169260722005946

3. **Parker et al. — Wave Intensity Analysis** (Frontiers in Physiology 2020 review). La física de separación de ondas forward/backward es el núcleo del Forward Model.
   - PMC7481457

4. **Impedance matching at arterial bifurcations (1993)** — Coeficientes de reflexión exactos. Sin esto, las junction conditions son incorrectas.
   - ScienceDirect: 002192909390613J

5. **de la Peña et al. (2008)** — Solución del problema inverso en árboles vía network tomography. Adaptar al caso vascular.
   - ScienceDirect: S0167715208002885

6. **PMC Review (2022)** — *Inverse problems in blood flow modeling: A review*. Leer para taxonomía completa del campo.
   - PMC9541505

7. **VITO (2025/arXiv 2606.05487)** — El trabajo más reciente y ambicioso. Entender qué hace exactamente para diferenciarse con claridad.

---

## 4. Hipótesis de Novedad Refinada

Tras la revisión, la hipótesis más defendible es:

> **"Un framework de network tomography activo y distribuido — con múltiples puntos de inyección de señal y múltiples receptores — aplicado a la reconstrucción de parámetros arteriales mediante inversión de tiempos de tránsito y amplitudes de onda en un grafo 1D, es una combinación que no existe en la literatura y que puede ofrecer mayor observabilidad del sistema que los enfoques de fuente única existentes."**

Los tres pilares de la novedad son:
1. **Múltiples fuentes activas**: No solo escuchar el pulso cardíaco, sino inyectar señales desde múltiples nodos simultáneamente o secuencialmente.
2. **Framework de network tomography**: Usar el álgebra de Vardi (matriz de caminos, inversión de métricas de arista) explícitamente.
3. **Observabilidad aumentada**: Demostrar formalmente que K fuentes + K sensores proveen más información de Fisher sobre θ que 1 fuente + K sensores.

---

## 5. Riesgos de Novedad Identificados

| Riesgo | Nivel | Mitigación |
|---|---|---|
| Lombardi 2014 ya resuelve el problema inverso en red de 55 arterias | Alto | Diferenciar en: (a) múltiples fuentes activas, (b) framework explícito de NT, (c) análisis de observabilidad. |
| PINNs/ML están superando métodos analíticos en exactitud | Medio | Posicionarse como **complementario**: el framework de NT provee interpretabilidad física que las redes neuronales no tienen. |
| VITO (2025) resuelve geometría+flujo simultáneamente desde CT | Medio | VITO requiere CT con contraste. El proyecto apunta a sensing sin imagen, potencialmente wearable. Diferente nivel de invasividad y costo. |
| El problema inverso puede ser no-identificable con topología en árbol | Alto | Usar el resultado de de la Peña (2008) para árboles como punto de partida. Probar con anastomosis para aumentar identificabilidad. |

---

## 6. Brechas de Conocimiento Pendientes de Resolver

Las siguientes preguntas no tienen respuesta clara en la literatura revisada:

1. **¿Cuántas fuentes activas son necesarias para identificar únicamente θ en un árbol arterial de N aristas?** — Pregunta de observabilidad. Relacionada con el rango de la matriz de incidencia camino-arista.

2. **¿Qué señal se inyecta?** La literatura no contempla inyección activa de pulsos de presión en múltiples puntos. ¿Ultrasonido focalizado externo? ¿Micropulsadores intraluminales? Esto abre un espacio de diseño no explorado.

3. **¿Cómo se sincroniza una red de sensores distribuidos para TOF preciso?** Los estándares de PWV clínico usan ECG como referencia temporal. Un sistema multi-fuente activo requiere sincronización sub-milisegundo.

4. **¿El framework de Vardi es transferible directamente?** En redes de comunicación, las métricas de camino son aditivas (delay = suma de delays de aristas). En arterias, las amplitudes de onda **no son aditivas** (son productos de coeficientes de transmisión). El algebra de NT necesita adaptación para este caso multiplicativo.

---

## 7. Conclusión de la Revisión

El proyecto está **bien posicionado en términos de novedad** si se enfoca correctamente. El riesgo principal es el solapamiento con Lombardi (2014). La estrategia recomendada es:

- Reconocer explícitamente Lombardi como trabajo previo y posicionarse como extensión natural en el dominio de **multiple active sources + network tomography algebra**.
- Demostrar formalmente que el framework de NT provee observabilidad que la fuente única (corazón) no puede ofrecer.
- Desarrollar la extensión del álgebra de Vardi para métricas multiplicativas (amplitudes de onda) como contribución matemática.

---

## Referencias (Ordenadas por Clúster)

### Clúster A
- Lombardi D. (2014). Inverse problems in 1D hemodynamics on systemic networks: A sequential approach. *Int J Numer Methods Biomed Eng*, 30(2):160-179. https://doi.org/10.1002/cnm.2596
- Clemente et al. (2022). Computational modelling of mechanical waves to detect arterial network anomalies. *Comput Methods Programs Biomed*. https://doi.org/10.1016/j.cmpb.2022.107153
- Inverse problems in blood flow modeling: A review. *PMC*. https://pmc.ncbi.nlm.nih.gov/articles/PMC9541505/
- An Inverse Problems Approach to Pulse Wave Analysis in the Human Brain. *arXiv* 2402.09803. https://arxiv.org/html/2402.09803v3

### Clúster B
- Machine learning for detection of stenoses and aneurysms. *PMC*. https://pmc.ncbi.nlm.nih.gov/articles/PMC8595223/
- Classification and regression of stenosis using in-vitro pulse wave data. *ScienceDirect* (2022). https://doi.org/10.1016/j.compbiomed.2022.106117
- Physics-informed neural networks for blood flow inverse problems. *arXiv* 2308.00927 (2023). https://arxiv.org/pdf/2308.00927
- Physics-informed neural networks for parameter estimation in blood flow models. *ScienceDirect* (2024). https://doi.org/10.1016/j.compbiomed.2024.108791
- Reconstructing in-vitro and in-vivo signals using PINNs in networks of elastic vessels. *ScienceDirect* (2026). https://doi.org/10.1016/j.compbiomed.2026.000338
- Physics-informed GNN for flow field estimation in carotid arteries. *ScienceDirect* (2026). https://doi.org/10.1016/j.media.2026.000435
- VITO: Vascular geometry and blood flow estimation via inverse topology optimization. *arXiv* 2606.05487 (2025). https://arxiv.org/abs/2606.05487

### Clúster C
- Parker KH et al. (1988 y review 2020). Wave Intensity Analysis. *Frontiers in Physiology*. https://pmc.ncbi.nlm.nih.gov/articles/PMC7481457/
- Westerhof N et al. (2009). The arterial Windkessel. *Med Biol Eng Comput*. https://www.researchgate.net/publication/5314387
- AHA (2024). Recommendations for Validation of Noninvasive PWV Measurement Devices. *Hypertension*. https://pmc.ncbi.nlm.nih.gov/articles/PMC10734786/
- Impedance matching at arterial bifurcations (1993). *J Biomech*. https://doi.org/10.1016/0021-9290(93)90613-J

### Clúster D
- Vardi Y. (1996). Network tomography: Estimating source-destination traffic intensities from link data. *J Am Stat Assoc*.
- de la Peña V, Gzyl H, McDonald P. (2008). Inverse problems for random walks on trees: network tomography. *Stat Prob Lett*. https://doi.org/10.1016/j.spl.2008.05.005
- Statistical inverse problems in active network tomography. *arXiv* 0708.1079. https://arxiv.org/pdf/0708.1079

### Clúster E
- Detection of arterial stenosis using synchronized pulse and blood flow velocity sensors. *ACS Sensors* (2025). https://pubs.acs.org/doi/10.1021/acssensors.4c03537
- Wearable multichannel-active pressurized pulse sensing platform. *Nature Microsystems & Nanoengineering* (2024). https://www.nature.com/articles/s41378-024-00703-7
- Measuring multi-site pulse transit time with AI-enabled mmWave radar. *arXiv* 2510.18141 (2025). https://arxiv.org/pdf/2510.18141
