# Análisis Comparativo de Versiones de Hardware
## Vascular Network Tomography — Version A vs Version B

---

## Tabla Comparativa Principal

| Criterio | Version A — Manguito (P-TDR) | Version B — Ultrasonido Guiado (A-GWT) |
|---|---|---|
| **Principio físico** | Onda de Moens-Korteweg (c ≈ 5 m/s) en lumen | Ondas guiadas en pared arterial (c ≈ 1500 m/s) |
| **Frecuencia de trabajo** | 10–500 Hz (subaudio) | 20–100 kHz (ultrasonido bajo) |
| **Resolución espacial** | ~5 mm | ~5 mm (comparables con mismo BW) |
| **Alcance máximo** | 50–100 cm (limitado por atenuación viscosa) | 30–50 cm (limitado por acoplamiento) |
| **Atenuación en ruta** | ~2–5 dB/m (onda de presión arterial) | ~0.025 dB/cm a 50 kHz (casi nula) |
| **SNR teórico del eco** | 47 dB (con 100 promedios) | 60+ dB (en práctica 40–60 dB reales) |
| **Resolución temporal** | 1 ms | 6.25 μs |
| **Costo del prototipo** | ~845 USD | ~1160 USD |
| **Invasividad** | No invasivo | No invasivo |
| **Incomodidad para el paciente** | Baja (familiar, como toma de tensión) | Media (gel + transductores fijados con vendaje) |
| **Complejidad de hardware** | Baja | Media-Alta |
| **Complejidad de software de señal** | Media (correlación en tiempo) | Alta (separación modal, dispersión) |
| **Tecnología base probada** | Sí (Arteriograph existe desde 2004) | Parcialmente (GWT existe en NDT industrial) |
| **Adaptación al modelo matemático** | Directa (mismo fenómeno físico que el solver) | Indirecta (modo distinto, necesita adaptación) |
| **Zonas anatómicas accesibles** | Extremidades (braquial, muslo, pantorrilla) | Cualquier zona con acceso cutáneo |
| **Limitación principal** | Ancho de banda neumático (<500 Hz) | Dispersión modal + acoplamiento |
| **Ruido fisiológico dominante** | Pulso cardíaco (superpuesto, mismo rango espectral) | Movimiento arterial (Doppler espurio) |
| **Tiempo de medición** | 10–15 min | 5 min |
| **Riesgo de artefactos** | Alto (interferencia cardíaca, movimiento) | Medio (reflexiones de tejido) |

---

## Análisis por Dimensión

### 1. Fidelidad física al modelo matemático

**Version A gana claramente.** El solver implementado en Phase I modela exactamente la onda de Moens-Korteweg. Los parámetros que el solver usa (impedancia Z = ρc/A, coeficiente de reflexión Γ, tiempo de tránsito T = L/c) son los mismos que se medirían con la Version A. No hay brecha entre el modelo y el hardware.

Version B excita un fenómeno distinto (onda guiada en la pared, no onda de presión en el lumen). Habría que extender el modelo para incluir la dispersión modal y el acoplamiento pared-lumen. Esto añade una capa de complejidad al problema inverso.

### 2. Resolución espacial en la práctica

Ambas versiones logran ~5 mm de resolución con la misma estrategia (chirp de 80 kHz BW equivalente en Version B, 500 Hz en Version A). Sin embargo, hay una diferencia importante:

- Version A: resolución limitada por el **tiempo de propagación** (c_lento × Δt)
- Version B: resolución limitada por el **ancho de banda** del chirp y la dispersión modal

En Version B, la dispersión modal puede degradar la resolución en la práctica si no se compensa correctamente (la frecuencia más baja viaja más lento que la más alta → el eco aparece "estirado" en el tiempo).

### 3. Accesibilidad anatómica

**Version B gana.** El manguito solo puede aplicarse en extremidades donde hay espacio físico para rodearlo. Muchas arterias clínicamente importantes (carótida, mesentérica, renal) no son accesibles con manguito pero sí con un transductor sobre la piel.

### 4. Riesgo de interferencia fisiológica

**Version B gana marginalmente.** La frecuencia de trabajo de la Version A (10–500 Hz) solapa con los armónicos del pulso cardíaco (1 Hz, 2 Hz, ..., hasta ~20–30 Hz). La excitación debe hacerse por encima de 20 Hz para evitar el solapamiento. En Version B, la banda de trabajo (20–100 kHz) está completamente libre de ruido fisiológico.

### 5. Ruta crítica hacia un prototipo funcional

**Version A es más directa:**
1. Componentes comerciales estándar (manguitos médicos, sensores MEMS)
2. Electrónica de señal en rango de audio (herramientas maduras)
3. Algoritmo de procesado simple (correlación cruzada)
4. Validación en modelo físico (tubo de silicona + agua) posible en pocas semanas

**Version B requiere pasos adicionales:**
1. Diseño o compra de transductores de 50 kHz (nicho de mercado)
2. Caracterización de modos de onda guiada en arterias (o phantom)
3. Algoritmo de separación modal
4. La validación en phantom es más compleja (tubo de pared delgada con fluido)

---

## Recomendación Estratégica

### Estrategia de dos fases

**Fase A (inmediata): Desarrollar prototipo Version A**
- Más barato, más rápido, y el modelo matemático ya está listo
- Permite validar la hipótesis central (tomografía de red vascular) sin desarrollar nuevo hardware de RF
- Suficiente para demostración de concepto y publicación inicial

**Fase B (posterior): Desarrollar Version B como mejora**
- Mayor resolución y accesibilidad anatómica
- Diferenciador tecnológico para patente (combinación de ondas guiadas + tomografía de red en hemodinámica = vacío confirmado en literatura)
- Requiere investigación adicional sobre acoplamiento y modos

### Criterio de decisión entre versiones para un estudio clínico

| Escenario clínico | Versión recomendada |
|---|---|
| Enfermedad arterial periférica (miembros inferiores) | **A** (manguitos en muslo/pantorrilla) |
| Estenosis carotídea | **B** (acceso transcutáneo en cuello) |
| Estenosis renal o mesentérica | **B** (acceso abdominal) |
| Screening masivo en bajo coste | **A** (hardware más barato, operación sencilla) |
| Alta resolución espacial requerida | **B** (mejor BW relativo) |

---

## Estimación de Esfuerzo de Desarrollo

| Hito | Version A | Version B |
|---|---|---|
| Prototipo bench (validación en tubo) | 4–6 semanas | 8–12 semanas |
| Prototipo clínico (validación en sujeto sano) | 3–4 meses | 5–8 meses |
| Publicación de resultados de hardware | 6 meses | 10–12 meses |
| Solicitud de patente | Al completar prototipo bench | Al completar caracterización modal |

---

## Apéndice: Phantom de Validación para Ambas Versiones

Antes de probar en pacientes, ambas versiones deben validarse en un **phantom vascular físico**:

```
Especificación del phantom:
  - Tubo de silicona (Sylgard 184): pared 1.5 mm, diámetro interno 8 mm
  - Fluido: agua con glicerina al 10% (ρ ≈ 1060 kg/m³, viscosidad ≈ 3.5 mPa·s ≈ sangre)
  - Bifurcación en Y impresa en 3D (resina biocompatible)
  - Estenosis intercambiables: segmentos con diámetro reducido al 30/50/70%
  - Longitud total: 40 cm (madre) + 15 cm por rama hija
  - Soporte: placa de acrílico con fijaciones para manguitos (Version A) o transductores (Version B)
```

Este phantom permite:
- Verificar que los tiempos de tránsito medidos coinciden con el modelo (within 5%)
- Verificar que los ecos de estenosis son detectables con el SNR esperado
- Calibrar los parámetros de la cadena de señal antes de pasar a estudios in vivo
