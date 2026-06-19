# Hardware Version B — Ultrasonido de Baja Frecuencia (Ondas Guiadas)
## Vascular Network Tomography — Especificación de Ingeniería

**Versión:** 1.0  
**Principio:** Inspección por Ondas Guiadas Arteriales (Arterial Guided Wave Testing, A-GWT)

---

## 1. Principio Físico

### Qué onda usamos (y por qué es distinta a Version A)

Esta versión NO excita la onda de Moens-Korteweg (c ≈ 5 m/s). En cambio, excita **ondas guiadas de alta frecuencia** que se propagan a lo largo de la pared arterial a velocidades mucho mayores (100–2000 m/s, dependiendo del modo).

La arteria es físicamente un **cilindro de pared delgada** relleno de fluido. Cuando se aplica vibración mecánica a su exterior, la pared soporta modos de onda guiada análogos a las **ondas de Lamb** en placas delgadas. Estos modos combinan deformación longitudinal y flexural de la pared.

Este principio es **idéntico** al que usan los ingenieros para inspeccionar tuberías de petróleo y gas con ultrasonido (Guided Wave Testing, GWT o Long Range Ultrasound Testing, LRUT) — salvo que aquí las "tuberías" son arterias de 5–25 mm de diámetro.

### Por qué 20–100 kHz y no MHz

El ultrasonido médico convencional (2–15 MHz) está diseñado para crear **imágenes transversales** con alta resolución espacial. Eso requiere frecuencias altas → longitudes de onda pequeñas → resolución milimétrica en corte transversal.

Nosotros queremos propagar la onda **a lo largo del árbol arterial** durante 10–50 cm para detectar estenosis remotas. A frecuencias altas:
- La atenuación en tejido es 0.5 dB/cm/MHz
- A 5 MHz, en 30 cm de tejido: 0.5 × 5 × 30 = **75 dB de atenuación** → señal indetectable

A frecuencias bajas:
- A 50 kHz (0.05 MHz), en 30 cm: 0.5 × 0.05 × 30 = **0.75 dB** → prácticamente sin pérdida
- A 100 kHz: **1.5 dB** → excelente

La baja frecuencia permite **propagación de largo alcance**. La contrapartida es menor resolución espacial, que se compensa con el procesado de señal (chirp + correlación).

---

## 2. Modos de Onda Guiada en Arterias

Para una pared arterial de grosor h ≈ 1 mm y radio R ≈ 4 mm, los modos más relevantes son:

| Modo | Tipo | Velocidad aprox. | Sensibilidad a estenosis |
|---|---|---|---|
| L(0,1) | Longitudinal axisimétrico | 1200–1800 m/s | Alta (compresión axial) |
| F(1,1) | Flexural (curvado) | 100–500 m/s | Media |
| T(0,1) | Torsional | 700–1200 m/s | Baja |

El modo **L(0,1)** es el más útil: es el de mayor velocidad, el menos dispersivo en el rango 20–100 kHz, y el más sensible a cambios de sección transversal (que es exactamente lo que hace una estenosis).

La velocidad del modo L(0,1) depende de la rigidez de la pared. Una estenosis crea:
1. Reducción local del diámetro → cambio de impedancia acústica de la pared
2. Aumento del grosor relativo de la pared → cambio de rigidez local
3. Ambos efectos producen **reflexión parcial** de la onda guiada

---

## 3. Arquitectura del Sistema

### 3.1 Diagrama de bloques

```
[PC / DAQ]
   |
   |-- [Generador de función arbitraria (AWG)]  →  chirp 20–100 kHz, 10 ciclos
   |        |
   |        v
   |   [Amplificador de potencia  ~5–20 W]
   |        |
   |        v
   |   [Transductor Tx (PZT-4, 50 kHz)]
   |        | ← gel de acoplamiento
   |        | ← sobre la piel, encima de la arteria
   |        |
   |        ~~~~~  ondas guiadas viajan por la pared arterial  ~~~~~
   |        |
   |   [Transductor Rx1]  [Rx2]  [Rx3]  (a distancias conocidas)
   |        |
   |        v
   |   [Amplificador de recepción + filtro BP 20–100 kHz]
   |        |
   |        v
   |   [ADC 16-bit, 1–10 MSa/s]
   |
   v
[Software: correlación cruzada → TOF → solver inverso]
```

### 3.2 Configuración de transductores

**Modo pitch-catch (recomendado para tomografía de red):**
- Transductor Tx fijo en la arteria proximal (equivalente al nodo fuente en el grafo)
- Múltiples Rx distribuidos a lo largo del árbol arterial
- La señal parte de Tx, viaja por la pared arterial, llega a cada Rx con un retardo determinado por la distancia y la velocidad del modo

**Modo pulso-eco:**
- Un solo transductor transmite y recibe
- Detecta ecos reflejados desde discontinuidades (estenosis, bifurcaciones)
- Más simple en hardware, pero no da toda la información del modo pitch-catch

Para la tomografía de red arterial, se usará **pitch-catch con múltiples pares Tx-Rx** — equivalente exacto a la configuración de múltiples fuentes y sensores del solver de Phase I.

---

## 4. Especificación de Componentes

### 4.1 Transductor piezoeléctrico (Tx/Rx)

| Parámetro | Especificación |
|---|---|
| Material | PZT-4 (alta potencia) para Tx, PZT-5A (alta sensibilidad) para Rx |
| Geometría | Disco circular, diámetro 20–25 mm, grosor ≈ λ/2 a frecuencia de resonancia |
| Frecuencia de resonancia | 50 kHz ± 5% |
| Sensibilidad de recepción | ~30–50 mV/Pa a resonancia |
| Impedancia eléctrica | ~50–500 Ω a resonancia |
| Modelos comerciales | Steminc SMD50T21F50 (50 kHz), Steiner & Martins SMBLTD55F50H |
| Precio unitario | 15–40 USD |

### 4.2 Acoplamiento a la piel

| Parámetro | Especificación |
|---|---|
| Gel de acoplamiento | Gel de ultrasonido médico estándar (impedancia acústica ≈ 1.65 MRayl) |
| Impedancia del PZT | ~35 MRayl |
| Impedancia de tejido blando | ~1.6 MRayl |
| Pérdida por desacoplamiento en interfaz | ~13 dB (se recupera parcialmente con capa de adaptación) |
| Capa de adaptación opcional | Resina epoxy con carga metálica, Z_adapt = sqrt(Z_PZT × Z_tejido) ≈ 7.5 MRayl |
| Pérdida con capa adaptadora | ~4 dB (mejora de 9 dB respecto a sin adaptación) |
| Presión de contacto necesaria | 0.5–2 N/cm² (aplicable con vendaje elástico o parche adhesivo) |

### 4.3 Generador de función arbitraria (AWG)

| Parámetro | Especificación |
|---|---|
| Frecuencia de salida | DC – 200 kHz |
| Formas de onda | Chirp, burst, PRBS, arbitraria |
| Resolución de amplitud | 14-bit |
| Modelos comerciales | Siglent SDG1032X (30 MHz, 2 ch, ~350 USD), o tarjeta NI USB-6366 (2 MSa/s) |
| Alternativa de bajo costo | STM32 DAC a 12-bit + amplificador operacional → hasta 200 kHz (~50 USD en total) |

### 4.4 Amplificador de potencia (Tx)

| Parámetro | Especificación |
|---|---|
| Potencia de salida | 5–20 W RMS (0.5–2 W/cm² sobre el transductor de 20 mm) |
| Ancho de banda | 10 Hz – 500 kHz |
| Impedancia de salida | 50 Ω o adaptada al transductor |
| Modelo comercial | EPA-104 (Piezo Systems), o amplificador de audio clase D + transformador de acoplamiento |
| Alternativa DIY | TDA2003 + LC de adaptación (~8 USD) |

### 4.5 Cadena de recepción

| Componente | Especificación |
|---|---|
| Preamplificador | LNA de bajo ruido, NF < 2 dB, ganancia 30–40 dB |
| Modelo comercial | AD8332 (Analog Devices), NF 1.9 dB, programable |
| Filtro paso-banda | 20–100 kHz, orden 6, Butterworth o Chebyshev |
| ADC | 16-bit, mínimo 500 kSa/s por canal (ADS9256 o similar) |
| Número de canales | 3–7 canales simultáneos |
| Aislamiento Tx/Rx | Switch de RF o limiter diode (protección del receptor durante transmisión) |

---

## 5. Señal de Excitación

### 5.1 Burst de chirp (recomendado)

```
x(t) = A · w(t) · sin(2π · (f0 + B/(2·T) · t) · t)

donde:
  f0 = 20 kHz   (frecuencia inicial)
  f1 = 100 kHz  (frecuencia final)
  B  = 80 kHz   (ancho de banda)
  T  = 500 μs   (duración del chirp)
  w(t) = ventana Hanning (reduce lóbulos laterales)
  A  = ajustado para intensidad de salida ≤ 500 mW/cm²
```

La correlación cruzada de la señal recibida con el chirp de referencia produce una respuesta impulsional del sistema con resolución temporal:

```
Δt = 1 / (2·B) = 1 / (2 × 80 000) = 6.25 μs
```

### 5.2 Resolución espacial resultante

| Modo de onda | Velocidad c | Resolución Δx = c·Δt/2 |
|---|---|---|
| L(0,1) en arteria | 1500 m/s | **4.7 mm** |
| Tejido blando (eco en directo) | 1540 m/s | 4.8 mm |
| Modo de Moens-Korteweg (si se excita) | 5–10 m/s | 0.016–0.031 mm (submilimétrico, pero requiere frecuencias < 1 kHz para esta resolución) |

La resolución de **~5 mm** para el modo L(0,1) es comparable o mejor que la Version A para la misma duración de señal.

---

## 6. Presupuesto de Ruido y SNR

### Estimación de señal recibida

```
Intensidad transmitida:         I_tx = 100 mW/cm²
Área del transductor:           S = π × (10mm)² = 3.14 cm²
Potencia total:                 P = 314 mW

Atenuación en tejido (50 kHz, 10 cm de trayecto):
  α = 0.5 dB/cm/MHz × 0.05 MHz = 0.025 dB/cm
  A_tejido = 0.025 × 10 = 0.25 dB → factor 0.97 (atenuación casi nula)

Pérdida de acoplamiento con adaptador:  4 dB → factor 0.63
Pérdida geométrica (divergencia):       ~6 dB a 10 cm → factor 0.5

Presión recibida en Rx (modo directo):
  P_rx ≈ sqrt(2 × ρ × c × I_tx) × 0.97 × 0.63 × 0.5
       = sqrt(2 × 1060 × 1540 × 0.01) × 0.31    [I en W/cm² = 0.01 W/cm² → 100 mW/cm²]
       ≈ 180 Pa × 0.31 ≈ 56 Pa
```

### Señal de eco desde la estenosis

```
Coeficiente de reflexión del modo L(0,1) en estenosis 70%:  Γ ≈ 0.15–0.30
Trayecto adicional (2 × distancia a estenosis, típ. 2×15cm = 30cm):  ~7.5 dB adicional
  P_eco ≈ 56 Pa × 0.25 × 0.97^(15) ≈ 56 × 0.25 × 0.70 ≈ 9.8 Pa
```

### Ruido del receptor

```
Densidad espectral de ruido del AD8332:  NF = 1.9 dB, en 50Ω
  NSD = sqrt(4 k T R) = sqrt(4 × 1.38e-23 × 310 × 50) ≈ 0.92 nV/√Hz
  Con NF 1.9 dB (factor 1.24):  ~1.14 nV/√Hz

Ruido RMS en ancho de banda 80 kHz:
  V_ruido = 1.14e-9 × sqrt(80000) = 322 nV = 0.32 μV

Sensibilidad del transductor a 50 kHz:  40 mV/Pa
  V_señal_eco = 9.8 Pa × 0.040 V/Pa = 392 mV

SNR del eco (sin promediado):
  SNR = 20 log10(392e-3 / 0.32e-6) = 20 log10(1.225e6) ≈ 122 dB
```

En la práctica, el SNR estará limitado por el ruido de acoplamiento, artefactos de movimiento y reflexiones espurias de tejido (~40–60 dB límite real). Pero el margen teórico es muy amplio, lo que significa que hay tolerancia para pérdidas de acoplamiento y accesorios de bajo costo.

---

## 7. Protocolo de Medición

1. Aplicar gel de ultrasonido y colocar transductores sobre piel, encima de las arterias diana
2. Fijar con vendaje elástico a 1 N/cm² de presión de contacto
3. Verificar acoplamiento: señal directa Tx→Rx más próximo debe superar umbral SNR > 30 dB
4. Transmitir 200 chirps por par Tx-Rx, registrar simultáneamente en todos los Rx
5. Aplicar correlación cruzada, extraer respuesta impulsional
6. Identificar picos de TOF y amplitudes → vector M_obs
7. Ejecutar solver inverso

Duración total: ~5 minutos por árbol arterial.

---

## 8. Costo Estimado del Prototipo

| Componente | Cantidad | Costo unit. (USD) | Total |
|---|---|---|---|
| Transductores PZT 50 kHz | 7 | 30 | 210 |
| Preamplificadores AD8332 + PCB | 6 | 20 | 120 |
| Filtros BP 20–100 kHz | 6 | 10 | 60 |
| ADC ADS9256 (8ch, 16-bit, 1 MSa/s) | 1 | 80 | 80 |
| Generador AWG (Siglent SDG1032X) | 1 | 350 | 350 |
| Amplificador de potencia Tx | 1 | 60 | 60 |
| FPGA/SoC para sincronización (Zynq Eval board) | 1 | 200 | 200 |
| Gel de acoplamiento, cables, carcasa | — | — | 80 |
| **Total** | | | **~1160 USD** |

---

## 9. Limitaciones Conocidas

1. **Dispersión modal**: los modos de ondas guiadas tienen velocidades que dependen de la frecuencia. Esto complica la extracción de TOF preciso y requiere compensación en el procesado de señal.
2. **Múltiples modos simultáneos**: el transductor puede excitar varios modos al mismo tiempo; separar sus ecos requiere procesado avanzado (análisis tiempo-frecuencia, 2D FFT).
3. **Acoplamiento indirecto**: la onda viaja desde el transductor a través de ~1–5 cm de tejido antes de alcanzar la pared arterial. Las reflexiones espurias de interfases de tejido (músculo/grasa/hueso) pueden enmascarar los ecos arteriales.
4. **Seguridad acústica**: a 500 mW/cm², el índice mecánico (MI) debe verificarse para no superar los límites de seguridad FDA (MI < 1.9 para diagnóstico). A 50 kHz y 100 mW/cm², se está en zona segura pero debe documentarse.
5. **Variabilidad anatómica**: la profundidad de la arteria varía entre pacientes (1–5 cm), afectando el acoplamiento y la atenuación geométrica.
