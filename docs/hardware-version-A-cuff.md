# Hardware Version A — Manguito de Presión con Excitación Activa
## Vascular Network Tomography — Especificación de Ingeniería

**Versión:** 1.0  
**Principio:** Reflectometría Pneumática en Dominio del Tiempo (P-TDR)

---

## 1. Principio Físico

### Qué onda usamos

Esta versión trabaja con la **onda de Moens-Korteweg** — exactamente el mismo fenómeno que modela el simulador de Phase I/III. Es una onda de presión que viaja a lo largo de la columna de sangre dentro de la arteria, impulsada por la deformación elástica de la pared arterial.

```
c_MK = sqrt(E·h / (2·ρ·R))  ≈  4–12 m/s (arterias mayores)
```

Esta velocidad es muy baja comparada con el sonido convencional (~1540 m/s en tejido). Eso es una ventaja: con un pulso de presión aplicado externamente, la onda creada se propaga mucho más lejos antes de atenuarse.

### Cómo acopla el manguito con la arteria

El manguito inflado a presión sub-sistólica (~40–60 mmHg, menor que la presión sistólica del paciente) comprime levemente los tejidos circundantes pero no ocluye la arteria. Cuando se aplica una perturbación de presión controlada a través del manguito, esa perturbación se transmite mecánicamente a través de los tejidos blandos hasta la pared arterial, donde se convierte en una onda de presión intraluminal.

Este acoplamiento ya fue demostrado clínicamente: el dispositivo Arteriograph (Tensiomed, Hungría) usa este principio a presión supra-sistólica para medir la rigidez aórtica desde un solo manguito braquial. Nuestra variante usa **excitación activa** (señal controlada, no solo el pulso cardíaco) y **múltiples manguitos** para hacer tomografía de red.

---

## 2. Arquitectura del Sistema

### 2.1 Diagrama de bloques

```
[PC / DAQ]
   |
   |-- [Generador de señal (chirp 10–500 Hz)]
   |        |
   |        v
   |   [Amplificador electro-neumático]
   |        |
   |        v
   |   [Manguito A (fuente)]  -- arteria -- [Manguito B (sensor)] -- [Manguito C (sensor)]
   |        |
   |        |-- [Transductor de presión bajo manguito]
   |
   |-- [ADC 16-bit, 2 kHz]
   |
   v
[Software de adquisición + forward solver + solver inverso]
```

### 2.2 Número de manguitos

Para reproducir el experimento de validación de Phase III (Y-bifurcación):
- **Mínimo**: 3 manguitos (fuente + 2 sensores)
- **Árbol de 5 generaciones**: 7 manguitos (en los sitios anatómicos estándar)

Sitios anatómicos para el árbol de miembro inferior:
| Manguito | Arteria | Posición |
|---|---|---|
| A (fuente) | Aorta abdominal | Epigastrio |
| B | Ilíaca común derecha | Fosa ilíaca |
| C | Ilíaca común izquierda | Fosa ilíaca |
| D | Femoral derecha | Muslo proximal |
| E | Femoral izquierda | Muslo proximal |
| F | Poplítea derecha | Fosa poplítea |
| G | Poplítea izquierda | Fosa poplítea |

---

## 3. Especificación de Componentes

### 3.1 Manguito (bladder)

| Parámetro | Especificación |
|---|---|
| Tipo | Manguito de PVC médico, cámara de aire sellada |
| Tamaño | 12×22 cm (braquial), 14×32 cm (muslo) |
| Presión estática de trabajo | 40–80 mmHg |
| Presión máxima | 300 mmHg |
| Volumen de cámara | 50–150 mL |

### 3.2 Transductor electro-neumático (excitación)

Convierte la señal eléctrica de chirp en variación de presión dentro del manguito.

| Parámetro | Especificación |
|---|---|
| Tipo | Altavoz de bobina móvil acoplado a cámara sellada, ó válvula proporcional piezoeléctrica |
| Alternativa comercial | Proporcinador de presión SMC ITV-0010 (rango 0–100 kPa, BW 10 Hz) |
| Alternativa DIY | Tweeter de 25 mm + cámara de 20 mL sellada en T con manguito |
| Amplitud de señal inyectada | 5–15 mmHg pico |
| Ancho de banda | DC – 500 Hz |
| Linealidad | < 1% THD en rango de trabajo |

**Nota crítica**: El ancho de banda del sistema neumático está limitado por el volumen de la cámara de aire y la compliance del tubo. Para 500 Hz se requiere minimizar el volumen muerto (< 5 mL entre transductor y pared del manguito).

### 3.3 Sensor de presión (recepción)

| Parámetro | Especificación |
|---|---|
| Tipo | MEMS piezorresistivo o piezoeléctrico |
| Modelos comerciales | Honeywell HSCMRNN005PGSAX3 (0–5 PSI, 0.05% FS ruido), NXP MPX5050 |
| Sensibilidad | 4 mV/mmHg (típico para Honeywell HSCM) |
| Rango de presión | 0–150 mmHg |
| Ruido de fondo | 0.05 mmHg RMS a 1 kHz BW |
| Tiempo de respuesta | < 1 ms |
| Montaje | Dentro del manguito, en contacto directo con la piel o en el tubo de conexión |

### 3.4 Cadena de adquisición

| Componente | Especificación |
|---|---|
| Amplificador de instrumentación | INA128 (Gain = 100, ruido 8 nV/√Hz), ó INA333 para bajo consumo |
| Filtro anti-aliasing | Butterworth 4° orden, fc = 800 Hz |
| ADC | 16-bit, mínimo 2 kSa/s por canal (ADS1115 o similar) |
| Número de canales simultáneos | Igual al número de manguitos (3–7) |
| Latencia de adquisición | < 5 ms (para sincronización entre canales) |
| Sincronización | GPIO hardware trigger, jitter < 100 μs |

---

## 4. Señal de Excitación

### 4.1 Señal recomendada: chirp lineal

```
x(t) = A · sin(2π · (f0 + (f1-f0)/(2T) · t) · t)
donde:
  f0 = 10 Hz   (frecuencia inicial)
  f1 = 500 Hz  (frecuencia final)
  T  = 0.5 s   (duración del barrido)
  A  = 5–10 mmHg pico
```

Ventajas del chirp:
- Espectro plano → igual energía en todas las frecuencias
- Correlación cruzada con señal de referencia → mejora SNR por factor sqrt(T·BW) = sqrt(0.5 × 490) ≈ 15.6 dB
- Resolución temporal por correlación: Δt = 1/(2·BW) = 1/1000 = 1 ms → resolución espacial = c·Δt/2 = 5 × 0.001 / 2 = 2.5 mm

### 4.2 Alternativa: PRBS (Pseudo-Random Binary Sequence)

Secuencia binaria pseudoaleatoria de máxima longitud (m-sequence). Tiene espectro prácticamente plano hasta la frecuencia de reloj y propiedades de correlación óptimas. Usada en radar y sonar de bajo costo.

---

## 5. Presupuesto de Ruido y SNR

### Señal transmitida

| Variable | Valor |
|---|---|
| Amplitud del pulso inyectado | 10 mmHg pico |
| Amplitud en sensor fuente (reflexión de unión, Γ ≈ 0.35) | ~3.5 mmHg |
| Ruido del sensor (Honeywell HSCM, 1 kHz BW) | 0.05 mmHg RMS |
| SNR señal directa (sin promediado) | 46 dB |

### Eco de estenosis

| Variable | Valor |
|---|---|
| Coeficiente de reflexión en estenosis 70% (calculado Phase I) | Γ_sten ≈ 0.35 |
| Amplitud del eco TDR esperada | ~1.2 mmHg |
| SNR del eco (sin promediado) | 27.6 dB |
| SNR con 100 promedios coherentes (mejora 20 dB) | **47.6 dB** ✓ |
| SNR mínimo aceptable para detección confiable | 30 dB |

**Conclusión**: SNR suficiente con 100 promedios coherentes (≈ 50 s de medición a 2 latidos/s de excitación).

### Fuentes de interferencia

| Fuente | Frecuencia | Estrategia de mitigación |
|---|---|---|
| Pulso cardíaco | 0.8–2 Hz y armónicos hasta 20 Hz | Excitar por encima de 20 Hz; sustracción de línea de base |
| Respiración | 0.15–0.4 Hz | Promediado coherente sincronizado con ECG |
| Movimiento del paciente | Broadband | Detección de artefactos por varianza; rechazo de latidos contaminados |
| Ruido EM (50/60 Hz) | 50/60 Hz y armónicos | Filtro notch; apantallamiento del amplificador |

---

## 6. Resolución Espacial

La resolución con la que se puede localizar una estenosis depende del ancho de banda de la señal:

```
Δx = c · (1 / (2·BW)) = 5 m/s · (1 / (2·500 Hz)) = 5 mm
```

Con BW = 500 Hz y c = 5 m/s → resolución de **5 mm**. Suficiente para localizar estenosis en segmentos de 7–35 cm.

---

## 7. Protocolo de Medición

1. Colocar manguitos en sitios anatómicos definidos
2. Inflar a presión estática de acoplamiento (40–60 mmHg, sub-diastólica)
3. Registrar línea de base de 10 s (solo pulso cardíaco, sin excitación activa)
4. Inyectar señal de chirp desde manguito fuente, registrar todos los sensores simultáneamente
5. Repetir 100–200 veces, promediar coherentemente
6. Aplicar correlación cruzada con señal de referencia
7. Extraer M_obs = [TOF, amplitudes] → alimentar al solver inverso

Duración total del estudio: ~10–15 minutos por extremidad.

---

## 8. Costo Estimado del Prototipo

| Componente | Cantidad | Costo unit. (USD) | Total |
|---|---|---|---|
| Manguitos médicos | 7 | 30 | 210 |
| Sensores Honeywell HSCM | 7 | 15 | 105 |
| Transductores electro-neumáticos (SMC ITV) | 1 | 180 | 180 |
| Amplificadores INA128 + filtros | 7 | 5 | 35 |
| ADC ADS1115 (4 ch, 16-bit) | 2 | 10 | 20 |
| Microcontrolador (Raspberry Pi 4 + HAT ADC) | 1 | 75 | 75 |
| Regulador de presión + compresor mini | 1 | 120 | 120 |
| Tubería, conectores, carcasa | — | — | 100 |
| **Total** | | | **~845 USD** |

---

## 9. Limitaciones Conocidas

1. **Acoplamiento variable**: la capa de tejido adiposo entre el manguito y la arteria atenúa y distorsiona la señal. Pacientes obesos tendrán peor SNR.
2. **Acceso anatómico limitado**: no es posible colocar manguitos en arteria carótida, coronarias, o aorta torácica con este método.
3. **Crosstalk entre manguitos adyacentes**: si dos manguitos están cerca, la señal puede propagarse a través de los tejidos sin pasar por la arteria (artefacto de tejido blando).
4. **Ancho de banda limitado por el sistema neumático**: superar 500 Hz requiere minimizar el volumen muerto y usar una cámara de aire muy compacta.
5. **No diferencia unívocamente modo de propagación**: la señal medida mezcla el modo arterial (Moens-Korteweg) con modos de tejido blando. Se requiere modelo del tejido para separar contribuciones.
