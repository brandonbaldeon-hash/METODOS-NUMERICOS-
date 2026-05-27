# Análisis Numérico de Impedancia Bioeléctrica para Monitoreo Fisiológico de Tejido 🩺📊

Este repositorio contiene el desarrollo e implementación computacional de herramientas de análisis numérico (Interpolación, Diferenciación y Métodos de Búsqueda de Raíces) aplicadas a un sistema portátil de monitoreo fisiológico. El objetivo principal es caracterizar la impedancia bioeléctrica de un tejido bajo condiciones controladas (37 °C) dentro de un espectro de frecuencias de 100 Hz a 2730 Hz.

## 📊 Estructura y Respuestas Técnicas del Informe

### 🔹 Parte A y B: Interpolación y Modelado del Tejido
Para la reconstrucción continua de la señal de impedancia a partir de las 30 mediciones experimentales se evaluaron y compararon tres enfoques metodológicos:
1. **Interpolación Matricial (Vandermonde):** Planteamiento del sistema lineal $V \cdot c = |Z|$. Presenta limitaciones severas de estabilidad (mal condicionamiento) si se intenta un ajuste polinomial global de alto grado.
2. **Polinomios de Lagrange:** Implementación matemática precisa para interpolaciones locales eficientes.
3. **Splines Cúbicos Naturales:** Modelo óptimo elegido debido a su continuidad garantizada hasta la segunda derivada ($C^2$), eliminando el fenómeno de Runge y suavizando las transiciones en zonas de alta curvatura.

---

### 🔹 Parte C y D: Identificación de la Banda de Operación Segura ($|Z| = 150\ \Omega$)
El módulo de transmisión inalámbrica sufre atenuación crítica si la magnitud de la impedancia supera el umbral $Z_{th} = 150\ \Omega$. Se programaron algoritmos para resolver la ecuación:

$$\ f(|Z|) = |Z|(f) - 150 = 0 \$$

#### Comparativa de Métodos de Búsqueda de Raíces:
* **Método de Bisección:** Destaca por su **robustez absoluta**. Al garantizar el cumplimiento del Teorema del Valor Intermedio (cambio de signo en el intervalo $[100, 2730]\text{ Hz}$), asegura la convergencia matemática hacia las raíces con un error controlado, requiriendo más iteraciones pero con riesgo cero de divergencia.
* **Método de Newton-Raphson:** Implementado con una aproximación inicial cercana a los puntos de cruce. Destaca por su **convergencia cuadrática**, localizando las raíces con al menos 4 cifras significativas en una fracción de los pasos requeridos por bisección, condicionado a la disponibilidad de la derivada del spline.

#### Análisis de Sensibilidad ($f \approx 2000\text{ Hz}$):
Se calculó numéricamente la derivada inversa (sensibilidad de la ubicación de la raíz frente a ruidos en la medición):

$$\ \frac{df}{d|Z|} = \frac{1}{\frac{d|Z|}{df}} \$$

*Conclusión del análisis:* Pequeñas variaciones o ruidos de truncamiento en el instrumento de lectura de impedancia modifican de forma directa la ventana de frecuencia estimada para la banda segura, justificando la necesidad de aplicar filtrado digital.

---

### 🔹 Parte E: Implicaciones en Ingeniería y Aplicaciones Técnicas

* **Ingeniería Biomédica:** La identificación del mínimo de impedancia permite interpretar el estado fisiológico y la viabilidad del tejido celular (detección de isquemia o cambios en la perfusión local).
* **Diseño Electrónico:** Los resultados numéricos definen las especificaciones de diseño para los filtros analógicos de calibración, ajustando las frecuencias de corte en los puntos donde la impedancia $|Z|$ experimenta tasas de cambio rápidas.
* **Telecomunicaciones / Transmisión:** Delimita con precisión el rango de frecuencias seguras para la transmisión de datos del biosensor implantable hacia el receptor exterior sin pérdida de potencia.

## 🛠️ Propuestas de Mejora en la Adquisición de Datos
1. **Muestreo Adaptativo Densificado:** Implementar un paso de frecuencia más fino exclusivamente en las regiones de alta curvatura (cerca de los puntos de inflexión y del mínimo global) para maximizar la precisión de las derivadas numéricas.
2. **Control Estricto de Temperatura:** Minimizar la deriva térmica del biosensor para evitar desplazamientos artificiales en las curvas de impedancia durante el ensayo.

## 🚀 Requisitos del Entorno
* Python 3.x
* NumPy
* SciPy
* Matplotlib

---
**Desarrollado por:** Brandon Baldeón Ore (Código: 20190249) — Facultad de Ingeniería Electrónica y Eléctrica, Universidad Nacional Mayor de San Marcos.
**Docente:** Josué Alfonso Miranda Fernández.
