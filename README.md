# Sistema de Monitoreo Fisiológico por Impedancia Bioeléctrica 🩺⚡

Este repositorio contiene el desarrollo analítico y numérico completo en Python para la calibración y el procesamiento de señales en un sensor portátil de impedancia bioeléctrica, simulado a una temperatura constante de 37 °C sobre muestras de tejido biológico.

## 📋 Estructura del Proyecto

El código implementa metodologías robustas de análisis numérico divididas en:
- **Análisis Exploratorio**: Visualización matemática del comportamiento en "U" asimétrico del tejido (Modelo Cole-Cole).
- **Interpolación avanzada**: Comparación empírica entre la interpolación polinómica local (Lagrange y Vandermonde), demostración de las oscilaciones numéricas del **Fenómeno de Runge** (con grado 29 global) y el ajuste óptimo mediante **Splines Cúbicos Naturales**.
- **Derivación Numérica**: Cálculo de la tasa de cambio y localización exacta del mínimo fisiológico local mediante cruce por cero ($f_{min} \approx 742.18$ Hz).
- **Algoritmos de Búsqueda de Raíces**: Implementación comparativa de convergencia entre los métodos de **Bisección** y **Newton-Raphson** para delimitar el ancho de banda seguro ($|Z| = 150$ Ω).
- **Análisis de Sensibilidad**: Evaluación de incertidumbres numéricas frente a ruidos experimentales instrumentales ($df/d|Z|$).

## 🚀 Instalación y Uso

1. Clona este repositorio:
   ```bash
   git clone [https://github.com/tu-usuario/sistema-monitoreo-impedancia.git](https://github.com/tu-usuario/sistema-monitoreo-impedancia.git)
   cd sistema-monitoreo-impedancia
