import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline
from metodos import (interpolacion_lagrange, matriz_vandermonde, 
                     evaluar_polinomio, biseccion, newton_raphson)

# =============================================================================
# DATOS EXPERIMENTALES (Tabla 1)
# =============================================================================
f = np.array([100, 120, 145, 170, 200, 235, 270, 310, 355, 405, 
              460, 520, 585, 655, 730, 810, 895, 985, 1080, 1180, 
              1290, 1410, 1540, 1680, 1830, 1990, 2160, 2340, 2530, 2730], dtype=float)

Z = np.array([152.3, 149.1, 146.8, 144.9, 142.0, 139.5, 137.9, 136.1, 134.8, 133.6, 
              132.7, 131.9, 131.4, 131.1, 130.9, 131.0, 131.3, 131.9, 132.7, 133.8, 
              135.2, 136.9, 138.9, 141.1, 143.5, 146.1, 149.0, 152.2, 155.6, 159.2], dtype=float)

print("=== PARTE A: ANÁLISIS EXPLORATORIO ===")
print(f"Número de muestras: {len(f)}")
print(f"Rango de frecuencias: {f.min()} Hz a {f.max()} Hz")
print(f"Mínimo visual directo en tabla: {Z.min()} Ω en f = 730 Hz\n")

# Gráfico Parte A
plt.figure(figsize=(10, 4.5))
plt.plot(f, Z, 'o-', color='#7c2d37', label='Datos medidos |Z|(f)')
plt.annotate('Mínimo ≈ (730 Hz, 130.9 Ω)', xy=(730, 130.9), xytext=(1200, 136),
             arrowprops=dict(facecolor='black', arrowstyle='->'))
plt.title('Parte A — Magnitud de Impedancia vs. Frecuencia')
plt.xlabel('Frecuencia f (Hz)')
plt.ylabel('Impedancia |Z| (Ω)')
plt.grid(True, linestyle=':', alpha=0.6)
plt.legend()
plt.savefig('grafico_parte_A.png', dpi=300)
plt.close()

# =============================================================================
# PARTE B1: INTERPOLACIÓN POLINÓMICA Y FENÓMENO DE RUNGE
# =============================================================================
print("=== PARTE B1: INTERPOLACIÓN POLINÓMICA ===")
# Demostración del fenómeno de Runge (Grado 29 global)
coeff_runge = np.polyfit(f, Z, 29)
f_fina = np.linspace(f.min(), f.max(), 1000)
Z_runge = np.polyval(coeff_runge, f_fina)

# Subconjunto de 8 puntos locales para evitar Runge (655 Hz a 1290 Hz)
idx_locales = (f >= 655) & (f <= 1290)
f_local = f[idx_locales]
Z_local = Z[idx_locales]

# Interpolación en f = 1000 Hz
Z_1000_lagrange = interpolacion_lagrange(f_local, Z_local, 1000)
coeff_vander = matriz_vandermonde(f_local, Z_local)
Z_1000_vander = evaluar_polinomio(coeff_vander, 1000)

print(f"Z(1000 Hz) por Lagrange (Local): {Z_1000_lagrange:.4f} Ω")
print(f"Z(1000 Hz) por Vandermonde (Local): {Z_1000_vander:.4f} Ω")

# Validación Leave-One-Out (LOO) con los 5 puntos específicos fijados en el informe
loo_points = [730, 1080, 655, 1290, 810]
errores_rel = []

print("\nTabla 2: Validación Leave-One-Out (LOO)")
print(f"{'f (Hz)':<10}{'Z_pred (Ω)':<15}{'Z_real (Ω)':<15}{'Error rel (%)':<15}")
for p in loo_points:
    # Excluir el punto actual del conjunto local
    f_loo = f_local[f_local != p]
    Z_loo = Z_local[f_local != p]
    
    z_pred = interpolacion_lagrange(f_loo, Z_loo, p)
    z_real = Z[f == p][0]
    err = abs(z_pred - z_real) / z_real * 100
    errores_rel.append(err)
    print(f"{p:<10}{z_pred:<15.4f}{z_real:<15.4f}{err:<15.4f}")

print(f"Error LOO promedio: {np.mean(errores_rel):.4f}%\n")

# Gráfico Fenómeno de Runge
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
ax1.plot(f, Z, 'ok', label='Datos')
ax1.plot(f_fina, Z_runge, '-', color='crimson', label='Polinomio Grado 29 (Runge)')
ax1.set_ylim(100, 300)
ax1.set_title('B1 — Grado 29: Oscilaciones de Runge')
ax1.set_xlabel('f (Hz)')
ax1.set_ylabel('|Z| (Ω)')
ax1.grid(True, alpha=0.4)
ax1.legend()

# Polinomios por tramos locales
ax2.plot(f, Z, 'ok', label='Datos')
f_fina_local = np.linspace(f_local.min(), f_local.max(), 200)
Z_local_eval = [interpolacion_lagrange(f_local, Z_local, xi) for xi in f_fina_local]
ax2.plot(f_fina_local, Z_local_eval, '-', color='navy', label='Polinomio Local Estables')
ax2.set_title('B1 — Polinomios Estables por Subconjunto Local')
ax2.set_xlabel('f (Hz)')
ax2.set_ylabel('|Z| (Ω)')
ax2.grid(True, alpha=0.4)
ax2.legend()
plt.savefig('grafico_parte_B1.png', dpi=300)
plt.close()

# =============================================================================
# PARTE B2 & C: SPLINE CÚBICO NATURAL Y DERIVACIÓN
# =============================================================================
print("=== PARTE B2 & C: SPLINE CÚBICO NATURAL Y DERIVACIÓN ===")
cs = CubicSpline(f, Z, bc_type='natural')
Z_1000_spline = cs(1000)
print(f"Z(1000 Hz) por Spline Cúbico Natural: {Z_1000_spline:.4f} Ω")

# Derivadas del Spline
cs_d1 = cs.derivative(1)
cs_d2 = cs.derivative(2)

# Encontrar el cero de la primera derivada en el intervalo del mínimo [655, 810]
f_min_calculado, _ = biseccion(cs_d1, 655, 810, tol=1e-6)
Z_min_calculado = cs(f_min_calculado)
d2_en_min = cs_d2(f_min_calculado)

print(f"Frecuencia exacta del mínimo encontrado: {f_min_calculado:.4f} Hz")
print(f"Impedancia mínima exacta calculada: {Z_min_calculado:.4f} Ω")
print(f"Segunda derivada en el mínimo: {d2_en_min:.6f} Ω/Hz² (Concavidad positiva)\n")

# Gráfico de Derivadas
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))
ax1.plot(f_fina, cs_d1(f_fina), color='royalblue', label='d|Z|/df (Primera Derivada)')
ax1.plot(f, cs_d1(f), 'or', alpha=0.5, label='Derivada en nodos')
ax1.axhline(0, color='black', linestyle='--', alpha=0.5)
ax1.axvline(f_min_calculado, color='red', linestyle=':', label=f'Mínimo f={f_min_calculado:.2f} Hz')
ax1.set_title('Parte C — Primera Derivada del Spline Cúbico')
ax1.set_ylabel('d|Z|/df (Ω/Hz)')
ax1.grid(True, alpha=0.4)
ax1.legend()

ax2.plot(f_fina, cs_d2(f_fina), color='forestgreen', label='d²|Z|/df² (Segunda Derivada)')
ax2.axhline(0, color='black', linestyle='--', alpha=0.5)
ax2.axvline(f_min_calculado, color='red', linestyle=':')
ax2.set_title('Parte C — Segunda Derivada del Spline Cúbico')
ax2.set_xlabel('Frecuencia f (Hz)')
ax2.set_ylabel('d²|Z|/df² (Ω/Hz²)')
ax2.grid(True, alpha=0.4)
ax2.legend()
plt.tight_layout()
plt.savefig('grafico_parte_C.png', dpi=300)
plt.close()

# =============================================================================
# PARTE D: BÚSQUEDA DE RAÍCES Y ANÁLISIS DE SENSIBILIDAD
# =============================================================================
print("=== PARTE D: BÚSQUEDA DE RAÍCES (|Z|(f) - 150 = 0) ===")
# Definición de la función objetivo y su derivada
def func_objetivo(x): return cs(x) - 150
def d_func_objetivo(x): return cs_d1(x)

# Raíz 1 (Baja frecuencia) en [100, 200]
f1_bis, iters_f1_bis = biseccion(func_objetivo, 100, 200, tol=1e-4)
f1_nr, iters_f1_nr = newton_raphson(func_objetivo, d_func_objetivo, 120, tol=1e-4)

# Raíz 2 (Alta frecuencia) en [2000, 2400]
f2_bis, iters_f2_bis = biseccion(func_objetivo, 2000, 2400, tol=1e-4)
f2_nr, iters_f2_nr = newton_raphson(func_objetivo, d_func_objetivo, 2200, tol=1e-4)

print("Tabla 4: Comparación de Métodos de Búsqueda de Raíces")
print(f"{'Raíz':<15}{'Bisección (Hz)':<18}{'Iters (Bis)':<15}{'Newton-R (Hz)':<18}{'Iters (NR)':<15}")
print(f"{'f1 (Baja)':<15}{f1_bis:<18.4f}{iters_f1_bis:<15}{f1_nr:<18.4f}{iters_f1_nr:<15}")
print(f"{'f2 (Alta)':<15}{f2_bis:<18.4f}{iters_f2_bis:<15}{f2_nr:<18.4f}{iters_f2_nr:<15}\n")

# Análisis de sensibilidad en f2
derivada_f2 = cs_d1(f2_nr)
sensibilidad_f2 = 1.0 / derivada_f2
print("Análisis de Sensibilidad en la Región Operativa Alta:")
print(f"d|Z|/df en f2: {derivada_f2:.4f} Ω/Hz")
print(f"Sensibilidad df/d|Z|: {sensibilidad_f2:.4f} Hz/Ω")
print(f"Una variación física de ±0.5 Ω altera la ubicación de la raíz f2 en ±{abs(sensibilidad_f2 * 0.5):.2f} Hz.")

# Gráfico Raíces Límites
plt.figure(figsize=(10, 5))
plt.plot(f_fina, cs(f_fina), color='blue', label='Spline Cúbico')
plt.axhline(150, color='red', linestyle='--', label='Z_th = 150.0 Ω')
plt.plot([f1_nr, f2_nr], [150, 150], 'X', color='orange', markersize=10, label='Raíces identificadas')
plt.title('Parte D — Identificación de Frecuencias Límite (Z_th = 150 Ω)')
plt.xlabel('f (Hz)')
plt.ylabel('|Z| (Ω)')
plt.grid(True, alpha=0.3)
plt.legend()
plt.savefig('grafico_parte_D.png', dpi=300)
plt.close()

print("\nSimulación completada con éxito. Gráficos exportados en la carpeta raíz.")
