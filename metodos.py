import numpy as np

def interpolacion_lagrange(x_nodos, y_nodos, x_eval):
    """
    Calcula el valor interpolado en x_eval usando el polinomio de Lagrange.
    """
    n = len(x_nodos)
    resultado = 0.0
    for i in range(n):
        termino = y_nodos[i]
        for j in range(n):
            if i != j:
                termino *= (x_eval - x_nodos[j]) / (x_nodos[i] - x_nodos[j])
        resultado += termino
    return resultado

def matriz_vandermonde(x_nodos, y_nodos):
    """
    Encuentra los coeficientes del polinomio resolviendo el sistema matricial V * c = y.
    """
    V = np.vander(x_nodos, increasing=True)
    coeficientes = np.linalg.solve(V, y_nodos)
    return coeficientes

def evaluar_polinomio(coeficientes, x_eval):
    """
    Evalúa un polinomio dado por sus coeficientes (orden creciente de grado).
    """
    resultado = 0.0
    for i, c in enumerate(coeficientes):
        resultado += c * (x_eval ** i)
    return resultado

def biseccion(func, a, b, tol=1e-5, max_iter=100):
    """
    Encuentra la raíz de una función en el intervalo [a, b] mediante Bisección.
    """
    if func(a) * func(b) >= 0:
        raise ValueError("El intervalo no contiene un cambio de signo (func(a)*func(b) >= 0).")
    
    iters = 0
    while (b - a) / 2.0 > tol and iters < max_iter:
        mid = (a + b) / 2.0
        iters += 1
        if func(mid) == 0:
            return mid, iters
        elif func(a) * func(mid) < 0:
            b = mid
        else:
            a = mid
    return (a + b) / 2.0, iters

def newton_raphson(func, dfunc, x0, tol=1e-5, max_iter=100):
    """
    Encuentra la raíz de una función usando el método de Newton-Raphson.
    """
    x = x0
    for iters in range(1, max_iter + 1):
        fx = func(x)
        dfx = dfunc(x)
        if abs(dfx) < 1e-12:
            raise ZeroDivisionError("Derivada cercana a cero en Newton-Raphson.")
        
        x_new = x - fx / dfx
        if abs(x_new - x) < tol:
            return x_new, iters
        x = x_new
    return x, max_iter
