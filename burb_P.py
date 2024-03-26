"""
Programa para calcular el punto burbuja de una mezcla usando la ecuación de Antoine		
burb_P: Calcula {y¡} y P, conocidas {x¡} y T			
"""
import numpy as np
import math

def presion_vapor_antoine(A, B, C, T):
  """
  Calcula la presión de vapor de un componente puro usando la ecuación de Antoine.

  Parámetros:
    A, B, C: Parámetros de la ecuación de Antoine para el componente (floats).
    T: Temperatura en Celsius (float).

  Retorno:
    Presión de vapor en mmHg (float).
  """
  P = 10**(A - B / (T + C))
  return P

def burb_P(T, mezcla):
  """
  Calcula la presión en el punto burbuja de una mezcla usando la ecuación de Antoine.

  Parámetros:
    T: Temperatura del sistema en °C (float).
    mezcla: Diccionario con los parametros de Antoine y fracciones molares

  Retorno:
    Presión total del sistema en mmHg (float).
  """

  Presiones = []
  x = [mezcla[componente]["x"] for componente in mezcla]
  A = [mezcla[componente]["A"] for componente in mezcla]
  B = [mezcla[componente]["B"] for componente in mezcla]
  C = [mezcla[componente]["C"] for componente in mezcla]
  nombre = [mezcla[componente]['nombre'] for componente in mezcla]

  n = len(x)
  for i in range(n):
    Pi = presion_vapor_antoine(A[i], B[i], C[i], T)
    Presiones.append(Pi)

  P_sat = np.dot(x, Presiones)

  return P_sat

def fraccion_vapor_P(T, mezcla, P_sat):
  """
  Calcula la fracción molar de una sustancia en la burbuja de vapor de una mezcla.

  Parámetros:
    T: Temperatura del sistema en °C (float).
    mezcla: Diccionario con los datos de cada componente

  Retorno:
    Fracción molar del componente en la burbuja de vapor (float).
  """

  y = []
  x = [mezcla[componente]["x"] for componente in mezcla]
  A = [mezcla[componente]["A"] for componente in mezcla]
  B = [mezcla[componente]["B"] for componente in mezcla]
  C = [mezcla[componente]["C"] for componente in mezcla]

  n = len(x)
  for i in range(n):
    Pi = presion_vapor_antoine(A[i], B[i], C[i], T)
    yi = (x[i]*Pi)/P_sat
    y.append(yi)
  return y

def seleccionar_compuestos(compuestos, seleccion):
    mezcla = {}
    for num in seleccion:
        compuesto = compuestos[num]
        x = float(input(f"Ingrese la fracción molar para {compuesto['nombre']}: "))
        mezcla[f"componente{num}"] = {
            "nombre": compuesto["nombre"],
            "A": compuesto["A"],
            "B": compuesto["B"],
            "C": compuesto["C"],
            "x": x
        }
    return mezcla

# Diccionario de compuestos numerados
compuestos = {
    1: {"nombre": "n-Hexano", "A": 6.88555, "B": 1175.817, "C": 224.867},
    2: {"nombre": "n-Heptano", "A": 6.90253, "B": 1267.828, "C": 216.823},
    3: {"nombre": "n-Octano", "A": 6.91874, "B": 1351.756, "C": 209.1},
    4: {"nombre": "n-Decano", "A": 6.95707, "B": 1503.568, "C": 194.738},
    5: {"nombre": "Acetona", "A": 7.31414, "B": 1315.67, "C": 240.479},
    6: {"nombre": "Agua", "A": 8.09126, "B": 1582.91, "C": 239.096},
    7: {"nombre": "n-Butano", "A": 6.82485, "B": 943.453, "C": 239.711}
}

# Selección de compuestos por número
seleccion = [1, 2, 3, 4, 7]  # Selecciona los compuestos 1, 2, 3, 4, 7
mezcla = seleccionar_compuestos(compuestos, seleccion)

"""Punto burbuja con temperatura conocida [Burb P]"""
T = 45.38084 # Temperatura total del sistema en °C

P_sat = burb_P(T, mezcla)
y = fraccion_vapor_P(T, mezcla, P_sat)

print("--------------------------------------------------")
print(f"Presión total del sistema: {P_sat:.4f} mmHg")
print("--------------------------------------------------")
print("Fracción de vapor cada componente")
for componente, valor_y in zip(mezcla.values(), y):
    print(f"{componente['nombre']}: {valor_y:.4f}")
print("--------------------------------------------------")