"""
Programa para calcular la temperatura de ebullición 
de una mezcla binaria usando la ecuación de Antoine.

"""

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

def temperatura_ebullicion_mezcla(P, x1, A1, B1, C1, A2, B2, C2):
  """
  Calcula la temperatura de ebullición de una mezcla binaria usando la ecuación de Antoine.

  Parámetros:
    P: Presión total del sistema en mmHg (float).
    x1: Fracción molar del primer componente (float).
    A1, B1, C1: Parámetros de la ecuación de Antoine para el primer componente (floats).
    A2, B2, C2: Parámetros de la ecuación de Antoine para el segundo componente (floats).

  Retorno:
    Temperatura de ebullición en Kelvin (float).
  """
  x2 = 1 - x1
  T1 = (B1/A1-math.log(P, 10))-C1  # Temperatura sat. sustancia 2
  T2 = (B2/A2-math.log(P, 10))-C2  # Temperatura sat. sustancia 1
  Tsup = x1*T1 + x2*T2

  P1 = presion_vapor_antoine(A1, B1, C1, Tsup)
  P2 = presion_vapor_antoine(A2, B2, C2, Tsup)

  if T1 > T2:
    Ti, Pi, Ai, Bi, Ci = T1, P1, A1, B1, C2
  else:
    Ti, Pi, Ai, Bi, Ci = T2, P2, A2, B2, C2

  volatilidad1 = P1/Pi
  volatilidad2 = P2/Pi
  Pcal = P/(x1*volatilidad1 + x2*volatilidad2)
  Tcal = (Bi/(Ai-math.log(Pcal, 10)))-Ci

  epsilon = 1**(-6)  # Tolerancia para la convergencia

  while abs(Tcal - Tsup) > epsilon:
    Tsup = Tcal

    P1 = presion_vapor_antoine(A1, B1, C1, Tsup)
    P2 = presion_vapor_antoine(A2, B2, C2, Tsup)

    if T1 > T2:
      Pi = P1
    else:
      Pi = P2

    volatilidad1 = P1/Pi
    volatilidad2 = P2/Pi 

    Pcal = P/(x1*volatilidad1 + x2*volatilidad2)
    Tcal = (Bi/(Ai-math.log(Pcal, 10)))-Ci
  return Tcal

def fraccion_molar_vapor(P, x1, A1, B1, C1, A2, B2, C2):
  """
  Calcula la fracción molar de una sustancia en la burbuja de vapor de una mezcla binaria.

  Parámetros:
    P: Presión total del sistema en mmHg (float).
    x1: Fracción molar del primer componente en el líquido (float).
    A1, B1, C1: Parámetros de la ecuación de Antoine para el primer componente (floats).
    A2, B2, C2: Parámetros de la ecuación de Antoine para el segundo componente (floats).

  Retorno:
    Fracción molar del componente en la burbuja de vapor (float).
  """
  x2 = 1 - x1
  p1 = presion_vapor_antoine(A1, B1, C1, T_ebullicion)
  p2 = presion_vapor_antoine(A2, B2, C2, T_ebullicion)
  return (x1*p1)/P

x1 = 0.5  # Fracción molar del primer componente
P = 760  # Presión total del sistema en mmHg

# Parámetros de Antoine para la sustancia 1
A1 = 7.31414
B1 = 1315.67
C1 = 240.479
# Parámetros de Antoine para la sustancia 2
A2 = 8.09126
B2 = 1582.91
C2 = 239.096

T_ebullicion = temperatura_ebullicion_mezcla(P, x1, A1, B1, C1, A2, B2, C2)
T_ebuliicion_Kelvin = T_ebullicion + 273.15
y1 = fraccion_molar_vapor(P, x1, A1, B1, C1, A2, B2, C2)
y2 = 1 - y1

print(f"Temperatura de ebullición: {T_ebullicion:.4f} °C    ||   {T_ebuliicion_Kelvin:.4f} K")
print(f"Fracción molar de la sustancia 1 en el vapor: {y1:.3f}")
print(f"Fracción molar de la sustancia 2 en el vapor: {y2:.3f}")
