"""
Programa para calcular la presión de saturación y entalpía de vaporización

"""

import math
import numpy as np

def presion_entalpia_sat(Psup, T, Tc, Pc):
  """
  Calcula la presión de saturación y entalpía de vaporización
  mediante la ec. de Van Der Wals y el factor de compresibilidad Z.

  Parámetros:
    Psup: Presión supuesta menor a la crítica (float).
    Tc: Temperatura crítica en Kelvin (float).
    Pc: Presión crítica en bar (float).
    T: Temperatura en Kelvin (float).

  Retorno:
    Presión de saturación en bar (float).
  """

  Tr = T / Tc
  Pr = Psup / Pc
  A = (27 * Pr) / (64 * (Tr**2))
  B = Pr / (8 * Tr)

  def calcular_z(A, B):
    b = -(B + 1)
    d = -A * B
    p = np.poly1d([1, b, A, d])
    raices = np.roots(p)
    return raices.tolist()

  z = calcular_z(A, B)
  zliq = min(z)
  zvap = max(z)

  def fi(z, A, B):
    return math.exp(z - 1 - math.log(z - B) - (A / z))

  fi_liq = fi(zliq, A, B)
  fi_vap = fi(zvap, A, B)

  Pi = Psup

  epsilon = 1e-6  # Tolerancia para la convergencia
  while abs(fi_vap - fi_liq) > epsilon:
    Pi_1 = Pi * (fi_liq / fi_vap)
    Pr = Pi_1 / Pc
    A = (27 * Pr) / (64 * (Tr**2))
    B = Pr / (8 * Tr)
    z = calcular_z(A, B)
    zliq = min(z)
    zvap = max(z)
    fi_liq = fi(zliq, A, B)
    fi_vap = fi(zvap, A, B)
    Pi = Pi_1
    # Cálculo de la entalpía residual
    h_liq_res = R * T * (zliq - 1 - A / zliq)
    h_vap_res = R * T * (zvap - 1 - A / zvap)
    # Cálculo de la entalpía de vaporización
    h_fusion = h_vap_res - h_liq_res

  return Pi, h_liq_res, h_vap_res, h_fusion

# Ejemplo de uso
Psup = 30  # Presión supuesta (bar)
T = 240  # Temperatura (K)
Tc = 304.2  # Temperatura crítica (K)
Pc = 73.83  # Presión crítica (bar)
R = 8.314  # Constante de los gases ideales (J/mol*K)

P_sat, h_liq_res, h_vap_res, h_fusion = presion_entalpia_sat(Psup, T, Tc, Pc)

print(f"La presión de saturación es de {P_sat:.4f} bar")
print(f"La entalpía residual del líquido es de {h_liq_res:.4f} J/mol")
print(f"La entalpía residual del vapor es de {h_vap_res:.4f} J/mol")
print(f"La entalpía de fusión es de {h_fusion:.4f} J/mol")