"""
vpa_model.py

Core functions for Virtual Population Analysis (VPA) using Pope's approximation and Baranov's catch equation.

"""

import numpy as np

def calculate_terminal(M, F_AY, C_ay):
  """
  Estimates terminal population-at-age (survivors) from Baranov's catch equation

  Parameters:
    M - natural mortality rate (assumed costant across year and age)
    F_AY - assumed terminal fishing mortality
    C_ay - catch-at-age matrix with shape [ages, years] (data)
  """
  den = 1 - np.exp(-(M + F_AY))
  terminal_catch = C_ay[:,-1]
  N_AY = (terminal_catch / den) * ((F_AY + M) / F_AY)

  return N_AY
  

  
  
