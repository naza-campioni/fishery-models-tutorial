"""
vpa_model.py

Core functions for Virtual Population Analysis (VPA) using Pope's approximation and Baranov's catch equation.

"""

import numpy as np

def calculate_terminal(M, F_AY, C_ay):
  """
  Estimates terminal population-at-age (survivors) from Baranov's catch equation

  Parameters:
    M (float) - natural mortality rate (assumed costant across year and age)
    F_AY (float) - assumed terminal fishing mortality
    C_ay (np.array) - catch-at-age matrix with shape [ages, years] (data)

  Returns:
    N_AY (np.array) - estimated terminal population-at-age, shape [ages]
  """
  if C_ay.ndim != 2:
    raise ValueError("C_ay must be a 2D array with shape [ages, years]")
    
  den = 1 - np.exp(-(M + F_AY))
  terminal_catch = C_ay[:,-1]
  N_AY = (terminal_catch / den) * ((F_AY + M) / F_AY)

  return N_AY




def calculate_N(M, C_ay, N_AY):
  """
  Estimates population-at-age matrix using Pope's approximation

  Parameters:
    M (float) - natural mortality rate (assumed costant across year and age)
    C_ay (np.array) - catch-at-age matrix with shape [ages, years] (data)
    N_AY (np.array) - estimated terminal population-at-age, shape [ages]

  Returns:
    N_ay (np.array) - estimated population-at-age, shape [ages, years]
  """

  a, y = C_ay.shape
  N_ay = np.zeros([a,y])
  N_ay[:,-1] = N_AY.flatten() # survivors
  
  # the last age group is tricky - it'd need a_MAX + 1 which doesn't
  # exist, so we assume all age classes beyond a_MAX are just a_MAX
  
  for i in reversed(range(y-1)):
    for j in reversed(range(a)):
      if j == a-1: # if at a_MAX
        N_ay[j,i] = N_ay[j,i+1]*np.exp(M) + C_ay[j,i]*np.exp(M/2)
      else:
        N_ay[j,i] = N_ay[j+1,i+1]*np.exp(M) + C_ay[j,i]*np.exp(M/2)




def calculate_N_vectorized(M, C_ay, N_AY):
  """
  Estimates population-at-age matrix using Pope's approximation but in
  a vectorized manner - just a different approach but same logic

  Parameters:
    M (float) - natural mortality rate (assumed costant across year and age)
    C_ay (np.array) - catch-at-age matrix with shape [ages, years] (data)
    N_AY (np.array) - estimated terminal population-at-age, shape [ages]

  Returns:
    N_ay (np.array) - estimated population-at-age, shape [ages, years]
  """

    a, y = C_ay.shape
    N_ay = np.zeros([a,y])
    N_ay[:,-1] = N_AY.flatten() # survivors
  
    # the last age group is tricky - it'd need a_MAX + 1 which doesn't
    # exist, so we assume all age classes beyond a_MAX are just a_MAX
    
    for i in reversed(range(y-1)):
      # calculate from penultimate age class backwards
      N_ay[:-1,i] = N_ay[1:,i+1] * np.exp(M) + C_ay[:-1,i].T * np.exp(M/2)
    
      # then at the last age we'd need N_ay[a_MAX + 1, y] so we just use N_ay[a_MAX]
      N_ay[-1,i] = N_ay[-1, i+1] * np.exp(M) + C_ay[-1,i].T * np.exp(M/2)
  

  

  
  
