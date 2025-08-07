import numpy as np
from tqdm import tqdm

def calculate_terminal_N(M, F_AY, C_ay):
  """
  Estimates terminal population-at-age (survivors) from Baranov's catch equation
  For details check the VPA folder
  """
  if C_ay.ndim != 2:
    raise ValueError("C_ay must be a 2D array with shape [ages, years]")

  den = 1 - np.exp(-(M + F_AY))
  terminal_catch = C_ay[:,-1]
  N_AY = (terminal_catch / den) * ((F_AY + M) / F_AY)

  return N_AY

def calculate_N_VPA(M, C_ay, N_AY):
  """
  VPA estimate of N_ay - more details in the VPA folder
  """

  a, y = C_ay.shape
  N_ay = np.zeros([a,y])
  N_ay[:,-1] = N_AY.flatten() # survivors

  # the last age group is tricky - it'd need a_MAX + 1 which doesn't
  # exist, so we assume all age classes beyond a_MAX are just a_MAX -> plus group

  for age in reversed(range(a)):
    for year in reversed(range(y-1)):
      if age == a-1:
        N_ay[age, year] = N_ay[age, year +1]*np.exp(M) + C_ay[age,year]*np.exp(M/2)
      else:
        N_ay[age, year] = N_ay[age +1, year +1]*np.exp(M) + C_ay[age,year]*np.exp(M/2)
  return N_ay

def calculate_F_ay(a, y, N_ay, M, F_AY):
  """
  Population dynamics equation solved for F_ay
  Clips low fishing mortality values to a minimum of 0.2.
  This avoids runaway growth in N_ay estimates due to weak exponential decay
  when Z = F + M is too small.
  """
  F_ay = np.zeros((a,y))
  for age in range(a-1):
    for year in range(y-1):
      num = N_ay[age, year]
      den = N_ay[age+1, year+1]
      if num > 0 and den > 0:
          F_ay[age,year] = max(1e-6, np.log(num / den) - M) 
      else:
          F_ay[age,year] = F_AY  # fallback

  F_ay[-1, :-1] = F_ay[-2, :-1]  # plus group
  F_ay[:, -1] = F_AY             # terminal year

  F_ay = np.where(F_ay < 0.1, 0.2, F_ay) # clip to avoid blowup of population estimates
  F_ay = np.where(F_ay > 2.0, 2.0, F_ay)
  return F_ay

def calculate_ln_r(N_ay, u_ay, w):
  """
  Equation 19
  """
  ln_ratio = np.log(N_ay/u_ay)
  num_r = np.sum(w*ln_ratio, axis=1) # sum over years
  den_r = np.sum(w, axis=1) # sum over years
  
  ln_r = num_r/den_r     # equation 19
  return ln_r 

def calculate_adjusted_weights(a, y, w, F_ay):
  """
  Weights adjustment as per after equation 22
  """
  ECF = np.zeros((a,y))
  for age in range(a):
    ECF[age,:] = np.exp(np.sum(F_ay[age:,:], axis=0))  # sum over ages    
  w_1 = w/ECF    
  return w_1

def calculate_cumZ(a, y, F_ay, M_ay):
  cumZ = np.zeros((a,y))
  Z_ay = F_ay + M_ay
  for age in range(a):
    cumZ[age,:] = np.sum(Z_ay[age:,:], axis=0)
  return cumZ


def estimate_Pk(ln_r, u_ay, w_1, cumZ):
  """
  Equation 22
  """
  num = np.sum(w_1 * (np.log(u_ay) + ln_r[:, np.newaxis] - cumZ), axis=0)
  den = np.sum(w_1, axis=0)
  ln_Pk = num / den
  return np.exp(ln_Pk)

def diagonal_N_ay(a, y, Pk, ECM, C_ay, M_ay):
  """
  Reconstruction of diagonal of N_ay through P_t(k)
  """
  N_ay = np.zeros((a, y))
  for age in range(a):
    for year in range(y):
      k = year - age
      if 0 <= k < len(Pk):
        N_ = Pk[k] * ECM[age, year]
        P_c = 0
        for i in range(age, a):
          if 0 <= k + i < y:
            P_c += ECM[i, k + i] * C_ay[i, k + i] * np.exp(-0.5 * M_ay[i, k + i])
            N_ay[age, year] = N_ + P_c
      else:
        N_ay[age,year] = 0

  return N_ay
      
def fill_last_age_zeros(a, y, N_ay, C_ay, M):
  """
  Fill in the zero entries from the last age through VPA
  """
  for year in reversed(range(y-(y-a)-1)):
    N_ay[-1, year] = N_ay[-1, year +1]*np.exp(M) + C_ay[-1,year]*np.exp(M/2)
    
  return N_ay

def fill_remaining_zeros(a, y, N_ay, C_ay, M):
  """
  Fill in the remaining zero entries through VPA
  """
  for age in reversed(range(1,a-1)):
    for year in reversed(range(y-(y-a)-1)):
      if N_ay[age,year] == 0 and N_ay[age+1,year+1]>0:
        N_ay[age,year] = N_ay[age +1, year +1]*np.exp(M) + C_ay[age,year]*np.exp(M/2)
        
  return N_ay

def reconstruct_Nay(a, y, Pk, ECM, M_ay, C_ay):
  N_ay = diagonal_N_ay(a, y, Pk, ECM, C_ay, M_ay)
  N_ay = fill_last_age_zeros(a, y, N_ay, C_ay, M_ay[0,0])  # M_ay is constant
  N_ay = fill_remaining_zeros(a, y, N_ay, C_ay, M_ay[0,0])
  return N_ay

def main_XSA(a, y, C_ay, u_ay, M, F_AY, w, iterations=5):
  """
    Main XSA function implementing Shepherd's Extended Survivors Analysis algorithm.

    Parameters:
    - C_ay: np.array, catch-at-age [a x y]
    - u_ay: np.array, abundance index matrix [a x y]
    - w:    np.array, weights matrix [a x y]
    - M:    float, natural mortality
    - F_AY: float, terminal fishing mortality
    - iterations: int, number of XSA iterations

    Returns:
    - N_ay: np.array, reconstructed population at age matrix [a x y]
    - F_ay: np.array, fishing mortality at age matrix [a x y]
    - N_AY_: np.array, survivors P_t(k)
    """
  
  # calculate terminal population-at-age (survivors)
  N_AY = calculate_terminal_N(M, F_AY, C_ay)
  
  # calculate full population-at-age matrix using Pope's approximation
  N_ay = calculate_N_VPA(M, C_ay, N_AY)
  first_Nay = N_ay
  
  # initialize quantities
  M_ay = np.ones((a,y))*M     # natural mortality matrix - constant through ages&years
  F_ay = np.zeros((a,y))      # fishing mortality matrix
  ln_r = np.zeros(a)          # reciprocal of catchability coefficient
  ECF = np.zeros((a,y))       # exponential cum fishing mortality
  cumZ = np.zeros((a,y))      # cum mortality
  ECM = np.zeros((a,y))       # exponential cum natural mortality
  
  for age in range(a):
    ECM[age, :] = np.exp(np.sum(M_ay[age:, :], axis=0))    

  for it in tqdm(range(iterations)):
    F_ay = calculate_F_ay(a, y, N_ay, M, F_AY)
    cumZ = calculate_cumZ(a, y, F_ay, M_ay)
    ln_r = calculate_ln_r(N_ay, u_ay, w)
    w_1 = calculate_adjusted_weights(a, y, w, F_ay)    
    Pk = estimate_Pk(ln_r, u_ay, w_1, cumZ)
    N_ay = reconstruct_Nay(a, y, Pk, ECM, M_ay, C_ay)

  return N_ay, F_ay, Pk, first_Nay
    
  

