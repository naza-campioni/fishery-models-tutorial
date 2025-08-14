import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize

def calculate_initial_N(s, F, M, C):
  """
  Calculates recruitments estimate from the catch data from the Baranov's equation.
  """
  Z = s[0]*F + M[0,:]
  num = C[0,:]*Z
  den = s[0]*F*(1 - np.exp(-Z))
  return num/den


def calculate_sa(a, a_50, a_95):
  """
  Assumes a logistic curve for the selectivity-at-age.
  """
  s_a = np.zeros(a)
  for i in range(a):
    s_a[i] = 1 / (1 + np.exp(-np.log(19) * (i-a_50) / (a_95-a_50)))
  return s_a


def calculate_N(a, y, s, F, M, N_init):
  """
  This function projects forward the numbers-at-age using the usual dynamics population equations.
  The special case a = a_max is treated following Equation 11.3, p.335 of 
  Haddon, Malcom, Modelling and Quantitative Methods in Fisheries, Chapman&Hall/CRC, 2001.
  """
  N = np.zeros((a,y))
  N[0,:] = N_init
  for age in range(a-1):
    if age+1 == a-1:
      N[age+1,0] = N[age, 0]*np.exp(-M[age, 0])/(1 - np.exp(-M[age,0]))
    else:
      N[age+1,0] = N[age, 0]*np.exp(-M[age, 0])
    for year in range(y-1):
      if age+1 == a-1:
        N[age+1, year+1] = N[age, year]*np.exp(-(s[age]*F[year] + M[age, year])) + N[age+1, year]*np.exp(-(s[age+1]*F[year] + M[age+1, year])) 
      else:
        N[age+1, year+1] = N[age, year]*np.exp(-(s[age]*F[year] + M[age, year]))
  return N


def calculate_C_hat(s, F, M, N):
  """
  Baranov's equation.
  """
  F_ay = s[...,np.newaxis]*F[...,np.newaxis].T
  Z = F_ay + M
  return (F_ay/Z)*N*(1-np.exp(-Z))


def calculate_C_tot_hat(C):
  """
  Total catch in year.
  """
  return np.sum(C,axis=0)


def calculate_C_prop_hat(C):
  """
  Proportions of catches in one year.
  """
  return C/np.sum(C, axis=0, keepdims=True)


def calculate_I_hat(q, N, ts, Z):
  """
  There are different ways to calculate I_ay - here it's simply 
  I = q.N.exp(-ts.Z)
  """
  return q*N*np.exp(-ts*(Z))


def ll_I(I, I_hat, sigma_I):
  """
  Lognormal likelihood.
  """
  ll = -0.5*(((np.log(I) - np.log(I_hat))**2)/sigma_I**2 + np.log(2*np.pi*sigma_I**2))
  return np.sum(ll)


def ll_C_tot(C_tot, C_tot_hat, sigma_tot):
  """
  Lognormal likelihood.
  """
  ll = -0.5*(((np.log(C_tot) - np.log(C_tot_hat))**2)/sigma_tot**2 + np.log(2*np.pi*sigma_tot**2))
  return np.sum(ll)


def ll_C_prop(n_y, C_prop, C_prop_hat):
  """
  Multinomial likelihood.
  """
  counts_obs = C_prop * n_y
  ll = (counts_obs * np.log(C_prop_hat)).sum()
  return ll



def neg_ll(theta, data):
  """
  Objective/loss function, returning the negative log-likelihood of the model.

  Parameters:
  - theta: list of parameters to optimize
  - data: dictionary of data

  Returns:
  - negative log-likelihood
  """
  
  a = data['a']
  y = data['y']
  C_ay = data['C_ay']
  C_tot = data['C_tot']
  C_prop = data['C_prop']
  I_ay = data['I_ay']
  M_ay = data['M_ay']
  ts = data['ts']
  n_y = data['n_y']

  a_50 = theta[0]                         # logistic curve parameter for selectivity
  a_95 = theta[1]                         # logistic curve parameter for selectivity
  F_y = np.exp(theta[2:2+y])              # total yearly fishing mortality
  N_init = np.exp(theta[2+y:2+2*y])       # recruits
  q = np.exp(theta[-3])                   # cathcability coefficient
  s_tot = np.exp(theta[-2])               # total catch variance
  s_I = np.exp(theta[-1])                 # index variance
  
  s = calculate_sa(a, a_50, a_95)

  N_ay = calculate_N(a, y, s, F_y, M_ay, N_init)
  C_hat = calculate_C_hat(s, F_y, M_ay, N_ay)

  F_ay = s[...,np.newaxis]*F_y[...,np.newaxis].T
  Z = F_ay + M

  I_hat = calculate_I_hat(q, N_ay, ts, Z)
  C_tot_hat = calculate_C_tot_hat(C_hat)
  C_prop_hat = calculate_C_prop_hat(C_hat)

  ll_Ind = ll_I(I_ay, I_hat, s_I)
  ll_tot = ll_C_tot(C_tot, C_tot_hat, s_tot)
  # ll_prop = ll_C_prop(a, C_prop, C_prop_hat, s_prop)
  ll_prop = ll_C_prop(n_y, C_prop, C_prop_hat)

  ll_total = ll_Ind + ll_tot + ll_prop
  return -ll_total

def fit_sca(data, theta_init):
  result = minimize(neg_ll, theta_init, args=(data,), method='L-BFGS-B')
  return result



"""
The following functions calculate the continuation-ratio logit transform that is
sometimes used to transform the proportions in case the actual number of fish aged per year n_y, ie the effective
sample size, is unknown.
Once transformed, the transformed proportions can be modelled normally.
For the notebook we assumed a known n_y for the sake of simplicity - the optimizer might behave badly
if n_y is to be estimated.
What proper SCA models would do is define a weak prior on weakly identifiable parameters, but for the moment
we'll stick to a simpler case just to showcase the core ideas of SCA models.
"""





# def cr_logit(a, C):
#   logit = []
#   for i in range(a-1):   # a-1 because the last proportion is left out (r/1-r is not defined since r=1)
#     num = C[i]
#     den = np.sum(C[i:])

#     r = num/den
#     logit.append(np.log(r/(1-r)))
#   return np.asarray(logit)


# def convert(a, y, C):
#   C_ = np.zeros((a-1,y))    # a-1 because we leave last proportion out
#   for i in range(y):
#     C_[:,i] = cr_logit(a,C[:,i])
#   return C_


# def cr_inverse(logits):
#     logits = np.asarray(logits)
#     A_minus1 = len(logits)
#     A = A_minus1 + 1
#     p = np.zeros(A)

#     S = 1.0  # running sum from current age to end
#     for i in range(A_minus1):
#         r = 1 / (1 + np.exp(-logits[i]))  # logistic
#         p[i] = r * S
#         S = S * (1 - r)
#     p[-1] = S
#     return p

# def inv_convert(a, y, logits):
#   log_ = np.zeros((a,y))
#   for i in range(y):
#     log_[:,i] = cr_inverse(logits[:,i])
#   return log_

# def ll_C_prop(a, C_prop, C_prop_hat, sigma_prop):
#   C_ = convert(a, C_prop)
#   C_hat = convert(a, C_prop_hat)

#   ll = -0.5*(((C_ - C_hat)**2)/sigma_prop**2 + np.log(2*np.pi*sigma_prop**2))
#   return np.sum(ll)
