import numpy as np

T = 1
E_ = -np.sinh(8)*8 / (3 + np.cosh(8))
Cv_ = 64/(T*T) * (np.cosh(8)*(3 + np.cosh(8)) - (np.sinh(8)**2))/(3 + np.cosh(8))**2
M2_ = (8*np.exp(8) + 8)/(3 + np.cosh(8))
absM_ = (2*np.exp(8) + 4)/(3 + np.cosh(8))
chi_ = M2_/T
print E_, Cv_, M2_, absM_, chi_
