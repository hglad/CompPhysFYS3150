import numpy as np
import matplotlib.pyplot as plt

E, M = np.loadtxt("ising_data.txt", usecols=(0,1), unpack=True, dtype='float')
n = len(E)
T = 1
#E = np.loadtxt("ising_data.txt")
#print E
E2 = np.dot(E,E)
M2 = np.dot(M,M)

Cv = (E2 - E)/(T*T*n)
Chi = (M2 - M)/(T*n)

"""
for line in infile:
    cols = line.split(' ')
    E.append(float(cols[0]))
    M.append(float(cols[1]))
"""

x = np.linspace(0, len(E)-1, len(E))        # MC-cycles

plt.plot(x, E)
plt.grid('on'); plt.ylabel('E')

plt.figure()
plt.plot(x, M)
plt.grid('on'); plt.ylabel('M')

plt.figure()
plt.plot(x, Cv)
plt.grid('on'); plt.ylabel('Cv')

plt.figure()
plt.plot(x, Chi)
plt.grid('on'); plt.ylabel('Chi')

# Analytical values
E_ = -np.sinh(8)*8 / (3 + np.cosh(8))
Cv_ = T**(-2) * (64*np.cosh(8)*(3 + np.cosh(8)) - 8*np.sinh(8)*8*np.sinh(8)) / (3 + np.cosh(8))**2
M2_ = (8*np.exp(8) + 8)/(3 + np.cosh(8))
absM_ = (2*np.exp(8) + 4)/(3 + np.cosh(8))
chi_ = M2_/T

print E_, absM_, M2_, Cv_, chi_

plt.show()
