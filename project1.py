import numpy as np
import matplotlib.pyplot as plt

n = 100
x = np.linspace(0,1,n)
h = (x[-1] - x[0])/float(n)

# Vectors representing matrix A
a = np.ones(n-1)*(-1)
b = np.ones(n)*2
c = np.ones(n-1)*(-1)
#
v = np.zeros(n)
f = np.zeros(n) # Discretization of f(x)

def func_f(x):
    return h**2*100*np.exp(-10*x)

def u(x):       # Analytical solution
    return 1 - (1-np.exp(-10))*x - np.exp(-10*x)

# Row reduction using vectors only
"""
for i in range(0, n-1):
    a[i] -= a[i]/b[i-1] * b[i-1]
    b[i+1] -= a[i]/b[i] * c[i]
"""
f_tilde = np.zeros(n)
b_tilde = np.zeros(n)
f = func_f(x)
for i in range(1, n-1):
    f_tilde[i] = f[i] - a[i-1]/b[i-1] * f[i-1]
    b_tilde[i] = b[i] - a[i-1]*c[i-1]/b[i-1]

for i in range(1, n-1):
    i = n-1-i
    #print (i)
    v[i] = f_tilde[i]/(c[i-1]+b_tilde[i])

plt.plot(x, v)
plt.plot(x, u(x))
plt.legend(['Numerical', 'Analytical'])
plt.grid('on')
plt.show()
