import numpy as np
import matplotlib.pyplot as plt

n = 10
x = np.linspace(0,1,n)
h = x[-1]/float(n-1)
print (x, h)

# Vectors representing matrix A
a = np.ones(n-1)*(-1)
b = np.ones(n)*2
c = np.ones(n-1)*(-1)
#
v = np.zeros(n)
f_mark = np.zeros(n) # Discretization of f(x)

def func_f(x):  #f_mark (b_tilde in project text)
    return h**2*100*np.exp(-10*x)

def u(x):       # Analytical solution
    return 1 - (1-np.exp(-10))*x - np.exp(-10*x)

f_mark = func_f(x)

"""
Set f_tilde and b_tilde to f and b respectively, so that the first step can be calculated properly (since the first step does not use tilde values).
"""
f_tilde = f_mark
b_tilde = b
for i in range(1, n):
    b_tilde[i] = b[i] - a[i-1]/b_tilde[i-1] * c[i-1]
    f_tilde[i] = f_mark[i] - a[i-1]/b_tilde[i-1] * f_tilde[i-1]

for i in range(n-2, 0, -1):
    v[i] = (f_tilde[i] - v[i+1]*c[i])/b_tilde[i]

plt.plot(x, v)
plt.plot(x, u(x))
plt.legend(['Numerical', 'Analytical'])
plt.grid('on')
plt.show()
