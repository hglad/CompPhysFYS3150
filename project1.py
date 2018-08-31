import numpy as np

n = 10
x = np.linspace(0,1,n)
h = (x[-1] - x[0])/float(n)

# Vectors representing matrix A
a = np.zeros(n)
b = np.zeros(n)
c = np.zeros(n)
#
v = np.zeros(n)
f = np.zeros(n) # Discretization of f(x)

def f(x):
    return 100*np.exp(-10*x)

def u(x):       # Analytical solution
    return 1 - (1-np.exp(-10))*x - np.exp(-10*x)

# Row reduction using vectors only
sub_1 = np.asarray([])
