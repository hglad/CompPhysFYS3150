import numpy as np
import matplotlib.pyplot as plt

# Open file containing data for task b), read lines and store data in arrays
infile = open("project1_b_data*.txt", 'r')
x = []
u = []
v = []
for line in infile:
    cols = line.split(' ')
    x.append(float(cols[0]))
    u.append(float(cols[1]))
    v.append(float(cols[2]))

n = len(x)-2        # Do not include end points

plt.plot(x, v)
plt.plot(x, u)
plt.title("Solution for n=%g" %n)
plt.xlabel('x'); plt.ylabel('v(x) and f(x)')
plt.legend(['Numerical solution v(x)', 'Analytical solution f(x)'])
plt.grid('on')
plt.show()
