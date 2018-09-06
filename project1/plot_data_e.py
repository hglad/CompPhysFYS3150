import numpy as np
import matplotlib.pyplot as plt

infile = open("project1_e_data.txt", 'r')
x = []
u = []
v = []
for line in infile:
    cols = line.split(' ')
    x.append(float(cols[0]))
    u.append(float(cols[1]))
    v.append(float(cols[2]))

#print (x)

plt.plot(x, v)
plt.plot(x, u)
plt.title('Task e) (ARMADILLO)')
plt.legend(['Numerical solution', 'Analytical solution'])
plt.grid('on')
plt.show()


"""infile = open("project1_c_data.txt", 'r')
x = []
u = []
v = []
for line in infile:
    cols = line.split(' ')
    x.append(float(cols[0]))
    u.append(float(cols[1]))
    v.append(float(cols[2]))

#print (x)

plt.plot(x, v)
plt.plot(x, u)
plt.title('Task c)')
plt.legend(['Numerical solution', 'Analytical solution'])
plt.grid('on')
plt.show()
"""
