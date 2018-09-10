import numpy as np
import matplotlib.pyplot as plt

# Open file containing data for task e), read lines and store data in arrays
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
