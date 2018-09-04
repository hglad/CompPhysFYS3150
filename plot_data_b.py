import numpy as np
import matplotlib.pyplot as plt

infile = open("project1_b_data.txt", 'r')
x = []
u = []
v = []
for line in infile:
    cols = line.split(' ')
    x.append(float(cols[0]))
    u.append(float(cols[1]))
    v.append(float(cols[2]))

print (x)

plt.plot(x, v)
plt.plot(x, u)
plt.grid('on')
plt.show()
