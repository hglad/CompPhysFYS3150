import numpy as np
import matplotlib.pyplot as plt

infile = open("project1_d_data.txt", 'r')
log10_h = []
eps = []

for line in infile:
    cols = line.split(' ')
    log10_h.append(float(cols[0]))
    eps.append(float(cols[1]))

#print (x)
plt.plot(log10_h, eps)
plt.title('Task d)')
plt.grid('on')
plt.show()
