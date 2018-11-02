import numpy as np
import matplotlib.pyplot as plt

infile = open("ising_data.txt", 'r')
y = []

for line in infile:
    cols = line.split(' ')
    y.append(float(cols[0]))

x = np.linspace(0, len(y)-1, len(y))        # MC-cycles

plt.plot(x, y)

plt.grid('on')
plt.show()
