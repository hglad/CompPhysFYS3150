import numpy as np
import matplotlib.pyplot as plt

# Open file containing data for task b), read lines and store data in arrays
infile = open("project3a.txt", 'r')
x = []
y = []
for line in infile:
    cols = line.split(' ')
    x.append(float(cols[0]))
    y.append(float(cols[1]))

plt.plot(x, y)
plt.grid('on')
plt.show()
