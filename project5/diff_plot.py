import numpy as np
import matplotlib.pyplot as plt

u = np.loadtxt('crank_dx=0.010.txt')
u = u[:, 1:]            # Remove first column as it only contains indices
print u.shape   # u[time, position]
"""
#print u
n = 20000
t = np.linspace(0, n, n)
m = len(u[0,:])
for i in range(m):
#    print len(t), len(u[:,i])
    plt.plot(t, u[:,i])

#u2 = np.loadtxt('diffusion_dx=0.01.txt', usecols=1)

#plt.plot(t, u)
#plt.plot(t, u2)
plt.axis([0, n, 0, 1])
plt.show()
"""

n = len(u[0,:])
x = np.linspace(0, 1, n)
print x
# Put mmatplotlib in interactive mode for animation
plt.ion()

# Setup the figure before starting animation
fig = plt.figure() # Create window
ax = fig.add_subplot(111) # Add axes
line, = ax.plot( x, u[0,:], label='u(x)' ) # Fetch the line object
ax.set_ylim([0,1])

for i in range(n):
#    print len(t), len(u[:,i])
#    ax.plot(x, u[i,:])
    line.set_ydata( u[i, :] ) # Update the y values
    plt.draw() # Update the plot
    plt.pause(0.01)

# Turn off interactive mode
plt.ioff()

# Add show so that windows do not automatically close
plt.show()
