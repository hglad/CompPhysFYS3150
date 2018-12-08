import numpy as np
import matplotlib.pyplot as plt

u = np.loadtxt('1D_forward_dx=0.100.txt')
exact_u = np.loadtxt('1D_analytical_dx=0.100.txt')
u = u[:, 1:]            # Remove first column as it only contains indices
exact_u = exact_u[:, 1:]
print u.shape   # u[time, position]

def time_compare(u):
    nx = len(u[0,:])
    t = len(u[:,0])
    x = np.linspace(0, 1, nx)

    for i in range(0, t, 10):
        plt.plot(x, u[i,:])
    plt.show()
    """
    plt.plot(x, u[0,:])
    plt.plot(x, u[10, :])
    plt.plot(x, u[-1,:])
    plt.show()
    """

def animate(u, exact_u):
    nx = len(u[0,:])
    t = len(u[:,0])
    x = np.linspace(0, 1, nx)

    # Put matplotlib in interactive mode for animation
    plt.ion()

    # Setup the figure before starting animation
    fig = plt.figure()
    ax = fig.add_subplot(111)
    line, = ax.plot( x, u[0,:], label='Numerical' ) # Fetch the line object
    line2, = ax.plot( x, exact_u[0,:], label='Exact' )

    ax.set_ylim([0,1])
    plt.title('Timestep: 0')
    ax.grid('on')
    ax.set_xlabel('x')
    ax.set_ylabel('u(x)')
    ax.legend(loc='best')

    print nx
    for i in range(t):
    #    print len(t), len(u[:,i])
    #    ax.plot(x, u[i,:])
        line.set_ydata( u[i, :] ) # Update the y values
        line2.set_ydata( exact_u[i, :])
        plt.title('Timestep: %s' % i)
        plt.draw() # Update the plot
        plt.pause(0.1)

    # Turn off interactive mode
    plt.ioff()
    # Add show so that windows do not automatically close
    plt.show()

time_compare(u)
animate(u, exact_u)

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
