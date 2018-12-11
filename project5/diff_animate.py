import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

u = np.loadtxt('1D_crank_dx=0.010.txt')
exact_u = np.loadtxt('1D_analytical_dx=0.010.txt')
u = u[:, 1:]            # Remove first column as it only contains indices
exact_u = exact_u[:, 1:]
print u.shape   # u[time, position]

def time_compare(u, exact_u):
    nx = len(u[0,:])
    nx_exact = len(exact_u[0,:])
    t = len(u[:,0])
    x = np.linspace(0, 1, nx)
    x_exact = np.linspace(0,1,nx_exact)
    """
    for i in range(0, t, 10):
        plt.plot(x, u[i,:])
    plt.show()
    """
    #plt.plot(x, u[0,:])
    plt.plot(x, u[-1,:])
    #plt.plot(x, u[0,:])
    plt.plot(x_exact, exact_u[-1,:])
    plt.grid('on')
    plt.show()

def update_line(*args):
    global i
    if i < t:
        i += 1
    else:
        i=0
    lines[0].set_ydata(u[i-1, :])
    lines[1].set_ydata(exact_u[i-1, :])
    ax.set_title("Temperature distribution over time (dx = %s)\nTimestep: %s/%s" % (1./(nx-2), i-1, t))
#    ax.set_title("t = %1.3f / %1.1f" % (i*dt, 1))
    return lines

nx = len(u[0,:])
t = len(u[:,0])
x = np.linspace(0, 1, nx)
dt = 1./t

# Setup the figure before starting animation
fig, ax = plt.subplots(figsize=(8,8))
ax.set(xlim=(0,1), ylim=(0,1))

line, = ax.plot(x, u[0,:], label='Numerical')
line2, = ax.plot(x, exact_u[0,:], label='Analytical')
lines = [line, line2]      # FuncAnimation requires list of lines to update

lines[0].set_ydata(u[0, :])
lines[1].set_ydata(exact_u[0, :])

i = 0   # Counter for updating plot
ani = animation.FuncAnimation(fig, update_line, interval=10, blit=False, frames=t)

# Figure formatting
ax.set_ylim([0,1])
ax.grid('on')
ax.set_xlabel('$x$')
ax.set_ylabel('$T(x)$')
ax.legend(loc='upper left')

save = False
if save == True:
    # Set up formatting for the movie files
    Writer = animation.writers['ffmpeg']      # requires ffmpeg to be installed
    writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=3000)
    ani.save('1D_anim.mp4', writer=writer)
else:
    plt.show()


#time_compare(u, exact_u)