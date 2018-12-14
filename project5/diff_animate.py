import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

u1 = np.loadtxt('1D_forward_dx=0.010_alpha=0.5.txt')
u2 = np.loadtxt('1D_backward_dx=0.010_alpha=0.5.txt')
u3 = np.loadtxt('1D_crank_dx=0.010_alpha=0.5.txt')
exact_u = np.loadtxt('1D_analytical_dx=0.010.txt')
exact_u = exact_u[:, 1:]
#print u.shape   # u[time, position]

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
    plt.plot(x, u[-1,:])
    plt.plot(x_exact, exact_u[-1,:])
    plt.grid('on')
    plt.show()

def update_line(*args):
    global i
    if i < t:
        i += 1
    else:
        i=0
    lines[0].set_ydata(u1[i-1, :])
    lines[1].set_ydata(u2[i-1, :])
    lines[2].set_ydata(u3[i-1, :])
    lines[3].set_ydata(exact_u[i-1, :])
    t_ = float(i)/float(t);
    ax.set_title("Temperature distribution at t = %1.2f\ndx = %s" % (t_, 1./(nx-2)), size=15)
#    ax.set_title("t = %1.3f / %1.1f" % (i*dt, 1))
    return lines

nx = len(u1[0,:])
t = len(u1[:,0])
x = np.linspace(0, 1, nx)
dt = 1./t

# Setup the figure before starting animation
fig, ax = plt.subplots(figsize=(8,8))
ax.set(xlim=(0,1), ylim=(0,1))

line, = ax.plot(x, u1[0,:], label='Forward Euler')
line2, = ax.plot(x, u2[0,:], label='Backward Euler')
line3, = ax.plot(x, u3[0,:], label='Crank-Nicolson')
line4, = ax.plot(x, exact_u[0,:], label='Analytical')
lines = [line, line2, line3, line4]      # FuncAnimation requires list of lines to update

lines[0].set_ydata(u1[0, :])
lines[1].set_ydata(u2[0, :])
lines[2].set_ydata(u3[0, :])
lines[3].set_ydata(exact_u[0, :])

i = 0   # Counter for updating plot
ani = animation.FuncAnimation(fig, update_line, interval=10, blit=False, frames=int(t/2))

# Figure formatting
ax.set_ylim([0,1])
ax.grid('on')
ax.set_xlabel('$x$', size=15)
ax.set_ylabel('$T(x)$', size=15)
ax.legend(loc='upper left', prop = {'size':15})

save = True
if save == True:
    # Set up formatting for the movie files
    Writer = animation.writers['ffmpeg']      # requires ffmpeg to be installed
    writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=3000)
    ani.save('1D_anim_dx=%1.2f.mp4' % (1./(nx-2)), writer=writer)
else:
    plt.show()


#time_compare(u, exact_u)
