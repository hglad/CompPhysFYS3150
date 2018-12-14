import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.animation as animation

fe_u = np.loadtxt('2D_forward_dx=0.100_alpha=0.25.txt')
be_u = np.loadtxt('2D_backward_dx=0.010_alpha=1.00.txt')
#an_u = np.loadtxt('2D_analytical_dx=0.111.txt')
an_u = np.loadtxt('2D_analytical_dx=0.010.txt')
"""
fe_u = np.loadtxt('2D_forward_dx=0.010_alpha=0.25.txt')
be_u = np.loadtxt('2D_backward_dx=0.010_alpha=0.25.txt')
an_u = np.loadtxt('2D_analytical_dx=0.010.txt')
"""
print fe_u.shape, be_u.shape, an_u.shape   # u[time, position]
u = abs(an_u - be_u)
#u = an_u
#u = be_u
#u = be_u

n = len(u[0,:])
x = np.linspace(0, 1, n)
dx = 1./n

t_steps = int(len(u[:,0])/len(u[0,:]))  # get number of time steps from matrix
print t_steps
mat = np.zeros((t_steps, n, n))
# Create matrices for individual time steps
for t in range(t_steps):
    start = n*t
    end = start+n
    mat[t] = u[start:end]

fig = plt.figure(figsize=(8,8))
im = plt.imshow(mat[0], cmap=cm.coolwarm, animated=True)
plt.clim(0,1)
# Figure formatting
cbar = plt.colorbar()
cbar.ax.get_yaxis().labelpad = 15
#cbar.ax.set_ylabel('$\\Delta T$', rotation=270, size=15)
cbar.ax.set_ylabel('$\\Delta T$', rotation=270, size=15)
plt.xlabel('$x$', size=15); plt.ylabel('$y$', size=15)
#plt.title('Temperature distribution at t = %s\nin a %s x %s grid (Backward Euler)' % (0,n,n), size=15)
plt.title('Difference between solutions from Backward Euler\nand analytical solution at t = %1.2f\nin a %s x %s grid' % (0,n,n), size=15)

i = 0
def updatefig(*args):
    global i
    if (i < t_steps/4):
        i += 1
    else:
        i=0         # reset animation
    im.set_array(mat[i-1])
    t = float(i)/float(t_steps);
#    plt.title('Temperature distribution at t = %1.2f\nin a %s x %s grid (Backward Euler)' % (t,n,n), size=15)
    plt.title('Difference between solutions from Backward Euler\nand analytical solution at t = %1.2f\nin a %s x %s grid' % (t,n,n), size=15)
    return im,

save = False
ani = animation.FuncAnimation(fig, updatefig, interval=100, blit=False, frames=int(t_steps/4))

if save == True:
    # Set up formatting for the movie files
    Writer = animation.writers['ffmpeg']      # requires ffmpeg to be installed
    writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)
    ani.save('2D_anim_be_error_dx=%s.mp4' % dx, writer=writer)
else:
    plt.show()
