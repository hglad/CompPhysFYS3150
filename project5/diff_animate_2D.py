import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.animation as animation

fe_u = np.loadtxt('2D_forward_dx=0.010.txt')
be_u = np.loadtxt('2D_backward_dx=0.010.txt')
print fe_u.shape, be_u.shape   # u[time, position]
u = abs(fe_u - be_u)
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
im = plt.imshow(mat[0], cmap=cm.terrain, animated=True)
# plt.clim(0,0.15)
# Figure formatting
cbar = plt.colorbar()
cbar.ax.get_yaxis().labelpad = 15
cbar.ax.set_ylabel('$\\Delta T$', rotation=270, size=15)
plt.xlabel('$x$', size=15); plt.ylabel('$y$', size=15)
#plt.title('Temperature distribution over\ntime in a %s x %s grid' % (n,n), size=15)

i = 0
def updatefig(*args):
    global i
    if (i < t_steps/3):
        i += 1
    else:
        i=0         # reset animation
    im.set_array(mat[i-1])
    plt.title(" (dx = %s)\nTimestep: %s/%s" % (1./(n), i-1, t))
    return im,

save = True
ani = animation.FuncAnimation(fig, updatefig, interval=1, blit=False, frames=int(t_steps/3))

if save == True:
    # Set up formatting for the movie files
    Writer = animation.writers['ffmpeg']      # requires ffmpeg to be installed
    writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)
    ani.save('2D_anim_deltaT_dx=%s.mp4' % dx, writer=writer)
else:
    plt.show()
