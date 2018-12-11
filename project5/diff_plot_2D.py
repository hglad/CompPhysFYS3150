import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.animation as animation

u = np.loadtxt('2D_forward_dx=0.010.txt')

print u.shape   # u[time, position]

n = len(u[0,:])
x = np.linspace(0, 1, n)

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
#plt.axis('equal')

#legend
cbar = plt.colorbar()
cbar.ax.get_yaxis().labelpad = 15
cbar.ax.set_ylabel('$T(x, y)$', rotation=270)
plt.xlabel('$x$'); plt.ylabel('$y$')
plt.title('Temperature distribution over time in a %s x %s grid' % (n,n))

i = 0
def updatefig(*args):
    global i
    if (i < t_steps):
        i += 1
    else:
        i=0         # reset animation
    im.set_array(mat[i-1])
    return im,

save = True
ani = animation.FuncAnimation(fig, updatefig, interval=1, blit=True, frames=int(t_steps/5))

if save == True:
    # Set up formatting for the movie files
    Writer = animation.writers['ffmpeg']      # requires ffmpeg to be installed
    writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)
    ani.save('2D_anim.mp4', writer=writer)
else:
    plt.show()
