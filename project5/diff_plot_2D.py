import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.animation as animation

u = np.loadtxt('2D_forward_dx=0.010.txt')

#print u
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
    """
    cmap = mpl.colors.ListedColormap(['blue','black','red'])
    bounds=[-6, 2.2, 6]
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

    # tell imshow about color map so that only set colors are used
    img = pyplot.imshow(mat[t],interpolation='nearest',
                        cmap = cmap,norm=norm)

    # make a color bar
    pyplot.colorbar(img,cmap=cmap,
                    norm=norm,boundaries=bounds,ticks=[0,0.5,1])

    pyplot.show()
    """
    #plt.imshow(mat[t], cmap= cm.coolwarm)
    #plt.show()

fig = plt.figure(figsize=(8,8))
im = plt.imshow(mat[0], cmap=cm.coolwarm, animated=True)
#plt.axis('equal')
plt.colorbar()

i = 0
def updatefig(*args):
    global i
    if (i < t_steps/2):
        i += 1
    else:
        i=0         # reset animation
    im.set_array(mat[i])
    return im,

ani = animation.FuncAnimation(fig, updatefig, interval=1, blit=True)
plt.show()
