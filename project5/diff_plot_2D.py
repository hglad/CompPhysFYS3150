import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.animation as animation

fe_u = np.loadtxt('2D_forward_dx=0.010.txt')
be_u = np.loadtxt('2D_backward_dx=0.010.txt')
#u = abs(fe_u - be_u)
u = be_u
print u.shape   # u[time, position]

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

# Plot specific points in time
fig = plt.figure(figsize=(8,8))
im = plt.imshow(mat[-1], cmap=cm.coolwarm)
cbar = plt.colorbar()
#plt.clim(0,0.15)
cbar.ax.get_yaxis().labelpad = 15
cbar.ax.set_ylabel('$T$', rotation=270, size=15)
plt.xlabel('$x$', size=15); plt.ylabel('$y$', size=15)
#plt.title('Difference between solutions from Forward Euler\nand Backward Euler schemes\nin a %s x %s grid' % (n,n), size=15)
plt.title('Initial temperature distribution\nin a %s x %s grid' % (n,n), size=15)
plt.savefig('figures/T(x,y)_2D_dx=%s_steady.eps' % (dx))
plt.savefig('figures/T(x,y)_2D_dx=%s_steady.pdf' % (dx))

plt.show()
