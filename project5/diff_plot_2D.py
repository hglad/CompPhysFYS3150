import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.animation as animation

an_u = np.loadtxt('2D_analytical_dx=0.010.txt')

#print be_u.shape, an_u.shape
t = 0.05
nt_an = int(len(an_u[:,0])/len(an_u[0,:]))
print nt_an
#nt_an = len(an_u[:,0])
t_an_ = int(nt_an*t)     # index
nx_exact = len(an_u[0,:])
x_exact = np.linspace(0,1,nx_exact)
mat_an = np.zeros((nt_an, nx_exact, nx_exact))

for t in range(nt_an):
    start = nx_exact*t
    end = start+nx_exact
    mat_an[t] = an_u[start:end]

#for s in ['0.25', '0.50', '1.00', '2.00', '4.00']:
for s in ['8.00']:
    be_u = np.loadtxt("2D_backward_dx=0.010_alpha=" + s + ".txt")
    print be_u.shape
#    u = abs(be_u - an_u)
    #u = an_u
    #print u.shape   # u[time, position]
    #t_ = int(nt*t)     # index for chosen time to plot for
    nx = len(be_u[0,:])
    nt = int(len(be_u[:,0])/len(be_u[0,:]))
    x = np.linspace(0, 1, nx)
    dx = 1./nx

    mat = np.zeros((nt, nx, nx))

    # Create matrices for individual time steps
    for t in range(nt):
        start = nx*t
        end = start+nx
        mat[t] = be_u[start:end]


    # Plot specific points in time
    #index = 0
    #t = float(index)/float(nt)
    t = 0.05
    index = int(t*nt)
    u = abs(mat[index] - mat_an[t_an_])
    #print index, t
    fig = plt.figure(figsize=(6,6))
    im = plt.imshow(u, cmap=cm.terrain)
    cbar = plt.colorbar()
    plt.clim(0,0.01)
    #plt.clim(0,1)
    cbar.ax.get_yaxis().labelpad = 15
    cbar.ax.set_ylabel('$\\Delta T$', rotation=270, size=15)
    plt.xlabel('$x$', size=15); plt.ylabel('$y$', size=15)
    plt.title('Difference between solutions from Backward Euler\nand analytical solution at t = %1.2f\nin a %s x %s grid (alpha = %s)' % (t,nx,nx,s), size=15)
    #plt.title('Temperature distribution at t = %1.2f\nin a %s x %s grid (analytic)' % (t,nx,nx), size=15)
    #plt.savefig('figures/fe_error_2D_dx=%s_t=%1.2f.eps' % (dx,t))
    plt.savefig('figures/be_error_2D_dx=%s_t=%1.2f_alpha=%s.pdf' % (dx,t,s))
plt.show()
