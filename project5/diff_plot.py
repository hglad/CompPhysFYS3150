import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

u = np.loadtxt('1D_crank_dx=0.010.txt')
u_exact = np.loadtxt('1D_analytical_dx=0.010.txt')
u = u[:, 1:]            # Remove first column as it only contains indices
u_exact = u_exact[:, 1:]
print u.shape   # u[time, position]

def time_compare(u, u_exact, t):
    nx = len(u[0,:])
    nx_exact = len(u_exact[0,:])

    x = np.linspace(0, 1, nx)
    x_exact = np.linspace(0,1,nx_exact)

    # Find indices for numerical and analytical arrays
    t_max = len(u[:,0])
    t_ = int(t_max*t)     # index
    t_max_exact = len(u_exact[:,0])
    t_exact_ = int(t_max_exact*t)     # index

    dx = 1./(nx-2)
    plt.figure()
    plt.plot(x, u[t_,:])
    plt.plot(x_exact, u_exact[t_exact_,:])
    plt.grid('on')
    plt.xlabel('$x$', size=15)
    plt.ylabel('$T(x)$', size=15)
    plt.legend(['Crank-Nicolson', 'Analytical'], loc='upper left', prop = {'size':15})
    #plt.title("Temperature distribution in rod (dx = %s)\nTimestep: %s/%s" % (dx, t, t_max))
    plt.title("Temperature distribution in 1D rod\ndx = %s" % dx, size=15)
    plt.savefig('figures/T(x)_CN_dx=%s_t=%s_199.eps' % (dx,t_))
    plt.savefig('figures/T(x)_CN_dx=%s_t=%s_199.pdf' % (dx,t_))
    print "dx = %s  t = %s/%s" % (dx, t_, t_max)


t1 = 0.05         # smooth, curved
t2 = 0.3        # almost linear
time_compare(u, u_exact, t1)
time_compare(u, u_exact, t2)

plt.show()
