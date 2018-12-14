import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

u1 = np.loadtxt('1D_forward_dx=0.010_alpha=0.5.txt')
u2 = np.loadtxt('1D_backward_dx=0.010_alpha=0.5.txt')
u3 = np.loadtxt('1D_crank_dx=0.010_alpha=0.5.txt')
"""
u1 = np.loadtxt('1D_forward_dx=0.100_alpha=0.5.txt')
u2 = np.loadtxt('1D_backward_dx=0.100_alpha=0.5.txt')
u3 = np.loadtxt('1D_crank_dx=0.100_alpha=0.5.txt')
"""

# Load arrays for calculating max error per alpha value
e_fe = []
e_be = []
e_cn = []
alpha_list = [0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
for alpha in alpha_list:
    e_fe.append(np.loadtxt('1D_forward_dx=0.010_alpha=%s.txt' % alpha))
    e_be.append(np.loadtxt('1D_backward_dx=0.010_alpha=%s.txt' % alpha))
    e_cn.append(np.loadtxt('1D_crank_dx=0.010_alpha=%s.txt' % alpha))


u_exact = np.loadtxt('1D_analytical_dx=0.010.txt')
# u1 = u1[:, 1:]            # Remove first column as it only contains indices
# u2 = u2[:, 1:]
# u3 = u3[:, 1:]
u_exact = u_exact[:, 1:]

def time_compare(u_arrays, u_exact, t):
    # Find indices for analytical array
    nx_exact = len(u_exact[0,:])
    x_exact = np.linspace(0,1,nx_exact)

    t_max_exact = len(u_exact[:,0])
    t_exact_ = int(t_max_exact*t)     # index

    plt.figure()
    # Loop through given solutions, finding correct time step for each and plotting
    for u_ in u_arrays:
        nx_ = len(u_[0,:])
        x_ = np.linspace(0, 1, nx_)
        t_max = len(u_[:,0])
        t_ = int(t_max*t)     # index for chosen time to plot for
        dx = 1./(nx_-2)
        plt.plot(x_, u_[t_,:])

    plt.plot(x_exact, u_exact[t_exact_,:])

    plt.grid('on')
    plt.xlabel('$x$', size=15)
    plt.ylabel('$Absolute error$', size=15)
    plt.legend(['Forward Euler', 'Backward Euler', 'Crank-Nicolson', 'Analytical'], loc='upper left', prop = {'size':15})

    plt.title("Temperature distribution in 1D rod\ndx = %s" % dx, size=15)
    plt.savefig('figures/T(x)_compare_dx=%s_t=%s_199.eps' % (dx,t_))
    plt.savefig('figures/T(x)_compare_dx=%s_t=%s_199.pdf' % (dx,t_))

#    print "Error: %s" % (np.std(u[t_,:] - u_exact[t_exact_,:]))
    print "dx = %s  t = %s/%s" % (dx, t_, t_max)

def compare_error(u_arrays, u_exact, t):
    # Find indices for analytical array
    nx_exact = len(u_exact[0,:])
    x_exact = np.linspace(0,1,nx_exact)
    print nx_exact

    t_max_exact = len(u_exact[:,0])
    t_exact_ = int(t_max_exact*t)     # index

    # Find correct indices for analytical array
    u1 = u_arrays[0]
    print u1.shape, u_exact.shape
    nx_ = len(u1[0,:])
    x_ = np.linspace(0, 1, nx_)
    t_max = len(u1[:,0])
    t_ = int(t_max*t)     # index for chosen time to plot for
    dx = 1./(nx_-2)

    # Compare length of analytical vs. numerical arrays to find correct indices
    fac = ((nx_exact-2)/(nx_-2))
    indices = range(0, nx_exact, fac)
    print indices

    if (nx_exact-1 not in indices) == True:       # get end point
        indices.append(nx_exact-1)
    # set u_exact to same shape as u
    u_exact = np.take(u_exact[t_,:], indices)

    print u1.shape, u_exact.shape

    plt.figure()
    # Loop through given solutions and plot
    for u_ in u_arrays:
        u_error = abs(u_[t_,:] - u_exact)
        plt.semilogy(x_[1:-2], u_error[1:-2])
        # exclude end points since they are identical to analytical solution

    # Exact solution

#    plt.plot(x_exact, u_exact[t_exact_,:])

    plt.grid('on')
    plt.xlabel('$x$', size=15)
    plt.ylabel('Absolute error', size=15)
    plt.legend(['Forward Euler', 'Backward Euler', 'Crank-Nicolson'], loc='best', prop = {'size':15})

    plt.title("Absolute error for numerical solvers at t = %1.2f\ndx = %s" % (t,dx), size=15)
    plt.savefig('figures/error_dx=%s_t=%1.2f.eps' % (dx,t))
    plt.savefig('figures/error_dx=%s_t=%1.2f.pdf' % (dx,t))

#    print "Error: %s" % (np.std(u[t_,:] - u_exact[t_exact_,:]))
    print "dx = %s  t = %s/%s" % (dx, t_, t_max)

# Plot one data point per array
def max_error_plot(e_fe, e_be, e_cn, u_exact, alpha_list):
    # Find indices for analytical array
    nx_exact = len(u_exact[0,:])
    nt_exact = len(u_exact[:,0])
    x_exact = np.linspace(0,1,nx_exact)
    print nx_exact

    n_arrays = len(e_fe)
    print n_arrays

    # Find correct indices for analytical array
    u1 = e_fe[0]

    print u1.shape, u_exact.shape
    nx_ = len(u1[0,:])
    x_ = np.linspace(0, 1, nx_)

    dx = 1./(nx_-2)

    # Compare length of analytical vs. numerical arrays to find correct indices
    fac = ((nx_exact-2)/(nx_-2))
    print nx_exact, nx_
    indices = range(0, nx_exact, fac)

    if (nx_exact-1 not in indices) == True:       # get end point
        indices.append(nx_exact-1)

    # set u_exact to same shape as u
    for i in range(nt_exact):
        u_exact[i,:] = np.take(u_exact[i,:], indices)
    print u_exact

    error_list = []
    for i in range(n_arrays):
        print e_fe[i].shape
    #    error = np.max(abs( e_fe[i] - u_exact))
    #    error_list.append(error)

    plt.plot(alpha_list, error_list)
    plt.show()

t1 = 0.05         # smooth, curved
t2 = 0.3        # almost linear
u_arrays = [u1, u2, u3]
#time_compare(u_arrays, u_exact, t1)
#time_compare(u_arrays, u_exact, t2)
#compare_error(u_arrays, u_exact, t1)
compare_error(u_arrays, u_exact, t2)

# not used
#max_error_plot(e_fe, e_be, e_cn, u_exact, alpha_list)

plt.show()
