import numpy as np
import matplotlib.pyplot as plt

"""
Running this file produces relevant results for this project.
"""
# Shows analytical results for the L = 2 case.
def print_analytical(T=1):
    E_ = -np.sinh(8)*8 / (3 + np.cosh(8))
    Cv_ = 64/(T*T) * (np.cosh(8)*(3 + np.cosh(8)) - (np.sinh(8)**2))/(3 + np.cosh(8))**2
    M2_ = (8*np.exp(8) + 8)/(3 + np.cosh(8))
    absM_ = (2*np.exp(8) + 4)/(3 + np.cosh(8))
    chi_ = M2_/T
    print "Analytical expectation values at T = %1.2f:" % T
    print "Energy: %1.5f\nHeat capacity: %1.5f\nMag. squared: %1.5f\nAbsolute mag.: %1.5f\nSusceptibility: %1.5f\n" % (E_, Cv_, M2_, absM_, chi_)

# Used for plotting a single file, used for testing
def plot_single_arrays(file):
    E, M = np.loadtxt(file, usecols=(0,1), unpack=True)
    numMC = len(E)
    x = np.linspace(0, numMC-1, numMC)
    plt.plot(x, E)
    plt.grid('on')
    plt.show()

# Plot single energies and absolute magnetisation as function of MC-cycles
def plot_arrays_c(files_T1, files_T2, L):
    for file1, file2 in zip(files_T1, files_T2):
        E1, M1 = np.loadtxt(file1, usecols=(0,1), unpack=True, dtype='float')

        E2, M2 = np.loadtxt(file2, usecols=(0,1), unpack=True, dtype='float')
        numMC = len(E1)

        x = np.linspace(0, numMC-1, numMC)
        M1 = abs(M1)
        M2 = abs(M2)
        plt.plot(x, E1)
        plt.plot(x, E2)
        plt.title('Total energy per MC-cycle, %1.0f x %1.0f lattice' % (L,L))
        plt.legend(["T = 1.0", "T = 2.4"])
        plt.grid('on'); plt.ylabel('E'); plt.xlabel('MC-cycle')
        plt.figure()

        plt.plot(x, M1)
        plt.plot(x, M2)
        plt.title('Total magnetisation per MC-cycle, %1.0f x %1.0f lattice' % (L,L))
        plt.legend(["T = 1.0", "T = 2.4"])
        plt.grid('on'); plt.ylabel('M'); plt.xlabel('MC-cycle')

        plt.show()

# Expectation values for different MC cycles and temperatures
def plot_expectation_c(files_T1, files_T2, L):
    n = len(files_T1)
    E = np.zeros((n,2)); absM = np.zeros((n,2))
    M2 = np.zeros((n,2)); C_V = np.zeros((n,2))
    chi = np.zeros((n,2)); counts = np.zeros((n,2))
    MC = np.zeros((n,2)); T = np.zeros((n,2))

    i = 0
    for file1, file2 in zip(files_T1, files_T2):
        E[i,0], absM[i,0], M2[i,0], C_V[i,0], chi[i,0], counts[i,0], T[i,0], MC[i,0] = np.loadtxt(file1, usecols=(0,1,2,3,4,5,6,7), unpack=True)
        E[i,1], absM[i,1], M2[i,1], C_V[i,1], chi[i,1], counts[i,1], T[i,0], MC[i,1] = np.loadtxt(file2, usecols=(0,1,2,3,4,5,6,7), unpack=True)
        i += 1

    plt.title('Number of accepted states, %1.0f x %1.0f lattice' % (L,L), size=15)
    plt.semilogy(MC, counts, '-o')
    plt.xlabel('Num. of MC-cycles', size=15); plt.ylabel('Num. of accepted states', size=15)
    plt.grid('on')
    plt.legend(["T = 1.0", "T = 2.4"], prop = {'size':15})
    plt.figure()

    plt.title('Mean energy per spin, %1.0f x %1.0f lattice' % (L,L), size=15)
    plt.plot(MC, E, '-o')
    plt.xlabel('Num. of MC-cycles', size=15); plt.ylabel('$\\langle E \\rangle$', size=15)
    plt.grid('on')
    plt.legend(["T = 1.0", "T = 2.4"], prop = {'size':15})
    plt.figure()

    plt.title('Mean magnetisation per spin, %1.0f x %1.0f lattice' % (L,L), size=15)
    plt.plot(MC, absM, '-o')
    plt.xlabel('Num. of MC-cycles', size=15); plt.ylabel('$\\langle |M| \\rangle$', size=15)
    plt.grid('on')

    plt.legend(["T = 1.0", "T = 2.4"], prop = {'size':15})

    plt.show()

# Histograms for the probability distributions
def plot_hist_d():
    n = 400
    file1 = 'results/arrays_L=20T=1.00MC=100000.txt'
    file2 = 'results/arrays_L=20T=2.40MC=100000.txt'
    E1, M1 = np.loadtxt(file1, usecols=(0,1), unpack=True, dtype='float')
    E2, M2 = np.loadtxt(file2, usecols=(0,1), unpack=True, dtype='float')

    weights1 = np.ones_like(E1)/float(len(E1))
    weights2 = np.ones_like(E2)/float(len(E2))

    plt.figure()
    plt.title('Probabilities of energy states, 20 x 20 lattice at T = 1', size=15)
    plt.hist(E1/n, weights=weights1, bins=100)
    plt.xlabel('E', size=15); plt.ylabel('Probability', size=15)

    plt.figure()
    plt.title('Probabilities of energy states, 20 x 20 lattice at T = 2.4', size=15)
    plt.hist(E2/n, weights=weights2, bins=190)
    plt.xlabel('E', size=15); plt.ylabel('Probability', size=15)
    plt.show()

    means1 = 'results/means_L=20T=1.00-1.00MC=1000000.txt'
    means2 = 'results/means_L=20T=2.40-2.40MC=1000000.txt'

    E_T1, absM_T1, M2_T1, C_V_T1, chi_T1, counts_T1, T_T1, MC_T1 = np.loadtxt(means1, usecols=(0,1,2,3,4,5,6,7), unpack=True)
    E_T2, absM_T2, M2_T2, C_V_T2, chi_T2, counts_T2, T_T2, MC_T2 = np.loadtxt(means2, usecols=(0,1,2,3,4,5,6,7), unpack=True)

    # Variance of energy: heat capacity x (temperature)**2
    var_E1 = C_V_T1*(T_T1)**2
    var_E2 = C_V_T2*(T_T2)**2

    print "Energy variance at T = %1.2f: %f\nEnergy variance at T = %1.2f: %f\n" % (T_T1, var_E1, T_T2, var_E2)

# Plot expectation values as function of temperature for different lattice.
# Also perform linear regression to approximate crit. temperature as L --> inf.
def plot_temps_e():
    T_range = np.arange(2.20, 2.40+0.01, 0.01)

    # numL: number of different lattice dimensions
    # numT:    number of mean values per temperature
    L_strings = ['40', '60', '80', '100', '120']       # values of L to plot for
    L = np.array([1./40, 1./60, 1./80, 1./100, 1./120])
    T = np.zeros(len(L_strings))
    legends = []

    for i in range(len(L_strings)):
        legends.append('L = %s' % L_strings[i])

    MC_string = '1000000'
    numL = len(L_strings)
    numT = len(T_range)

    E = np.zeros((numT,numL)); absM = np.zeros((numT,numL))
    M2 = np.zeros((numT,numL)); C_V = np.zeros((numT,numL))
    chi = np.zeros((numT,numL)); counts = np.zeros((numT,numL))
    MC = np.zeros((numT,numL))

    T_start = '2.20'; T_final = '2.40'
    for i in range(len(L_strings)):
        file = ('results/means_L=%sT=%s-%sMC=%s.txt' % (L_strings[i], T_start, T_final, MC_string) )
        E[:,i], absM[:,i], M2[:,i], C_V[:,i], chi[:,i], counts[:,i], MC[:,i] = np.loadtxt(file, usecols=(0,1,2,3,4,5,6), unpack=True)
        T[i] = T_range[np.argmax(chi[:,i])] - L[i]

    plt.figure()
    plt.plot(T_range, chi)
    plt.title('Susceptibility of different lattices, %s MC-cycles' % MC_string, size=15)
    plt.legend(legends, prop={'size':15} )
    plt.xlabel('T', size=15); plt.ylabel('$\\langle \\chi \\rangle$', size=15)
    plt.grid('on')

    plt.figure()
    plt.plot(T_range, C_V)
    plt.title('Heat capacity of different lattices, %s MC-cycles' % MC_string, size=15)
    plt.legend(legends, prop={'size':15} )
    plt.xlabel('T', size=15); plt.ylabel('$\\langle C_V \\rangle$', size=15)
    plt.grid('on')

    plt.figure()
    plt.plot(T_range, absM)
    plt.title('Mean magnetisation of different lattices, %s MC-cycles' % MC_string, size=15)
    plt.legend(legends, prop={'size':15} )
    plt.xlabel('T', size=15); plt.ylabel('$\\langle |M| \\rangle$', size=15)
    plt.grid('on')

    plt.figure()
    plt.plot(T_range, E)
    plt.title('Mean energy of different lattices, %s MC-cycles' % MC_string, size=15)
    plt.legend(legends, prop={'size':15} )
    plt.xlabel('T', size=15); plt.ylabel('$\\langle E \\rangle$', size=15)
    plt.grid('on')

    plt.show()
    # Generate coefficients for 1st order polynomial to perform linear regression
    z = np.polyfit(L,T,1)
    zz = np.poly1d(z)

    print "Linear regression coefficients: %s, %s" % (zz[0], zz[1])

    plt.plot(L, T,'o')
    plt.plot(L, zz(L))
    plt.legend(['$T_C$', 'Fitted linear curve: ax + b = y\na = %1.5f, b = %1.5f ' % (zz[1], zz[0])], prop={'size':15})
    plt.title('Linear regression of critical temperature for different L',size=15)
    plt.grid('on')
    plt.xlabel('$L^{-1}$', size=15); plt.ylabel('$T_C$')
    plt.show()

#----------------------- c ----------------------
# L = 20
# files with different MC-cycles for means
MC1 = ['5000', '10000', '15000', '20000', '25000', '30000', '35000', '40000', '45000', '50000']
means1 = []
means2 = []
arrays1 = []
arrays2 = []
for MC in MC1:
    means1.append('results/means_L=20T=1.00-1.00MC=%s.txt' % MC)
    means2.append('results/means_L=20T=2.40-2.40MC=%s.txt' % MC)
    #arrays_T1.append('results/arrays_L=20T=1.00MC=%s.txt' % MC)

# files with different MC-cycles for energy and magnetisation arrays
MC2 = ['5000', '10000', '15000', '20000']

for MC in MC2:
    arrays1.append('results/arrays_L=20T=1.00MC=%s.txt' % MC)
    arrays2.append('results/arrays_L=20T=2.40MC=%s.txt' % MC)

print_analytical()
plot_expectation_c(means1, means2, L=20)
plot_arrays_c(arrays1, arrays2, L=20)

#----------------------- d ----------------------
plot_hist_d()

#----------------------- e ----------------------
plot_temps_e()


#plot_single_arrays('results/arrays_L=20T=1.00MC=10000.txt')
#plot_single_means('means/arrays_L=20T=1.00MC=10000.txt')
