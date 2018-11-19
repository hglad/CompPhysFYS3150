import numpy as np
import matplotlib.pyplot as plt
import numpy.linalg
import glob

def plot_single_arrays(file):
    E, M = np.loadtxt(file, usecols=(0,1), unpack=True)
    numMC = len(E)
    x = np.linspace(0, numMC-1, numMC)
    plt.plot(x, E)
    plt.grid('on')
    plt.show()

def plot_arrays_c(files_T1, files_T2, L):
    print files_T1

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

def plot_expectation_c(files_T1, files_T2, L):
    n = len(files_T1)
    E = np.zeros((n,2)); absM = np.zeros((n,2))
    M2 = np.zeros((n,2)); C_V = np.zeros((n,2))
    chi = np.zeros((n,2)); counts = np.zeros((n,2))
    MC = np.zeros((n,2)); T = np.zeros((n,2))
    #plot_array(filename)
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
    # plt.legend(legends, prop={'size':15} )
    plt.legend(["T = 1.0", "T = 2.4"], prop = {'size':15})

    plt.show()

def plot_hist_d():
    n = 400
    file1 = 'results/arrays_L=20T=1.00MC=100000.txt'
    file2 = 'results/arrays_L=20T=2.40MC=100000.txt'
    E1, M1 = np.loadtxt(file1, usecols=(0,1), unpack=True, dtype='float')
    E2, M2 = np.loadtxt(file2, usecols=(0,1), unpack=True, dtype='float')

    weights1 = np.ones_like(E1)/float(len(E1))
    weights2 = np.ones_like(E2)/float(len(E2))
#    plt.hist(myarray, weights=weights)
    plt.figure()
    plt.title('States per energy, 20 x 20 lattice at T = 1')
#    plt.hist(E1, bins=100, density=True, stacked=True)
    plt.hist(E1/n, weights=weights1, bins=100)
    plt.xlabel('E'); plt.ylabel('Number of states')

    plt.figure()
    plt.title('States per energy, 20 x 20 lattice at T = 2.4')
#    plt.hist(E2, bins=190, density=True, stacked=True)
    plt.hist(E2/n, weights=weights2, bins=190)
    plt.xlabel('E'); plt.ylabel('Number of states')
    plt.show()

    means1 = 'results/means_L=20T=1.00-1.00MC=1000000.txt'
    means2 = 'results/means_L=20T=2.40-2.40MC=1000000.txt'

    E_T1, absM_T1, M2_T1, C_V_T1, chi_T1, counts_T1, T_T1, MC_T1 = np.loadtxt(means1, usecols=(0,1,2,3,4,5,6,7), unpack=True)
    E_T2, absM_T2, M2_T2, C_V_T2, chi_T2, counts_T2, T_T2, MC_T2 = np.loadtxt(means2, usecols=(0,1,2,3,4,5,6,7), unpack=True)

    # Variance of energy: heat capacity x (temperature)**2
    var_E1 = C_V_T1*(T_T1)**2
    var_E2 = C_V_T2*(T_T2)**2

    print "Energy variance at T = %1.2f: %f\nEnergy variance at T = %1.2f: %f" % (T_T1, var_E1, T_T2, var_E2)

# Plot expectation values as function of temperature for different lattices
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

    #T = np.array([T_range[np.argmax(chi[4])] - L[0], T_range[np.argmax(chi[3])]- L[1], T_range[np.argmax(chi[2])]- L[2], T_range[np.argmax(chi[1])]- L[3], T_range[np.argmax(chi[0])]- L[4]])

    #L = np.array([1.0/120, 1.0/100, 1.0/80, 1.0/60, 1.0/40])
    #T = np.array([T_range[np.argmax(chi[4])] - L[4], T_range[np.argmax(chi[3])] - L[3], T_range[np.argmax(chi[2])]- L[2], T_range[np.argmax(chi[1])]- L[1], T_range[np.argmax(chi[0])]- L[0]])

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

    z = np.polyfit(L,T,1)
    zz = np.poly1d(z)

    print(zz[0], zz[1])
#    lin_L = [0, L[0]]
    plt.plot(L, T,'o')
    plt.plot(L, zz(L))
#    plt.plot(lin_L, zz(lin_L) )
    plt.legend(['$T_C$', 'Fitted linear curve'], prop={'size':15})
    plt.title('Linear regression of critical temperature for different L',size=15)
    plt.grid('on')
    plt.xlabel('$L^{-1}$', size=15); plt.ylabel('$T_C$')
    plt.show()

#----------------------- c ----------------------
# L = 20
# files with different MC-cycles with T = 1
T1 = ['5000', '10000', '15000', '20000', '25000', '30000', '35000', '40000', '45000', '50000']
means_T1 = []
arrays_T1 = []
for MC in T1:
    means_T1.append('results/means_L=20T=1.00-1.00MC=%s.txt' % MC)
    arrays_T1.append('results/arrays_L=20T=1.00MC=%s.txt' % MC)

# files with different MC-cycles with T = 2.4
T2 = T1
means_T2 = []
arrays_T2 = []

for MC in T2:
    means_T2.append('results/means_L=20T=2.40-2.40MC=%s.txt' % MC)
    arrays_T2.append('results/arrays_L=20T=2.40MC=%s.txt' % MC)

#----------------------- e ----------------------       1 mill. cycles
#T_L20 = ['2.20', '2.21', '2.22', '2.23', '2.24', '2.25', '2.26', '2.27', '2.28', '2.29', '2.30', '2.31', '2.32', '2.33', '2.34', '2.35', '2.36', '2.37', '2.38', '2.39', '2.40']
"""
T_L40 = ['2.20', '2.21', '2.22', '2.23', '2.24', '2.25', '2.26', '2.27', '2.28', '2.29', '2.30']
T_L20 = T_L40

T_L20_range = np.arange(2.20, 2.30+0.01, 0.01)
T_L40_range = np.arange(2.20, 2.30+0.01, 0.01)

means_T_L20 = []
means_T_L40 = []
for i in range(len(T_L40)):
    means_T_L40.append('results/means_L=40T=%sMC=100000.txt' % T_L40[i])

for i in range(len(T_L20)):
    means_T_L20.append('results/means_L=20T=%sMC=100000.txt' % T_L20[i])
"""
#plot_expectation_c(means_T1, means_T2, L=20)
#plot_arrays_c(arrays_T1, arrays_T2, L=20)

plot_hist_d()
#plot_temps_e()


#plot_single_arrays('results/arrays_L=20T=1.00MC=10000.txt')
#plot_single_means('means/arrays_L=20T=1.00MC=10000.txt')


"""
T_L20 = np.arange(2.20, 2.40+0.01, 0.01)
for i in range(len(T_L20)):
    #T_L20[i] = str(T_L20[i])
    T_L20[i] = "%.2f" % T_L20[i]
    means_T_L20.append('results/means_L=20T=%sMC=1000000.txt' % T_L20[i])
    #temp = T_L20[i].substr(0,4);
"""
