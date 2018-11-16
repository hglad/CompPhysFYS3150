import numpy as np
import matplotlib.pyplot as plt
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
        plt.grid('on'); plt.ylabel('E [$JL^2$]'); plt.xlabel('MC-cycle')
        plt.figure()

        plt.plot(x, M1)
        plt.plot(x, M2)
        plt.title('Total magnetisation per MC-cycle, %1.0f x %1.0f lattice' % (L,L))
        plt.legend(["T = 1.0", "T = 2.4"])
        plt.grid('on'); plt.ylabel('M [$s$]'); plt.xlabel('MC-cycle')

        plt.show()

def plot_expectation_c(files_T1, files_T2, L):
    n = len(files_T1)
    E = np.zeros((n,2)); absM = np.zeros((n,2))
    M2 = np.zeros((n,2)); C_V = np.zeros((n,2))
    chi = np.zeros((n,2)); counts = np.zeros((n,2))
    MC = np.zeros((n,2))
    #plot_array(filename)
    i = 0
    for file1, file2 in zip(files_T1, files_T2):
        E[i,0], absM[i,0], M2[i,0], C_V[i,0], chi[i,0], counts[i,0], MC[i,0] = np.loadtxt(file1, usecols=(0,1,2,3,4,5,6), unpack=True)
        E[i,1], absM[i,1], M2[i,1], C_V[i,1], chi[i,1], counts[i,1], MC[i,1] = np.loadtxt(file2, usecols=(0,1,2,3,4,5,6), unpack=True)
        i += 1

    plt.title('Number of accepted states, %1.0f x %1.0f lattice' % (L,L))
    plt.semilogy(MC, counts, '-o')
    plt.xlabel('Num. of MC-cycles'); plt.ylabel('Num. of accepted states')
    plt.grid('on')
    plt.legend(["T = 1.0 kT/J", "T = 2.4 kT/J"])
    plt.figure()

    plt.title('Mean energy per spin, %1.0f x %1.0f lattice' % (L,L))
    plt.plot(MC, E, '-o')
    plt.xlabel('Num. of MC-cycles'); plt.ylabel('Mean energy [$J$]')
    plt.grid('on')
    plt.legend(["T = 1.0", "T = 2.4"])
    plt.figure()

    plt.title('Mean magnetisation per spin, %1.0f x %1.0f lattice' % (L,L))
    plt.plot(MC, absM, '-o')
    plt.xlabel('Num. of MC-cycles'); plt.ylabel('Mean magnetisation [s]')
    plt.grid('on')
    plt.legend(["T = 1.0 kT/J", "T = 2.4 kT/J"])

    plt.show()

# Plot expectation values as function of temperature for different lattices
def plot_temps_e():
    #T_strings = ['2.20', '2.21', '2.22', '2.23', '2.24', '2.25', '2.26', '2.27', '2.28', '2.29', '2.30', '2.31', '2.32', '2.33', '2.34', '2.35', '2.36', '2.37', '2.38', '2.39', '2.40']
    #T_strings = ['2.20', '2.21', '2.22', '2.23', '2.24', '2.25', '2.26', '2.27', '2.28', '2.29', '2.30']

    T_range = np.arange(2.20, 2.40+0.01, 0.01)

    # numL: number of different lattice dimensions
    # numT:    number of mean values per temperature
    L_strings = ['20', '40', '60', '80']       # values of L to plot for
    numL = len(L_strings)
    numT = len(T_range)

    E = np.zeros((numT,numL)); absM = np.zeros((numT,numL))
    M2 = np.zeros((numT,numL)); C_V = np.zeros((numT,numL))
    chi = np.zeros((numT,numL)); counts = np.zeros((numT,numL))
    MC = np.zeros((numT,numL))

    T_start = '2.20'; T_final = '2.40'
    for i in range(len(L_strings)):
        file = ('results/means_L=%sT=%s-%sMC=100000.txt' % (L_strings[i], T_start, T_final) )
        E[:,i], absM[:,i], M2[:,i], C_V[:,i], chi[:,i], counts[:,i], MC[:,i] = np.loadtxt(file, usecols=(0,1,2,3,4,5,6), unpack=True)
#    T_L20 = np.arange(2.20, 2.40+0.01, 0.01)

#    plt.title('$\\chi$, %1.0f x %1.0f lattice' % (L,L))
    legends = ['L = 20', 'L = 40', 'L = 60', 'L = 80', 'L = 100']
    plt.figure()
    plt.plot(T_range, chi)
    plt.legend(legends)
    plt.xlabel('T [kT/J]'); plt.ylabel('$\\chi$')
    plt.grid('on')

    plt.figure()
    plt.plot(T_range, C_V)
    plt.legend(legends)
    plt.xlabel('T [kT/J]'); plt.ylabel('$C_V$')
    plt.grid('on')

    plt.figure()
    plt.plot(T_range, absM)
    plt.legend(legends)
    plt.xlabel('T [kT/J]'); plt.ylabel('$|M|$')
    plt.grid('on')

    plt.figure()
    plt.plot(T_range, E)
    plt.legend(legends)
    plt.xlabel('T [kT/J]'); plt.ylabel('$E$')
    plt.grid('on')

    plt.show()

#----------------------- c ----------------------
# L = 20
# files with different MC-cycles with T = 1
T1 = ['5000', '10000', '15000', '20000', '25000', '30000', '35000', '40000']
means_T1 = []
arrays_T1 = []
for MC in T1:
    means_T1.append('results/means_L=20T=1.00-1.00MC=%s.txt' % MC)
    arrays_T1.append('results/arrays_L=20T=1.00-100MC=%s.txt' % MC)

# files with different MC-cycles with T = 2.4
T2 = T1
means_T2 = []
arrays_T2 = []

for MC in T2:
    means_T2.append('results/means_L=20T=2.40-2.40MC=%s.txt' % MC)
    arrays_T2.append('results/arrays_L=20T=2.40-2.40MC=%s.txt' % MC)

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

plot_temps_e()
#

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
