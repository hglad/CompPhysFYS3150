import numpy as np
import matplotlib.pyplot as plt
import glob

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
# L = 20
# files with different MC-cycles with T = 1
T1 = ['1000', '2500', '5000', '6500', '8000', '10000']
means_T1 = []
arrays_T1 = []
for MC in T1:
    means_T1.append('results/means_L=20T=1.00MC=%s.txt' % MC)
    arrays_T1.append('results/arrays_L=20T=1.00MC=%s.txt' % MC)

# files with different MC-cycles with T = 2
T2 = ['1000', '2500', '5000', '6500', '8000', '10000']
means_T2 = []
arrays_T2 = []

for MC in T2:
    means_T2.append('results/means_L=20T=2.40MC=%s.txt' % MC)
    arrays_T2.append('results/arrays_L=20T=2.40MC=%s.txt' % MC)

#plot_expectation_c(means_T1, means_T2, L=20)
plot_arrays_c(arrays_T1, arrays_T2, L=20)

#
