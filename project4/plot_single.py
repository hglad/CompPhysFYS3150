import numpy as np
import matplotlib.pyplot as plt
import glob

def plot_array(filename):
    E = np.loadtxt(filename, usecols=0, dtype='float')
    numMC = len(E)

    x = np.linspace(0, len(E)-1, numMC)

    plt.plot(x, E)
    plt.title('T = 1.0 [kT/J], L = 20')
    plt.grid('on'); plt.ylabel('E [$JL^2$]'); plt.xlabel('MC-cycle')
    plt.show()

def plot_c(files):
    n = len(files)
    print n
    E = np.zeros(n); absM = np.zeros(n)
    M2 = np.zeros(n); C_V = np.zeros(n)
    chi = np.zeros(n); counts = np.zeros(n)
    MC = np.zeros(n)
    #plot_array(filename)
    i = 0
    for file in files:
        E[i], absM[i], M2[i], C_V[i], chi[i], counts[i], MC[i] = np.loadtxt(file, usecols=(0,1,2,3,4,5,6), unpack=True)
        i += 1

    plt.title('T=1')
    plt.plot(MC, counts, 'o')
    plt.xlabel('MC-cycles'); plt.ylabel('Accepted states')
    plt.grid('on')

# files with different MC-cycles with T = 1
T1 = ['1000', '5000', '10000']
files_T1 = []
for MC in T1:
    files_T1.append('results/means_L=20T=1.00MC=%s.txt' % MC)
#files_T1 = ['results/means_L=20T=1.00MC=1000.txt', 'results/means_L=20T=1.00MC=5000.txt']
plot_c(files_T1)
plt.show()
