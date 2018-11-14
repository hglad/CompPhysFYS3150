import numpy as np
import matplotlib.pyplot as plt
import glob

# Find files with results for different temperatures
files = glob.glob('ising_arrays_*.txt')
temps = []

for file in files:
    temps.append(float(file[13:16]))    # extract temperature from file name

# Sort temperatures because the 'files' list is not sorted correctly
temps = np.sort(temps)

for T in temps:
    file = ('ising_arrays_%s00000.txt' % str(T))
    E, M = np.loadtxt(file, usecols=(0,1), unpack=True, dtype='float')
    numMC = len(E)
    #E = np.loadtxt("ising_data.txt")
    #print E
    E2 = np.dot(E,E)
    M2 = np.dot(M,M)

    #Cv = (E2 - E)/(T*T*numMC)
    #Chi = (M2 - M)/(T*numMC)
    """
    for line in infile:
        cols = line.split(' ')
        E.append(float(cols[0]))
        M.append(float(cols[1]))
    """
    x = np.linspace(0, len(E)-1, numMC)    # array for MC-cycles

    plt.figure()
    plt.title('T = %1.2f [kT/J]' % T)
    plt.plot(x, E)
    plt.grid('on'); plt.ylabel('E [$Js^2$]'); plt.xlabel('MC-cycle')
    """
    plt.figure()
    plt.plot(x, M)
    plt.grid('on'); plt.ylabel('M'); plt.xlabel('MC-cycle')

    plt.figure()
    plt.plot(x, Cv)
    plt.grid('on'); plt.ylabel('Cv')

    plt.figure()
    plt.plot(x, Chi)
    plt.grid('on'); plt.ylabel('Chi')
    """
    # Analytical values
    E_ = -np.sinh(8)*8 / (3 + np.cosh(8))
    #Cv_ = T**(-2) * (64*np.cosh(8)*(3 + np.cosh(8)) - 8*np.sinh(8)*8*np.sinh(8)) / (3 + np.cosh(8))**2
    Cv_ = 64/(T*T) * (np.cosh(8)*(3 + np.cosh(8)) - (np.sinh(8)**2))/(3 + np.cosh(8))**2
    M2_ = (8*np.exp(8) + 8)/(3 + np.cosh(8))
    absM_ = (2*np.exp(8) + 4)/(3 + np.cosh(8))
    chi_ = M2_/T

    #print E_, absM_, M2_, Cv_, chi_

plt.show()
