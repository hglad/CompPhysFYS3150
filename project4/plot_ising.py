import numpy as np
import matplotlib.pyplot as plt
import glob

def read_filename(filename):# read dimension, temperature and return as strings
    cols = filename.split('=')
    L_str = (cols[1])[:-1]
    T_str = (cols[2])[:-4]

    return L_str, T_str

def plot_c():
    numMC = len(E)
    x = np.linspace(0, len(E)-1, numMC)    # array for MC-cycles
    plt.plot()

def analytical_values(T):
    E_ = -np.sinh(8)*8 / (3 + np.cosh(8))
    Cv_ = 64/(T*T) * (np.cosh(8)*(3 + np.cosh(8)) - (np.sinh(8)**2))/(3 + np.cosh(8))**2
    M2_ = (8*np.exp(8) + 8)/(3 + np.cosh(8))
    absM_ = (2*np.exp(8) + 4)/(3 + np.cosh(8))
    chi_ = M2_/T
    print E_, Cv_, M2_, absM_, chi_

def prob_state(L):       # Find probable state for L=20
    files = glob.glob('results/ising_arrays_*.txt')
    temps = []

    for file in files:
        L_str, T_str = read_filename(file)
        temps.append(float(T_str))

    # Sort temperatures because the 'files' list is not sorted correctly
    temps = np.sort(temps)
    for T in temps:
        file = ('results/ising_arrays_%s.txt' % str(T))
        E = np.loadtxt(file, usecols=0, dtype='float')

        numMC = len(E)
        cut_off = int(numMC/10)     # du not plot first 10% of MC-cycles
        plt.hist(E[cut_off:-1], bins=10)
        print len(E[cut_off:-1])
        plt.title("T=%1.1f, %1.0f MC-cycles" % (T, numMC))
        plt.show()

def plot_expect_vals(L, T_start, T_end, T_step=0):
    i = 0
    files = glob.glob('results/ising_means_L=%sT*.txt' % str(L))
    temps = []
    temp_dict = {}

    for file in files:
        L_str, T_str = read_filename(file)
        if ( T_start <= float(T_str) <= T_end):
            temps.append(float(T_str))
            temp_dict[float(T_str)] = T_str # dictionary to get corresponding string
    if len(files) == 0:
        print "No files found"

    temps = np.sort(temps)
    # no step value specified, use all available files
    if T_step == 0:
        # Sort temperatures because the 'files' list is not sorted correctly

        E = np.zeros(len(temps)); absM = np.zeros(len(temps))
        M2 = np.zeros(len(temps)); C_V = np.zeros(len(temps))
        chi = np.zeros(len(temps));

        for T in temps:
            file = ('results/ising_means_L=%sT=%s.txt' % (L_str, temp_dict[T]))
            E[i], absM[i], M2[i], C_V[i], chi[i] = np.loadtxt(file, usecols=(0,1,2,3,4), unpack=True)
            i += 1

    # Use files with a specific temperature step length
    if T_step != 0:
        temps_ = np.arange(T_start, T_end+T_step, T_step)
        i = 0
        E = np.zeros(len(temps_)); absM = np.zeros(len(temps_))
        M2 = np.zeros(len(temps_)); C_V = np.zeros(len(temps_))
        chi = np.zeros(len(temps_))

        for T in temps_:
            T = round(T, 4)             # round float to match dict. entry
            file = ('results/ising_means_L=%sT=%s.txt' % (L_str, temp_dict[T]))
            E[i], absM[i], M2[i], C_V[i], chi[i] = np.loadtxt(file, usecols=(0,1,2,3,4), unpack=True)
            i += 1

    plt.plot(temps, E)
    plt.grid('on'); plt.xlabel('T'); plt.ylabel('$\\chi$')
    plt.show()


def plot_energy_mag(L, T_start, T_end):
    # Find files with results for different temperatures
    files = glob.glob('results/ising_arrays_L=%sT*.txt' % str(L))
    temps = []
    temp_dict = {}
    if len(files) == 0:
        print "No files found"

    for file in files:
        L_str, T_str = read_filename(file)
        if ( T_start <= float(T_str) <= T_end):
            temps.append(float(T_str))
            temp_dict[float(T_str)] = T_str # dictionary to get corresponding string

    # Sort temperatures because the 'files' list is not sorted correctly
    temps = np.sort(temps)
    for T in temps:
        file = ('results/ising_arrays_L=%sT=%s.txt' % (L_str, temp_dict[T]))
        numMC = np.loadtxt(file, usecols=0)
        E, M = np.loadtxt(file, usecols=(0,1), skiprows=1, unpack=True, dtype='float')

        numMC = len(E)
        x = np.linspace(0, len(E)-1, numMC)    # array for MC-cycles
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

        plt.figure()
        plt.title('T = %1.2f [kT/J], L = %1.0f' % (T, L))
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

        #print E_, absM_, M2_, Cv_, chi_

    plt.show()

#read_filename('ising_means_L=20T=1.400000.txt')

#plot_expect_vals(40, 2.2, 2.4, 0.01)
plot_energy_mag(40,1,1)         # dimension, T_start, T_end to plot for
#prob_state(20)
