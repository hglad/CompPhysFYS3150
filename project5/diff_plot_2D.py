import numpy as np
import matplotlib.pyplot as plt

u = np.loadtxt('2D_forward_dx=0.100.txt')


#print u
print u.shape   # u[time, position]
#print u

n = len(u[0,:])
x = np.linspace(0, 1, n)
#print x
print u[-2, :]
t_steps = len(u[:,0])/len(u[0,:])    # get number of time steps from matrix
#print t_steps
print 'loop'
for t in range(t_steps):
    start = t_steps*t
    end = start+n
    mat = u[start:end]
    
    print start, end
    print mat

    #plt.plot(u[])

plt.plot(u)
#plt.show()
