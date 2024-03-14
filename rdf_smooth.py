import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

data = np.loadtxt(fname = 'rad1')
with open('atoms') as f:
    lines = f.readlines()

atoms = []
for i,x in enumerate(lines):
    atoms.append(lines[i].strip('\n').strip('   '))

r = data[:,0]
for i in range(1,np.size(data[0,:])):
    plt.plot(r,data[:,i],color=cm.turbo(i/np.size(data[0,:])))

plt.legend(atoms)
plt.xlabel('Distance r, ($\AA$)',size=13)
plt.ylabel('g(r)',size=13)

plt.show()
