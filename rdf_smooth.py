import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt(fname = 'rad1')
with open('atoms') as f:
    lines = f.readlines()

atoms = []
for i,x in enumerate(lines):
    atoms.append(lines[i].strip('\n').strip('   '))

r = data[:,0]
g1 = data[:,1]
g2 = data[:,3]
g3 = data[:,5]
g4 = data[:,7]
g5 = data[:,9]
g6 = data[:,11]

plt.plot(r,g1)
plt.plot(r,g2)
plt.plot(r,g3)
plt.plot(r,g4)
plt.plot(r,g5)
plt.plot(r,g6)
#plt.legend(atoms,ncol=3,loc='center',bbox_to_anchor=(0.5, 1.05),fancybox=True, shadow=True)
plt.legend(atoms)
plt.xlabel('Distance ($\AA$)')
plt.ylabel('g(r)')

plt.show()

# data2 = np.loadtxt(fname = 'rad2')
# r2 = data2[:,0]
# g21 = data2[:,1]
# g22 = data2[:,3]
# g23 = data2[:,5]
# g24 = data2[:,7]
# g25 = data2[:,9]
# g26 = data2[:,11]

# plt.plot(r2,g21)
# plt.plot(r2,g22)
# plt.plot(r2,g23)
# plt.plot(r2,g24)
# plt.plot(r2,g25)
# plt.plot(r2,g26)
# #plt.legend(atoms,ncol=3,loc='center',bbox_to_anchor=(0.5, 1.05),fancybox=True, shadow=True)
# plt.legend(atoms)
# plt.xlabel('Distance ($\AA$)')
# plt.ylabel('g(r)')
# #plt.title('2200 K')
# #plt.savefig('/users/leslie/desktop/FeOOH-paper/rdf2000.eps',dpi=1000,format='eps')

# plt.show()
