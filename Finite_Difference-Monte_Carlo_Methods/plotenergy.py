import numpy as np
import matplotlib.pyplot as plt


fname = np.loadtxt("/home/quantum-monkey/workspace/CPAcodes/ps9/data/p2energy.dat", unpack = True, usecols = (1,2))

t = fname[0]
eps = fname[1]

plt.scatter(t, eps)
plt.ylabel(r'$\epsilon$')
plt.title('Energy Conservation')
plt.xlabel(r'$t$')
plt.xlim(-1,10)
plt.ylim(-0.01,0.2)
plt.savefig('fig/p2energy.png')
plt.close()

