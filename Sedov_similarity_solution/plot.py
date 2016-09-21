### Plotting.

import numpy as np
import glob
import matplotlib.pyplot as plt

infile = np.loadtxt('g_5_3.dat', unpack = True)
eta = infile[0]
rho = infile[1]
v = infile[2]
P = infile[3]
T = infile[4]
M = infile[5]

plt.plot(eta, rho, color = 'k', linestyle = '-.')
plt.plot(eta, v, color = 'k', linestyle = '--')
plt.plot(eta, P, color = 'k', linestyle = '-')
plt.plot(eta, T, color = 'k', linestyle = ':')
plt.ylim((-0.5,2))
plt.title(r'$\gamma = 5/3$')
plt.xlabel(r'$\eta = \xi / \xi_0$', fontsize = 16)
plt.legend((r'$\rho / \rho_1$', r'$v / v_1$', r'$P / P_1$', r'$T / T_1$'), loc = 'best')
plt.savefig('g_5_3_1.png')
plt.close()

plt.plot(eta, M, color = 'k')
plt.title(r'$\gamma = 5/3$')
plt.xlabel(r'$\eta = \xi / \xi_0$', fontsize = 16)
plt.ylabel(r'$M(\eta)$')
plt.ylim((-0.5,2))
plt.savefig('g_5_3_mass.png')
plt.close()