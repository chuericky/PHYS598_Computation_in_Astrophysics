### Plot of Q.4.

import numpy as np
import matplotlib.pyplot as plt

data_link = "/Users/rickyccy/Documents/Urbana-Champaign/Courses/PHYS598/ps10/q_4/data/"
fig_link = "/Users/rickyccy/Documents/Urbana-Champaign/Courses/PHYS598/ps10/q_4/fig/"

N = 32
N_str = str(N)
dt = 20. / N

founame = np.loadtxt(data_link + N_str + "_fourier.dat", unpack = True, usecols = (1,2,))

fou_freq = founame[0]
fou_value = np.abs(founame[1])

ana_freq = fou_freq

ana_value = np.pi**0.5 * np.exp(-np.power((np.pi * ana_freq), 2))

plt.subplot(1,1,1)
plt.scatter(fou_freq, fou_value)
plt.plot(ana_freq, ana_value, color = 'r')
plt.legend(['Analytic', 'Numerical'])
plt.title('Sample # = ' + N_str + ', ' + r'$\Delta$ = ' + str(dt), fontsize = 16)
plt.xlabel('Frequency, f', fontsize = 12)
plt.ylabel(r'$H(f)$', fontsize = 16)

"""
plt.subplot(2,1,2)
plt.scatter(fou_freq, fou_value)
plt.plot(ana_freq, ana_value, color = 'r')
plt.legend(['Analytic', 'Numerical'])
#plt.title('Sample # = ' + N_str, fontsize = 16)
plt.xlabel('Frequency, f', fontsize = 12)
plt.xlim((-1,1))
#plt.ylabel(r'$H(f)$', fontsize = 16)
"""
plt.savefig(fig_link + N_str + "_fourier.pdf")
plt.close()



invfouname = np.loadtxt(data_link + N_str + "_invfourier.dat", unpack = True, usecols = (1,2,))

inv_x = invfouname[0]
inv_value = invfouname[1]

ana_v = np.exp(-np.power(inv_x, 2))

plt.scatter(inv_x, inv_value)
plt.plot(inv_x, ana_v, color = 'r')
plt.legend(['Analytic', 'Numerical'])
plt.title('Sample # = ' + N_str + ', ' + r'$\Delta$ = ' + str(dt), fontsize = 16)
plt.xlabel('Time, t', fontsize = 12)
plt.ylabel(r'$h(t)$', fontsize = 16)

plt.savefig(fig_link + N_str + "_invfourier.pdf")
plt.close()
