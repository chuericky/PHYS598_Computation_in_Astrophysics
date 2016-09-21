import numpy as np
import matplotlib.pyplot as plt

datalink = "/Users/rickyccy/Documents/Urbana-Champaign/Courses/PHYS598/ps10/q_6/data/"
figlink = "/Users/rickyccy/Documents/Urbana-Champaign/Courses/PHYS598/ps10/q_6/fig/"

pendata = np.loadtxt(datalink + "1_pen_v.dat", unpack = True, usecols = (2,4,))

y = pendata[1]
x = pendata[0]

plt.subplot(2,1,1)
plt.plot(x, y)
plt.grid()
plt.title((r'One-sided power spectra of $\frac{d\theta}{dt}$, ' + r'$(q, b, \omega_0) = (0.5, 0.9, 2/3)$'))
plt.ylabel(r'$\ln P(\omega)$', fontsize = 16)
plt.subplot(2,1,2)
plt.plot(x, y)
plt.grid()
plt.xlim(0,5)
plt.ylabel(r'$\ln P(\omega)$', fontsize = 16)
plt.xlabel(r'$\omega / \omega_0$', fontsize = 16)
plt.savefig(figlink + "1_pen_v.pdf")
plt.close()

pendata = np.loadtxt(datalink + "1_pen_x.dat", unpack = True, usecols = (2,4,))

y = pendata[1]
x = pendata[0]


plt.subplot(2,1,1)
plt.plot(x, y)
plt.grid()
plt.title((r'One-sided power spectra of $\theta$, ' + r'$(q, b, \omega_0) = (0.5, 0.9, 2/3)$'))
plt.ylabel(r'$\ln P(\omega)$', fontsize = 16)
plt.subplot(2,1,2)
plt.plot(x, y)
plt.grid()
plt.xlim(0,5)
plt.ylabel(r'$\ln P(\omega)$', fontsize = 16)
plt.xlabel(r'$\omega / \omega_0$', fontsize = 16)
plt.savefig(figlink + "1_pen_x.pdf")
plt.close()