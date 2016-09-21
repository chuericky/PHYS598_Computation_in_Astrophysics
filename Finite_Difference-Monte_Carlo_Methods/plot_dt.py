import numpy as np
import matplotlib.pyplot as plt
import math
import scipy.special as sp

fname = np.loadtxt("data_dt.dat", unpack = True, usecols = (2,3))
fname2 = np.loadtxt("data_dt0.5.dat", unpack = True, usecols = (2,3))
fname3 = np.loadtxt("data_dt0.25.dat", unpack = True, usecols = (2,3))
fname4 = np.loadtxt("data_dt0.125.dat", unpack = True, usecols = (2,3))

x = fname[0]
u = fname[1]
x2 = fname2[0]
u2 = fname2[1]
x3 = fname3[0]
u3 = fname3[1]
x4 = fname4[0]
u4 = fname4[1]

dx = 0.02
dx2 = 0.02
dx3 = 0.02
dx4 = 0.02
dt = 0.02 * 0.02
dt2 = 0.5 * 0.02 * 0.02
dt3 = 0.25 * 0.02 * 0.02
dt4 = 0.125 * 0.02 * 0.02

xarr2 = np.arange(-1, -0.000001, dx)
xarr3 = np.arange(0, 1.000001, dx)
xarr4 = np.arange(0.02, 1.000001, dx)

uarr3 = np.power(xarr3, 1./3)
uarr2 = -np.power(xarr4, 1./3)[::-1]

xarr = np.append(xarr2, xarr3)
uarr = np.append(uarr2, uarr3)


xarr22 = np.arange(-1, -0.000001, dx2)
xarr23 = np.arange(0, 1.000001, dx2)
xarr24 = np.arange(0.01, 1.000001, dx2)

uarr23 = np.power(xarr23, 1./3)
uarr22 = -np.power(xarr24, 1./3)[::-1]

xarr_2 = np.append(xarr22, xarr23)
uarr_2 = np.append(uarr22, uarr23)


xarr32 = np.arange(-1, -0.000001, dx3)
xarr33 = np.arange(0, 1.000001, dx3)
xarr34 = np.arange(0.005, 1.000001, dx3)

uarr33 = np.power(xarr33, 1./3)
uarr32 = -np.power(xarr34, 1./3)[::-1]

xarr_3 = np.append(xarr32, xarr33)
uarr_3 = np.append(uarr32, uarr33)


xarr42 = np.arange(-1, -0.000001, dx4)
xarr43 = np.arange(0, 1.000001, dx4)
xarr44 = np.arange(0.0025, 1.000001, dx4)

uarr43 = np.power(xarr43, 1./3)
uarr42 = -np.power(xarr44, 1./3)[::-1]

xarr_4 = np.append(xarr42, xarr43)
uarr_4 = np.append(uarr42, uarr43)


plt.scatter(x, np.abs(np.divide((u - uarr),uarr)), color = 'g')
plt.scatter(x2, np.abs(np.divide((u2 - uarr_2),uarr_2)), color = 'k')
plt.scatter(x3, np.abs(np.divide((u3 - uarr_3),uarr_3)), color = 'r')
plt.scatter(x4, np.abs(np.divide((u4 - uarr_4),uarr_4)), color = 'c')
plt.legend(['Numerical','Numerical0.5','Numerical0.25','Numerical0.125'], loc = 'upper right')
plt.ylabel(r'$|\epsilon|$')
plt.title(r'Time = $100\Delta t$')
plt.xlabel(r'$x$')
plt.xlim((-1, 1))
#plt.ylim((-0.1,1.1))
#plt.show()
plt.savefig('fig/100_dtdifferr.png')
plt.close()