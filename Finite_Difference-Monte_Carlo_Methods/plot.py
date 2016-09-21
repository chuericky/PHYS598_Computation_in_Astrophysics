import numpy as np
import matplotlib.pyplot as plt
import math
import scipy.special as sp

fname1 = np.loadtxt("/home/quantum-monkey/workspace/CPAcodes/ps9/data/p3data1.dat", unpack = True, usecols = (2,3))
fname2 = np.loadtxt("/home/quantum-monkey/workspace/CPAcodes/ps9/data/p3data2.dat", unpack = True, usecols = (2,3))
fname3 = np.loadtxt("/home/quantum-monkey/workspace/CPAcodes/ps9/data/p3data3.dat", unpack = True, usecols = (2,3))

x1 = fname1[0]
u1 = fname1[1]
x2 = fname2[0]
u2 = fname2[1]
x3 = fname3[0]
u3 = fname3[1]


dx = 0.002
dt = 0.1 * dx * dx

xarr = np.arange(-1.0, 1.0,dx)
uarr = []

time1 = 100 * dt
#uarr1 = 1/np.sqrt(1 + 4 * time1) * np.exp(-np.power(xarr, 2) / (1 + 4 * time1))
#uarr1 = sp.erf(xarr/(2.0*np.sqrt(time1)))+ np.exp(xarr+time1)*sp.erfc(np.sqrt(time1)+xarr/(2.0*np.sqrt(time1)))
time2 = 5000 * dt
#uarr2 = 1/np.sqrt(1 + 4 * time2) * np.exp(-np.power(xarr, 2) / (1 + 4 * time2))
uarr2 = sp.erf(xarr/(2.0*np.sqrt(time2)))+ np.exp(xarr+time2)*sp.erfc(np.sqrt(time2)+xarr/(2.0*np.sqrt(time2)))
time3 = 7000 * dt
#uarr3 = 1/np.sqrt(1 + 4 * time3) * np.exp(-np.power(xarr, 2) / (1 + 4 * time3))
uarr3 = sp.erf(xarr/(2.0*np.sqrt(time3)))+ np.exp(xarr+time3)*sp.erfc(np.sqrt(time3)+xarr/(2.0*np.sqrt(time3)))

#uarr = np.power(abs(xarr),1.0/3.0)
for x in xarr:
    if x<0 : 
	uarr.append(-np.power(abs(x),1.0/3.0))
    else :
	uarr.append(np.power(x,1.0/3.0))
#print uarr

plt.plot(x1, u1, color = 'b')
plt.scatter(xarr, uarr, color = 'g')
plt.legend(['Numerical','Analytic'])
plt.ylabel(r'$u(x, t)$')
plt.title(r'Time = $100\Delta t$')
plt.xlabel(r'$x$')
plt.xlim((-1, 1))
plt.savefig('fig/p3_100.png')
plt.close()

plt.plot(x2, u2, color = 'b')
plt.scatter(xarr, uarr, color = 'g')
plt.legend(['Numerical','Analytic'])
plt.ylabel(r'$u(x, t)$')
plt.title(r'Time = $5000\Delta t$')
plt.xlabel(r'$x$')
plt.xlim((-1, 1))
plt.savefig('fig/p3_5000.png')
plt.close()

plt.plot(x3, u3,color = 'b')
plt.scatter(xarr, uarr, color = 'g')
plt.legend(['Numerical','Analytic'])
plt.ylabel(r'$u(x, t)$')
plt.title(r'Time = $50000\Delta t$')
plt.xlabel(r'$x$')
plt.xlim((-1, 1))
plt.savefig('fig/p3_50000.png')
plt.close()


