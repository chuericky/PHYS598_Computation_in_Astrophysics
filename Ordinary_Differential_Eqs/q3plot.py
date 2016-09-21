import numpy as np
import glob, os
import matplotlib.pyplot as plt

datalink = "/Users/rickyccy/Documents/Urbana-Champaign/Courses/PHYS598/ps6/data/q3/"
figlink = "/Users/rickyccy/Documents/Urbana-Champaign/Courses/PHYS598/ps6/ps6_tex/q3/"

os.chdir(datalink)

x = np.linspace(0, 1, 1001)
y = np.cos(np.pi * x / 2.) + 2 * np.sin(np.pi * x / 2.) - 1

fname = glob.glob('*')

for i in range(len(fname)):
    infile = np.loadtxt(fname[i], unpack = True, usecols = (0, 1,))
    t = infile[0]
    u = infile[1]
    plt.scatter(t, u, marker = 'x', color = 'r')
    plt.scatter(x, y, marker = '.', color = 'k')
    plt.legend(['Numerical','Analytic'], loc = 'lower right')
    plt.ylabel(r'$u(x)$', fontsize = 16)
    plt.xlabel(r'$x$', fontsize = 16)
    plt.xlim((0, 1))
    plt.ylim((0,1.4))
    plt.grid()
    plt.plot()
    plt.savefig(figlink + 'q3.png')
    plt.close()

    z = np.divide((u - y), y)
    plt.scatter(t, z, marker = 'x', color = 'b')
    plt.ylabel(r'$\epsilon$', fontsize = 16)
    plt.xlabel(r'$x$', fontsize = 16)
    plt.xlim((-0.1, 1.1))
    plt.ylim((-0.00005, 0.00005))
    plt.grid()
    plt.tight_layout()
    plt.plot()
    plt.savefig(figlink + 'q3err1.png')
    plt.close()