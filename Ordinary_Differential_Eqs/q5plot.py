import numpy as np
import glob, os
import matplotlib.pyplot as plt

datalink = "/Users/rickyccy/Documents/Urbana-Champaign/Courses/PHYS598/ps6/data/q5/"
figlink = "/Users/rickyccy/Documents/Urbana-Champaign/Courses/PHYS598/ps6/ps6_tex/q5/"

os.chdir(datalink)

fname = glob.glob('*')

for i in range(len(fname) - 1, len(fname) - 0):
    infile = np.loadtxt(fname[i], unpack = True, usecols = (0, 1, 3, 4))
    t = infile[0]
    u = infile[1]
    E = infile[2]
    V = infile[3]
    plt.plot(t, V, linestyle = '-.', color = 'r')
    plt.plot(t, E, linestyle = '--', color = 'b')
    plt.plot(t, u, linestyle = '-', color = 'k')
    plt.legend([r'$V(x)$',r'$E_n$', r'$u(x)$'], loc = 'lower right')
    plt.xlabel(r'$x$', fontsize = 16)
    plt.xlim((-5, 5))
    plt.ylim((-3.1,3.1))
    plt.title(r'$E_2 = \frac{5}{2}$', fontsize = 20, y = 1.03)
    plt.grid()
    plt.plot()
    #plt.show()
    plt.savefig(figlink + 'q5E3.png')
    plt.close()