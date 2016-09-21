### Plotting of Q.1

import numpy as np
import glob, os
import matplotlib.pyplot as plt

link = "/Users/rickyccy/Documents/Urbana-Champaign/Courses/PHYS598/ps6/data/q2/"
figlink = "/Users/rickyccy/Documents/Urbana-Champaign/Courses/PHYS598/ps6/ps6_tex/q2/"

os.chdir(link)

fname = glob.glob("stiff*")

for i in range(1,2):
    print fname[i]
    infile = np.loadtxt(fname[i], unpack = True)
    x = infile[0]
    y1 = infile[1]
    y2 = infile[2]
    y3 = infile[3]
    plt.scatter(x, y1, marker = '.')
    plt.scatter(x, y2, marker = '+')
    plt.scatter(x, y3, marker = 'x')
    plt.legend([r'$y_1(x)$', r'$y_2(x)$', r'$y_3(x)$'])
    plt.xlabel('x', fontsize = 18)
    plt.xlim((0,1000))
    plt.ylim((-0.2, 3.0))
    plt.yscale('log')
    plt.title(r'$\epsilon = 0.05$', fontsize = 18)
    plt.grid()
    plt.plot()
    plt.savefig(figlink + '0.05_log.png')
    plt.close()
