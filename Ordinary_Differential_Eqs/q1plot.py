### Plotting of Q.1

import numpy as np
import glob, os
import matplotlib.pyplot as plt

link = "/Users/rickyccy/Documents/Urbana-Champaign/Courses/PHYS598/ps6/data/q1/"
figlink = "/Users/rickyccy/Documents/Urbana-Champaign/Courses/PHYS598/ps6/figure/q1/"

os.chdir(link)

fname = glob.glob("*")

for i in range(1):
    infile = np.loadtxt(fname[i], unpack = True)
    t = infile[0]
    theta = infile[1]
    omega = infile[2]
    plt.scatter(t, theta, marker = 'x')
    plt.xlabel('t', fontsize = 18)
    plt.ylabel(r'$\theta$', fontsize = 18)
    plt.xlim((0,100))
    plt.title(r'$(q, b, \omega_0) = (0.5, 0.9, 2/3). t = [0, 100]$', fontsize = 18)
    plt.grid()
    plt.plot()
    plt.savefig(figlink + 'a1.png')
    plt.close()
    
    plt.scatter(t, omega, marker = 'x')
    plt.xlabel('t', fontsize = 18)
    plt.ylabel(r'$\frac{d\theta}{dt}$', fontsize = 18)
    plt.xlim((0,100))
    plt.title(r'$(q, b, \omega_0) = (0.5, 0.9, 2/3). t = [0, 100]$', fontsize = 18)
    plt.grid()
    plt.plot()
    plt.savefig(figlink + 'a2.png')
    plt.close()
    
    plt.scatter(theta, omega, marker = 'x')
    plt.xlabel(r'$\theta$', fontsize = 18)
    plt.ylabel(r'$\frac{d\theta}{dt}$', fontsize = 18)
    plt.title(r'$(q, b, \omega_0) = (0.5, 0.9, 2/3). t = [0, 100]$', fontsize = 18)
    plt.grid()
    plt.plot()
    plt.savefig(figlink + 'a3.png')
    plt.close()

    infile2 = np.loadtxt(fname[i + 2], unpack = True)
    t2 = infile2[0]
    theta2 = infile2[1]
    omega2 = infile2[2]
    plt.scatter(theta2, omega2, marker = 'x')
    plt.xlabel(r'$\theta$', fontsize = 18)
    plt.ylabel(r'$\frac{d\theta}{dt}$', fontsize = 18)
    plt.title(r'$(q, b, \omega_0) = (0.5, 0.9, 2/3). t = [0, 1000]$', fontsize = 18)
    #plt.xlim((-30,40))
    plt.grid()
    plt.plot()
    plt.savefig(figlink + 'a4.png')
    plt.close()