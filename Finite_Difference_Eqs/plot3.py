import matplotlib.pyplot as plt
import numpy as np
import glob, os

datalink = "/Users/rickyccy/Documents/Urbana-Champaign/Courses/PHYS598/ps5/data/"
figlink = "/Users/rickyccy/Documents/Urbana-Champaign/Courses/PHYS598/ps5/figure/"

os.chdir(datalink)

fname = glob.glob('c*')

g = (5./3)**0.5

for i in range(2,3):
    
    infile = np.loadtxt(fname[i], unpack = True)
    infile2 = np.loadtxt(fname[i + 1], unpack = True)

    r = infile[0]
    R = infile[1]
    rho = infile[2]
    P = infile[3]
    u = infile[4]
    t = infile[5][0]
    R2 = infile2[1]
    rho2 = infile2[2]
    P2 = infile2[3]
    u2 = infile2[4]
    t2 = infile2[5][0]
    
    plt.scatter(R, u, marker = '*', color = 'k')
    plt.scatter(R2, u2, marker = 'x', color = 'b')
    plt.legend(['t = 0.125040', 't = 0.239119'], loc = 'lower left')
    #plt.ylim((-0.1, 0.5))
    plt.vlines(x = 1.9493736000000002, ymin = 0, ymax = 11.62, color = 'k')
    plt.hlines(y = 0.0, xmin = 1.9493736000000002, xmax = 6, color = 'k')
    plt.hlines(y = 11.62, xmin = 0., xmax = 1.9493736000000002, color = 'k')
    plt.vlines(x = 3.72786521, ymin = 0, ymax = 11.62, color = 'b')
    plt.hlines(y = 0.0, xmin = 3.72786521, xmax = 6, color = 'b')
    plt.hlines(y = 11.62, xmin = 0, xmax = 3.72786521, color = 'b')
    plt.xlim((0, 5))
    plt.ylabel('u', fontsize = 20)
    plt.xlabel('R', fontsize = 20)
    plt.title(r'$U = 10C_0$')
    #plt.show()
    plt.savefig(figlink + 'C10_u.png')
    plt.close()

    
    plt.scatter(R, rho, marker = '*', color = 'k')
    plt.scatter(R2, rho2, marker = 'x', color = 'b')
    plt.legend(['t = 0.125040', 't = 0.239119'], loc = 'lower left')
    #plt.ylim((-0.1, 0.5))
    plt.vlines(x = 1.9493736000000002, ymin = 1, ymax = 4, color = 'k')
    plt.hlines(y = 1.0, xmin = 1.9493736000000002, xmax = 6, color = 'k')
    plt.hlines(y = 4, xmin = 0., xmax = 1.9493736000000002, color = 'k')
    plt.vlines(x = 3.72786521, ymin = 1, ymax = 4, color = 'b')
    plt.hlines(y = 1.0, xmin = 3.72786521, xmax = 6, color = 'b')
    plt.hlines(y = 4, xmin = 0, xmax = 3.72786521, color = 'b')
    plt.xlim((0, 5))
    plt.ylabel(r'$\rho$', fontsize = 20)
    plt.xlabel('R', fontsize = 20)
    plt.title(r'$U = 10C_0$')
    #plt.show()
    plt.savefig(figlink + 'C10_rho.png')
    plt.close()
    
    

    plt.scatter(R, P, marker = '*', color = 'k')
    plt.scatter(R2, P2, marker = 'x', color = 'b')
    plt.legend(['t = 0.125040', 't = 0.239119'], loc = 'lower left')
    #plt.ylim((-0.1, 0.5))
    plt.vlines(x = 1.9493736000000002, ymin = 1, ymax = 182, color = 'k')
    plt.hlines(y = 1.0, xmin = 1.9493736000000002, xmax = 6, color = 'k')
    plt.hlines(y = 182, xmin = 0., xmax = 1.9493736000000002, color = 'k')
    plt.vlines(x = 3.72786521, ymin = 1, ymax = 182, color = 'b')
    plt.hlines(y = 1.0, xmin = 3.72786521, xmax = 6, color = 'b')
    plt.hlines(y = 182, xmin = 0, xmax = 3.72786521, color = 'b')
    plt.xlim((0, 5))
    plt.ylabel('P', fontsize = 20)
    plt.xlabel('R', fontsize = 20)
    plt.title(r'$U = 10C_0$')
    #plt.show()
    plt.savefig(figlink + 'C10_P.png')
    plt.close()