CC = gcc
C++ = g++
CFLAGS = -O3
HOME_DIR = /Users/rickyccy/Documents/Urbana-Champaign/Courses/PHYS598/ps2/
LIBS = -lm

all: 1_bisection 1_newton 1_false 2_quadratic 2_perturb

1_bisection: 1_bisection.c Makefile
	$(CC) $(CFLAGS) -o 1_bisection 1_bisection.c $(LIBS)

1_newton: 1_newton.c Makefile
	$(CC) $(CFLAGS) -o 1_newton 1_newton.c $(LIBS)

1_false: 1_false.c Makefile
	$(CC) $(CFLAGS) -o 1_false 1_false.c $(LIBS)

2_quadratic: 2_quadratic.c Makefile
	$(CC) $(CFLAGS) -o 2_quadratic 2_quadratic.c $(LIBS)

2_perturb: 2_perturb.c Makefile
	$(CC) $(CFLAGS) -o 2_perturb 2_perturb.c $(LIBS)

3_mnewt: 3_mnewt.c nr3.h mnewt.h ludcmp.h Makefile
	$(C++) $(CFLAGS) -o 3_mnewt 3_mnewt.c $(LIBS)
