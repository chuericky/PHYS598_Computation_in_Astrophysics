CC = gcc
C++ = g++
CFLAGS = -O3
HOME_DIR = /Users/rickyccy/Documents/Urbana-Champaign/Courses/PHYS598/ps3/
LIBS = -lm

all: 1_trapezoidal 1_simpson 1_romberg 2_simpson 3_simpson2 2_simpson_cov 3_simpson 4_midpoint

1_trapezoidal: 1_trapezoidal.c Makefile
	$(CC) $(CFLAGS) -o 1_trapezoidal 1_trapezoidal.c $(LIBS)

1_simpson: 1_simpson.c Makefile
	$(CC) $(CFLAGS) -o 1_simpson 1_simpson.c $(LIBS)

1_romberg: 1_romberg.c Makefile
	$(CC) $(CFLAGS) -o 1_romberg 1_romberg.c $(LIBS)

2_simpson: 2_simpson.c Makefile
	$(CC) $(CFLAGS) -o 2_simpson 2_simpson.c $(LIBS)

2_simpson_cov: 2_simpson_cov.c Makefile
	$(CC) $(CFLAGS) -o 2_simpson_cov 2_simpson_cov.c $(LIBS)

3_simpson: 3_simpson.c Makefile
	$(CC) $(CFLAGS) -o 3_simpson 3_simpson.c $(LIBS)

3_simpson2: 3_simpson2.c Makefile
	$(CC) $(CFLAGS) -o 3_simpson2 3_simpson2.c $(LIBS)

4_midpoint: 4_midpoint.c Makefile
	$(CC) $(CFLAGS) -o 4_midpoint 4_midpoint.c $(LIBS)