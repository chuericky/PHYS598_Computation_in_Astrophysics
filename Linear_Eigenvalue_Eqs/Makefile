CC = gcc
CFLAGS = -O3
HOME_DIR = /Users/rickyccy/Documents/Urbana-Champaign/Courses/PHYS598/ps8/
LIBS = -lm

all: q1 q2

q1: q1.c nrutil.c nrutil.h nr.h tqli.c dpythag.c eigsrt.c tred2.c Makefile
	$(CC) $(CFLAGS) -o q1 q1.c nrutil.c tqli.c eigsrt.c dpythag.c tred2.c $(LIBS)

q2: q2.c nrutil.c nrutil.h nr.h tqli.c dpythag.c eigsrt.c tred2.c Makefile
	$(CC) $(CFLAGS) -o q2 q2.c nrutil.c tqli.c eigsrt.c dpythag.c tred2.c $(LIBS)

