# This makes both a sequential and parallel program to display the Mandelbrot Set
# Use module load Langs/Intel/14 MPI/OpenMPI

CC = icc
CFLAGS = -g -O3 -xHost -fno-alias
EXECUTABLES = mandseq mand mandi
all: $(EXECUTABLES)
mandseq: mandseq.c /home/fas/hpcprog/ahs3/cpsc424/utils/timing/timing.o
	$(CC) $(CFLAGS) -o $@ $^ -lX11 -lm
mand: mand.c /home/fas/hpcprog/ahs3/cpsc424/utils/timing/timing.o
	mpicc $(CFLAGS) -o $@ $^ -lX11 -lm
mandi: mandi.c /home/fas/hpcprog/ahs3/cpsc424/utils/timing/timing.o
	mpicc $(CFLAGS) -o $@ $^ -lX11 -lm
clean:
	rm -f $(EXECUTABLES) *.o
