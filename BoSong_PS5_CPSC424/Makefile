# This makes both a sequential and parallel program to display the Mandelbrot Set
# Use module load Langs/Intel/14 MPI/OpenMPI

TIMING = /home/fas/hpcprog/ahs3/cpsc424/utils/timing/timing.o
EXECUTABLES = loop serial collapse_loop task
CC = icc
CFLAGS = -g -O3 -openmp -xHost -fno-alias -std=c99


# YOU HAVE TO PROVIDE YOUR OWN SOURCE FILE HERE, OF COURSE! (SORRY!!)
all: $(EXECUTABLES)

drand.o: drand.c
	$(CC) $(CFLAGS) -c $<

.c.o:
	$(CC) $(CFLAGS) -c $<

serial: serial.o /home/fas/hpcprog/ahs3/cpsc424/utils/timing/timing.o
	$(CC) $(CFLAGS) -openmp -o $@ $^ -lm
collapse_loop: collapse_loop.o drand.o $(TIMING)
	$(CC) $(CFLAGS) -openmp -o $@ $^ -lm
task: task.o drand.o $(TIMING)
	$(CC) $(CFLAGS) -openmp -o $@ $^ -lm
loop: loop.o drand.o /home/fas/hpcprog/ahs3/cpsc424/utils/timing/timing.o
	$(CC) $(CFLAGS) -openmp -o $@ $^ -lm
clean:
	rm -f $(EXECUTABLES) *.o
