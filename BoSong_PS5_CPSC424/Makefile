# This makes both a sequential and parallel program to display the Mandelbrot Set
# Use module load Langs/Intel/14 MPI/OpenMPI

EXECUTABLES = mandomp
CC = icc
CFLAGS = -g -O3 -openmp -xHost -fno-alias -std=c99

mandomp: mandomp.o drand.o /home/fas/hpcprog/ahs3/cpsc424/utils/timing/timing.o
	$(CC) $(CFLAGS) -openmp -o $@ $^ -lm

# YOU HAVE TO PROVIDE YOUR OWN SOURCE FILE HERE, OF COURSE! (SORRY!!)
mandomp.o: mandomp.c
	$(CC) $(CFLAGS) -c $<

drand.o: drand.c
	$(CC) $(CFLAGS) -c $<

.c.o:
	$(CC) $(CFLAGS) -c $<

serial: serial.o /home/fas/hpcprog/ahs3/cpsc424/utils/timing/timing.o
	gcc -g -O3 -std=c99 -o $@ $^
serial.o: serial.c
	gcc -g -O3 -std=c99 -c $<

clean:
	rm -f $(EXECUTABLES) *.o
