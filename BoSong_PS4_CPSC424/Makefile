CC = mpicc
CFLAGS = -g -O3 -xHost -fno-alias -std=c99
FC = ifort
FFLAGS = -g -O3 -xHost -fno-alias

all: serial parallel
serial:	serial.o /home/fas/hpcprog/ahs3/cpsc424/utils/timing/timing.o
	$(CC) -o $@ $(CFLAGS) $^

parallel: parallel.o /home/fas/hpcprog/ahs3/cpsc424/utils/timing/timing.o
	$(CC) -o $@ $(CFLAGS) $^

fserial: fserial.o /home/fas/hpcprog/ahs3/cpsc424/utils/timing/timing.o
	$(FC) -o $@ $(FFLAGS) $^

.c.o:
	$(CC) $(CFLAGS) -c $<

.f.o:
	$(FC) $(FFLAGS) -c $<

clean:
	rm -f serial parallel fserial *.o
