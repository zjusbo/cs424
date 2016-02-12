TIMINGDIR = /home/fas/hpcprog/ahs3/cpsc424/utils/timing
CC = mpicc
CFLAGS = -g -O3 -xHost -fno-alias -std=c99 -I$(TIMINGDIR) -Wall

all: serial parallel
parallel: parallel.c ${TIMINGDIR}/timing.o
serial:	serial.o matmul.o $(TIMINGDIR)/timing.o
	$(CC) -o $@ $(CFLAGS) $^

.c.o:
	$(CC) $(CFLAGS) -c $<

clean:
	rm -f serial parallel *.o
