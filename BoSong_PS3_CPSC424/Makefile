TIMINGDIR = /home/fas/hpcprog/ahs3/cpsc424/utils/timing
CC = mpicc
CFLAGS = -g -O3 -xHost -fno-alias -std=c99 -I$(TIMINGDIR) -Wall

all: serial task2 task3s
task2: task2.c ${TIMINGDIR}/timing.o
task3: task3.c ${TIMINGDIR}/timing.o
serial:	serial.o matmul.o $(TIMINGDIR)/timing.o
	$(CC) -o $@ $(CFLAGS) $^

.c.o:
	$(CC) $(CFLAGS) -c $<

clean:
	rm -f serial parallel *.o
