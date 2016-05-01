CC = icc
ICC = icc
CFLAGS = -g -O3 -std=c99 -openmp -xHost -fno-alias
FC = ifort
IFLAGS = -g -O3 -xHost -fno-alias -openmp -std=c99

all: serial parallel 

serial:	serial.o queue.o /home/fas/hpcprog/ahs3/cpsc424/utils/timing/timing.o
	$(CC) -o $@ $(CFLAGS) $^

parallel: parallel.o queue.o /home/fas/hpcprog/ahs3/cpsc424/utils/timing/timing.o
	$(ICC) -o $@ $(IFLAGS) $^ -lm

parallel_source: parallel_source.o queue.o /home/fas/hpcprog/ahs3/cpsc424/utils/timing/timing.o
	$(ICC) -o $@ $(IFLAGS) $^ -lm

parallel_2: parallel_2.o queue.o /home/fas/hpcprog/ahs3/cpsc424/utils/timing/timing.o
	$(ICC) -o $@ $(IFLAGS) $^ -lm
parallel_array: parallel_array.o queue.o /home/fas/hpcprog/ahs3/cpsc424/utils/timing/timing.o
	$(ICC) -o $@ $(IFLAGS) $^ -lm
parallel_lock: parallel_lock.o queue.o /home/fas/hpcprog/ahs3/cpsc424/utils/timing/timing.o
	$(ICC) -o $@ $(IFLAGS) $^ -lm
parallel_lock_single: parallel_lock_single.o queue.o /home/fas/hpcprog/ahs3/cpsc424/utils/timing/timing.o
	$(ICC) -o $@ $(IFLAGS) $^ -lm

fserial: fserial.o /home/fas/hpcprog/ahs3/cpsc424/utils/timing/timing.o
	$(FC) -o $@ $(FFLAGS) $^

.c.o:
	$(CC) $(CFLAGS) -c $<

clean:
	rm -f serial parallel fserial parallel_lock parallel_loop parallel_2 parallel_lock_single parallel_source parallel_array *.o
