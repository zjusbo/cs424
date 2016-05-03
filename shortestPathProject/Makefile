CC = icc
ICC = icc
CFLAGS = -g -O3 -std=c99 -openmp -xHost -fno-alias
FC = ifort
IFLAGS = -g -O3 -xHost -fno-alias -openmp -std=c99

BINARY=serial parallel_critical_loop parallel_critical_task parallel_node parallel_lock_dual parallel_lock_single
all: ${BINARY}

serial:	serial.o queue.o /home/fas/hpcprog/ahs3/cpsc424/utils/timing/timing.o
	$(CC) -o $@ $(CFLAGS) $^

parallel_node: parallel_node.o queue.o /home/fas/hpcprog/ahs3/cpsc424/utils/timing/timing.o
	$(ICC) -o $@ $(IFLAGS) $^ -lm

parallel_critical_task: parallel_critical_task.o queue.o /home/fas/hpcprog/ahs3/cpsc424/utils/timing/timing.o
	$(ICC) -o $@ $(IFLAGS) $^ -lm

parallel_critical_loop: parallel_critical_loop.o queue.o /home/fas/hpcprog/ahs3/cpsc424/utils/timing/timing.o
	$(ICC) -o $@ $(IFLAGS) $^ -lm

parallel_lock_dual: parallel_lock_dual.o queue.o /home/fas/hpcprog/ahs3/cpsc424/utils/timing/timing.o
	$(ICC) -o $@ $(IFLAGS) $^ -lm

parallel_lock_single: parallel_lock_single.o queue.o /home/fas/hpcprog/ahs3/cpsc424/utils/timing/timing.o
	$(ICC) -o $@ $(IFLAGS) $^ -lm

.c.o:
	$(CC) $(CFLAGS) -c $<

clean:
	rm -f ${BINARY} *.o
