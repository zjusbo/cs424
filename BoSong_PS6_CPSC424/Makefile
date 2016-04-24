TIMING = /home/fas/hpcprog/ahs3/cpsc424/utils/timing/timing.o
CUDAPATH = /home/apps/fas/GPU/cuda_7.0.28

NVCC = $(CUDAPATH)/bin/nvcc

NVCCFLAGS = -I$(CUDAPATH)/include -O3

LFLAGS = -L$(CUDAPATH)/lib64 -lcuda -lcudart -lm

# Compiler-specific flags (by default, we always use sm_20)
GENCODE_SM20 = -gencode=arch=compute_20,code=\"sm_20,compute_20\"
GENCODE = $(GENCODE_SM20)

.SUFFIXES : .cu .ptx

BINARIES = matmul

matmul: matmul.o
	$(NVCC) $(GENCODE) $(LFLAGS) -o $@ $<

task1: task1.o
	$(NVCC) $(GENCODE) $(LFLAGS) -o $@ $<
task2: task2.o
	$(NVCC) $(GENCODE) $(LFLAGS) -o $@ $<
task3: task3.o
	$(NVCC) $(GENCODE) $(LFLAGS) -o $@ $<
multiGPU: multiGPU.o
	$(NVCC) $(GENCODE) $(LFLAGS) -o $@ $<
serial: serial.o $(TIMING)
	icc -std=c99 -g -O3 -fno-alias -o $@ $^
.c.o:
	icc -std=c99 -g -O3 -fno-alias -c $<
.cu.o:
	$(NVCC) $(GENCODE) $(NVCCFLAGS) -o $@ -c $<

clean:	
	rm -f *.o $(BINARIES)
