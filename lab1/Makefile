all: q1 q2

q1: q1.c timing.o
	icc main.c dummy.c timing.o q1.c -o q1.1 -g -O0 -fno-alias -std=c99
	icc main.c dummy.c timing.o q1.c -o q1.2 -g -O1 -fno-alias -std=c99
	icc main.c dummy.c timing.o q1.c -o q1.3 -g -O3 -no-vec -no-simd -fno-alias -std=c99
	icc main.c dummy.c timing.o q1.c -o q1.4 -g -O3 -xHost -fno-alias -std=c99

q2: q2.c timing.o
	icc main.c dummy.c timing.o q2.c -o q2.1 -g -O0 -fno-alias -std=c99
	icc main.c dummy.c timing.o q2.c -o q2.2 -g -O1 -fno-alias -std=c99
	icc main.c dummy.c timing.o q2.c -o q2.3 -g -O3 -no-vec -no-simd -fno-alias -std=c99
	icc main.c dummy.c timing.o q2.c -o q2.4 -g -O3 -xHost -fno-alias -std=c99

timing.o:
	cp ~ahs3/cpsc424/utils/timing/timing.o ./
clean:
	rm -f *.1 *.2 *.3 *.4
