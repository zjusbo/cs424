#!/bin/bash
#PBS -l procs=8,tpn=4,mem=68gb,walltime=30:00
#PBS -q cpsc424
#PBS -j oe

cd $PBS_O_WORKDIR

module load Langs/Intel/14 MPI/OpenMPI/1.6.5

for input in actualdata1 actualdata2 actualdata3 actualdata4
do
  mpiexec -n 8 ./parallel < ./data/$input > ./${input}_c_parallel.out
done

#./fserial < /home/fas/hpcprog/ahs3/cpsc424/assignment4/data/testdata1 > ./testdata1_f.out
