#!/bin/bash
#PBS -l procs=8,tpn=8,mem=34gb,walltime=30:00
#PBS -q cpsc424
#PBS -j oe

cd $PBS_O_WORKDIR

module load Langs/Intel/14

./serial < /home/fas/hpcprog/ahs3/cpsc424/assignment4/data/actualdata4 > ./actualdata4_c.out

#./fserial < /home/fas/hpcprog/ahs3/cpsc424/assignment4/data/testdata1 > ./testdata1_f.out
