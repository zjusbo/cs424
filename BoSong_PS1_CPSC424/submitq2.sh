#!/bin/bash
#PBS -l procs=8,tpn=8,mem=34gb,walltime=15:00
#PBS -q cpsc424
for i in {3..24}
do
~/cs424/lab1/q2 $i
done 
