#!/bin/bash
#PBS -l procs=8,tpn=8,mem=34gb,walltime=15:00
#PBS -q cpsc424
echo q1.1
./lab1/q1.1
echo q1.2
./lab1/q1.2
echo q1.3
./lab1/q1.3
echo q1.4
./lab1/q1.4

for filename in q2.1 q2.2 q2.3 q2.4
do
  echo $filename
  for i in {3..24}
  do 
    ./lab1/$filename $i
  done
done
