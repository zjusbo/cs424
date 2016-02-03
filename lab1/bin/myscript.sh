#!/bin/bash
#PBS -l procs=8,tpn=8,mem=34gb,walltime=15:00
#PBS -N bo.song_lab1_q1
#PBS -q cpsc424
output=./lab1/output

echo test.1 >> $output
./lab1/test.1 >> $output

echo test.2 >> $output
./lab1/test.2 >> $output

echo test.3 >> $output
./lab1/test.3 >> $output

echo test.4 >> $output
./lab1/test.4 >> $output
exit
