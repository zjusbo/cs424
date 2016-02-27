#!/bin/bash
#PBS -l procs=8,tpn=2,mem=136gb
#PBS -l walltime=30:00
#PBS -N task3
#PBS -r n
#PBS -j oe
#PBS -q cpsc424

# Load necessary module files
module load Langs/Intel/14 MPI/OpenMPI/1.6.5

# Print initial working directory
pwd

# Change to submission directory
cd $PBS_O_WORKDIR

# Now print the current working directory
pwd

# Print the node list
cat $PBS_NODEFILE

# Run the program 3 times
echo 4 nodes, proc = 8, tpn = 2
for p in 8
do
  for N in 8000 12000
  do
    echo p = $p, N = $N
    time mpiexec -n $p task4 $N
  done
done

exit
