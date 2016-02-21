#!/bin/bash
#PBS -l procs=8,tpn=8,mem=68gb
#PBS -l walltime=30:00
#PBS -N Mulmatrix
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
echo 1 node
for p in 8
do
  for N in 8000 12000
  do
    echo p = $p, N = $N
    time mpiexec -n $p parallel $N
  done
done

exit
