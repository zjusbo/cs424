#!/bin/bash
#PBS -l procs=8,tpn=8,mem=34gb
#PBS -l walltime=30:00
#PBS -N single_lock
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
export OMP_NUM_THREADS=1
./parallel_lock_single Random4-n.10.0.gr Random4-n.10.0.ss
export OMP_NUM_THREADS=2
./parallel_lock_single Random4-n.10.0.gr Random4-n.10.0.ss
export OMP_NUM_THREADS=4
./parallel_lock_single Random4-n.10.0.gr Random4-n.10.0.ss
export OMP_NUM_THREADS=8
./parallel_lock_single Random4-n.10.0.gr Random4-n.10.0.ss
exit
