#!/bin/bash
#PBS -l procs=8,tpn=8,mem=34gb
#PBS -l walltime=30:00
#PBS -N Lab2_task2
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
# echo 1024 1024 1024
# ./serial 1024 1024 1024

# echo 2048 2048 2048
# ./serial 2048 2048 2048

# echo 4096 4096 4096
# ./serial 4096 4096 4096

# echo 8192 8192 8192
# ./serial 8192 8192 8192


# echo 1024 1024 8192
# ./serial 1024 8192 1024

# echo  8192 8192 1024
# ./serial 8192 1024 8192

echo  8192 1024 8192
./serial 8192 8192 1024

exit
