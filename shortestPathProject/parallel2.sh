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

for i in 1 2 4 8
do
	export OMP_NUM_THREADS=$i
	echo USA-road-d.NE
	./serial ch9-1.1/inputs/USA-road-d/USA-road-d.NE.gr ch9-1.1/inputs/USA-road-d/USA-road-d.NE.ss
done
exit

echo USA-road-d.NY
./serial ch9-1.1/inputs/USA-road-d/USA-road-d.NY.gr ch9-1.1/inputs/USA-road-d/USA-road-d.NY.ss
exit