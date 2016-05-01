#!/bin/bash
#PBS -l procs=8,tpn=8,mem=34gb
#PBS -l walltime=30:00
#PBS -N parallel_source
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

export OMP_NUM_THREADS=8
echo Long-C.10.0
./parallel_source_dynamic ch9-1.1/inputs/Long-C/Long-C.10.0.gr ch9-1.1/inputs/Long-C/Long-C.10.0.ss
echo Long-n.20.0
./parallel_source_dynamic ch9-1.1/inputs/Long-n/Long-n.20.0.gr ch9-1.1/inputs/Long-n/Long-n.20.0.ss
echo Random4-C.12.0
./parallel_source_dynamic ch9-1.1/inputs/Random4-C/Random4-C.12.0.gr ch9-1.1/inputs/Random4-C/Random4-C.12.0.ss
echo Random4-n.20.0
./parallel_source_dynamic ch9-1.1/inputs/Random4-n/Random4-n.20.0.gr ch9-1.1/inputs/Random4-n/Random4-n.20.0.ss
echo Square-C.12.0
./parallel_source_dynamic ch9-1.1/inputs/Square-C/Square-C.12.0.gr ch9-1.1/inputs/Square-C/Square-C.12.0.ss
echo Square-n.20.0
./parallel_source_dynamic ch9-1.1/inputs/Square-n/Square-n.20.0.gr ch9-1.1/inputs/Square-n/Square-n.20.0.ss
echo USA-road-d.NE
./parallel_source_dynamic ch9-1.1/inputs/USA-road-d/USA-road-d.NE.gr ch9-1.1/inputs/USA-road-d/USA-road-d.NE.ss
echo USA-road-d.NY
./parallel_source_dynamic ch9-1.1/inputs/USA-road-d/USA-road-d.NY.gr ch9-1.1/inputs/USA-road-d/USA-road-d.NY.ss
exit
