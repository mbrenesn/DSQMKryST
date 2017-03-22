#!/bin/bash
#PBS -q regular
#PBS -l nodes=1:ppn=1
#PBS -l walltime=01:00:00
#PBS -T flush_cache

cd $PBS_O_WORKDIR

# Load your Intel MPI/MKL modules!
module load intel/pe-xe-2017--binary
module load intelmpi/2017--binary 
module load mkl/2017--binary

mpirun --map-by ppr:1:socket ./aubry_NC.x > test_NC.dat
