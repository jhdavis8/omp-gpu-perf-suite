#!/bin/bash
module purge
module load esslurm
module load cdt/19.11
module load PrgEnv-cray
module switch cce cce/9.1.3
module load craype-x86-skylake
module unload cray-libsci
module load cudatoolkit craype-accel-nvidia70
/bin/bash compiler-tests/run_tests.sh cray openmp
