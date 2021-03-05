#!/bin/bash
module purge
module load esslurm
module load pgi/20.4
module load cuda
mkdir -p data
/bin/bash compiler-tests/run_tests.sh nvcc cuda
