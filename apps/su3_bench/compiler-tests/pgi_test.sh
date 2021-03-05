#!/bin/bash
module purge
module load esslurm
module load pgi/20.4
module load cuda
/bin/bash compiler-tests/run_tests.sh pgi openacc
