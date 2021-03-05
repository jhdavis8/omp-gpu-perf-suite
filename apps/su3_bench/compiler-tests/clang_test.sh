#!/bin/bash
module purge
module load esslurm
module load PrgEnv-llvm/11.0.0-git_20200409
/bin/bash compiler-tests/run_tests.sh clang openmp
