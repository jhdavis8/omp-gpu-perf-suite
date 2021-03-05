#!/bin/bash
module purge
module load esslurm
module load PrgEnv-llvm/11.0.0-git_20200409
make -f Makefile.openmp clean; make -f Makefile.openmp COMPILER=clang VERSION=1 all
for i in 32 64 96 128 160; do srun bench_f32_openmp.exe -t $i -n 1600; done
for i in 400 800 1200 1600 2000; do srun bench_f32_openmp.exe -t 64 -n $i; done
