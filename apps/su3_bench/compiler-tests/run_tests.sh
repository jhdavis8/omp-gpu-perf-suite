#!/bin/bash
mkdir -p data
if [ "$2" == "openmp" ] ; then
    for i in 0 1 2 3 ; do
	make -f Makefile.$2 clean
	make -f Makefile.$2 COMPILER=$1 VERSION=$i all
	# srun bench_f32_$2.exe &>> data/$1_test_f32.v$i.out
	# srun bench_f64_$2.exe &>> data/$1_test_f64.v$i.out
	# srun nvprof --metrics dram_read_transactions --metrics dram_write_transactions ./bench_f32_$2.exe 2>> data/$1_test_f32.v$i.out
	# srun nvprof --metrics dram_read_transactions --metrics dram_write_transactions ./bench_f64_$2.exe 2>> data/$1_test_f64.v$i.out
	srun nvprof --print-gpu-trace ./bench_f32_$2.exe 2>> data/$1_test_f32.v$i.out
	srun nvprof --print-gpu-trace ./bench_f64_$2.exe 2>> data/$1_test_f64.v$i.out
	srun nvprof --metrics achieved_occupancy ./bench_f32_$2.exe 2>> data/$1_test_f32.v$i.out
	srun nvprof --metrics achieved_occupancy ./bench_f64_$2.exe 2>> data/$1_test_f64.v$i.out
	echo "Results written to data/$1_test_f32.v$i.out and data/$1_test_f64.v$i.out"
    done
else
    make -f Makefile.$2 clean
    make -f Makefile.$2 all
    # srun bench_f32_$2.exe &>> data/$1_test_f32.out
    # srun bench_f64_$2.exe &>> data/$1_test_f64.out
    # srun nvprof --metrics dram_read_transactions --metrics dram_write_transactions ./bench_f32_$2.exe 2>> data/$1_test_f32.out
    # srun nvprof --metrics dram_read_transactions --metrics dram_write_transactions ./bench_f64_$2.exe 2>> data/$1_test_f64.out
    srun nvprof --print-gpu-trace ./bench_f32_$2.exe 2>> data/$1_test_f32.out
    srun nvprof --print-gpu-trace ./bench_f64_$2.exe 2>> data/$1_test_f64.out
    srun nvprof --metrics achieved_occupancy ./bench_f32_$2.exe 2>> data/$1_test_f32.out
    srun nvprof --metrics achieved_occupancy ./bench_f64_$2.exe 2>> data/$1_test_f64.out
    echo "Results written to data/$1_test_f32.out and data/$1_test_f64.out"
fi
