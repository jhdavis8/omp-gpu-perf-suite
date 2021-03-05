## su3_bench: Lattice QCD SU(3) matrix-matrix multiply microbenchmark  
*This is a fork of Doug Doerfler's su3_bench repository, created to document the commands and modules used to collect data for comparing compiler OpenMP offload performance using the benchmark. These scripts can be found in the compiler-tests folder.*

The purpose of this microbenchmark is to provide a means to explore different programming methodologies using a simple, but not trivial, mathematical kernel. The kernel is based on the mult_su3_nn() SU(3) matrix-matrix multiply routine in the MILC Lattice Quantum Chromodynamics(LQCD) code. Matrix-matrix (and matrix-vector) SU(3) operations are a fundamental building block of LQCD applications. Most LQCD applications use custom implementations of these kernels, and they are usually written in machine specific languages and/or  intrinsics. 

### Design
The code is written in standard C and C++. The main driver routine is used for all programming model implementations, with programming model specific implementations self contained in respective C++ include files. Programming methods implemented include: OpenCL, OpenMP, OpenACC, SYCL, CUDA and HIP.

The code is the documentation. It's simple enough, so dive in.

Various makefiles are also included, one for each of the respective compile environments I've tried so far.

### Usage
*To test a particular compiler, run the appropriate script in the compiler-tests folder. For example, clang_test.sh will run the OpenMP benchmark in 32- and 64-bit precision, for all code versions. The results will be output to a data folder in the repository root. Each compiler's script loads all the appropriate modules, and they are designed to be run from within an interactive allocation on Cori.*

*If you build a particular code verison, you can see the command line arguments using the following:*
*bench_xxx.exe -h*

### Contact info
Doug Doerfler
dwdoerf@lbl.gov

Josh Davis (for evaluation scripts, etc.)
jhdavis@udel.edu