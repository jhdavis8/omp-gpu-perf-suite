# OpenMP GPU Offload Performance Mini-app Suite

This repository stores source codes for the benchmark apps used in "Performance Assessment of OpenMP Compilers Targeting NVIDIA V100 GPUs" (WACCPD 2020). For details on how each app is used for performance evaluation and what modifications (if any) were made in the process, see the original paper.

[Link to the paper on arXiv.](https://arxiv.org/abs/2010.09454)

Apps included:
* SU(3) (MILC)
* BabelStream
* laplace **(Pending)**
* GPP (BerkeleyGW) **(Pending)**
* ToyPush (XGC)

## Upstream sources:

### SU(3)
* Repo: `https://gitlab.com/NERSC/nersc-proxies/su3_bench`
* Branch: `master`
* Commit: `74ce86e5`

### BabelStream
* Repo: `https://github.com/UoB-HPC/BabelStream`
* Branch: `main`
* Commit: `9025afe`

### laplace
* N/A, see `Chandrasekaran, S., Juckeland, G.: OpenACC for Programmers: Concepts and Strategies. Addison-Wesley Professional (2017)`

### GPP
* Repo: `https://github.com/UoB-HPC/BabelStream`
* Branch: `master`
* Commit: `95c15073`

### ToyPush
* Repo: `https://gitlab.com/NERSC/toypush`
* Branch: `oacc-omp`
* Commit: `1870a860`

## Licensing/Availability

Licensing information for each app is included in its subfolder under `/apps`. The laplace and GPP apps have not yet been added as some additional work is needed to ensure the codes can be provided with the appropriate licensing.

The upstream version of the GPP app (without the modifications described in the paper) is available [here](https://gitlab.com/NERSC/nersc-proxies/BerkeleyGW-Kernels-CPP) from NERSC.

## Building and Running

At the moment an integrated build system across all of the apps has not been created. To build and run any of the apps, see the Makefile and README information provided in the relevant subfolder.
