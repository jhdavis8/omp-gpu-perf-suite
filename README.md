# OpenMP GPU Offload Performance Mini-app Suite

This repository stores source codes for the benchmark apps used in "Performance Assessment of OpenMP Compilers Targeting NVIDIA V100 GPUs" (WACCPD 2020). For details on how each app is used for performance evaluation and what modifications (if any) were made in the process, see the original paper.

[Link to the paper on arXiv.](https://arxiv.org/abs/2010.09454)

Apps included:
* su3 (MILC)
* BabelStream
* laplace **(Pending)**
* GPP (BerkeleyGW) **(Pending)**
* ToyPush (XGC)

## Licensing/Availability

Licensing information for each app is included in its subfolder under `/apps`. The laplace and GPP apps have not yet been added as some additional work is needed to ensure the codes can be provided with the appropriate licensing. The upstream version of the GPP app is available [here](https://gitlab.com/NERSC/nersc-proxies/BerkeleyGW-Kernels-CPP) from NERSC.

## Building and Running

At the moment an integrated build system across all of the apps has not been created. To build and run any of the apps, see the Makefile and README information provided in its subfolder.
