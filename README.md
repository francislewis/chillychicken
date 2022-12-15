# chillychicken
Parallel implementations of the Cooley-Turkey Fast Fourier Transform (FFT) algorithm in C for my University of Bristol masters unit in [Advanced Computational Physics](https://www.bris.ac.uk/unit-programme-catalogue/UnitDetails.jsa?unitCode=PHYSM0032), and as a way to teach myself C.

There is a [serial](fft_serial.c) approach along with three parallel approaches. A [shared memory approach](fft_omp.c) which uses OpenMP is included, along with a [distributed memory approach](fft_mpi.c) using MPI. A [combination of these approaches](fft_mpi_omp.c) is also included, although this is limited to  1D.   

The effectiveness of these approaches are investigated using the [Blue Crystal 4](https://www.acrc.bris.ac.uk/acrc/phase4.htm) Supercomputer.

Data on the speed of approach as a function of dimension, problem size and number of threads/proceses is included, alongside a [notebook for guided analysis](notebooks/Graph_Plotting.ipynb).

My full report can be found under [`reports/`](reports/)

### Requirements

- The [GNU Scientific Library (GSL)](https://www.gnu.org/software/gsl/) is required for the handling of 2D matrices.
- [FFTW](https://fftw.org/) is required to test the accuracy in the files under [`accuracy/`](accuracy/) 
- [OpenMP](https://www.openmp.org/) and [OpenMPI](https://www.open-mpi.org/) are also required for the parallel implementations.

### Compilation

Compilation instructions are included in [`compilation_instructions.txt`]([`compilation_instructions.txt). 
Instructions are provided for gcc and icc, although any compiler supporting the c99 standard will suffice.

