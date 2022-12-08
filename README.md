# chillychicken
A parallel implementation of the Cooley-Turkey FFT algorithm in C for my University of Bristol masters unit in [Advanced Computational Physics](https://www.bris.ac.uk/unit-programme-catalogue/UnitDetails.jsa?unitCode=PHYSM0032), and a way to teach myself C.

There is a [serial](fft_serial.c) approach along with two parallel approaches. A [shared memory approach](fft_omp.c) which uses OpenMP is included, along with a [distributed memory approach](fft_mpi.c) using MPI.   

The effectiveness of these approaches are investigated using the [Blue Crystal 4](https://www.acrc.bris.ac.uk/acrc/phase4.htm) Supercomputer.

### Requirements

- The [GNU Scientific Library (GSL)](https://www.gnu.org/software/gsl/) is required for the handling of 2D matrices.
- [FFTW](https://fftw.org/) is required to test the accuracy in the files under [`accuracy/`](accuracy/) 
- [OpenMP](https://www.openmp.org/) and [OpenMPI](https://www.open-mpi.org/) are also required for the parallel implementations.

TODO:

- Check best compiler flags for serial and MPI 
- Create MAKE files for easier building
- Create header version for integration into other scripts?
- Cleanup includes
- Combine MPI and OMP approaches
- Write test to compare results against known libraries
