chillychicken
=============

Parallel implementations of the Cooley-Turkey Fast Fourier Transform
(FFT) algorithm in C for my University of Bristol masters unit in
Advanced Computational Physics, and as a way to teach myself C.

There is a serial approach along with three parallel
approaches. A shared memory approach which uses OpenMP is
included, along with a distributed memory approach using
MPI. A combination of these approaches is also
included, although this is limited to 1D.

The effectiveness of these approaches are investigated using the Blue
Crystal Phase 4 Supercomputer.

Data on the speed of approach as a function of dimension, problem size
and number of threads/proceses is included, alongside a notebook for
guided analysis.

My full report can be found under reports/

##### Requirements ##### 

-   The GNU Scientific Library(GSL) is required for the handling of 2D matrices.
-   FFTW is required to test the accuracy in the files under accuracy/
-   OpenMP and OpenMPI are also required for the parallel implementations.

##### Compilation Instructions #####

Compilation instructions are included in compilation_instructions.txt

Instructions are provided for gcc and icc, although any compiler 
supporting the c99 standard will suffice.
