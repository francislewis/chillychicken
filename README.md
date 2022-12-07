# chillychicken
A parallel implementation of the Cooley-Turkey FFT algorithm in C for my University of Bristol masters unit in Advanced Computational Physics.

Two approaches are included, a shared memory approach using OpenMP and a distributed approach using MPI.

The effectiveness of these approaches are investigated using the [Blue Crystal 4](https://www.acrc.bris.ac.uk/acrc/phase4.htm) Supercomputer.

TODO:

- Check best compiler flags for serial and MPI 
- Create MAKE files for easier building
- Create header version for integration into other scripts?
