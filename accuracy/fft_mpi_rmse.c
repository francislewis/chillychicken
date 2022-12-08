// Aiming to combine 1D and 2D approaches here
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include <fftw3.h>
#include<complex.h>
#include<time.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_complex_math.h>
#include <mpi.h>

#define PI acos(-1)

int main(int argc, char *argv[]){
	int rank,comm_size;
	MPI_Init(&argc, &argv); // Start MPI
	MPI_Comm_size(MPI_COMM_WORLD,&comm_size); // Number of MPI Processes
	MPI_Comm_rank(MPI_COMM_WORLD,&rank); // Current rank

    if (rank==0){
        // Ensure correct number of arguments
        if (argc > 3){
            printf("Usage: mpiexec -np <PROCESSES> %s <SIZE> <DIMENSION>\n",argv[0]);
            return 1;
        }
        // Check if dimension in [1,2] using DeMorgan
        if (!(atoi(argv[2])>0 && atoi(argv[2])<3)){
            printf("Number of dimensions can only be 1 or 2\n");
            return 1;
        }
    }

    // Get command line arguments
    int size = atoi(argv[1]);
    int dimension = atoi(argv[2]);

    // Print info
    if (rank==0){
        printf("\nSize: %d, Dimension: %d, Number of Processes: %d", size, dimension, comm_size);
    }

    // Initialise arrays
    double complex sub_even[(size / comm_size / 2)];
    double complex sub_odd[(size / comm_size / 2)];
    double complex master_even[ (size / comm_size / 2) * comm_size];
    double complex master_odd[ (size / comm_size / 2) * comm_size];
    double complex results[size];

    // [0] - index, [1] - value, shared between master and sub array(s) to keep track when gathering results
    double complex sub_array[(size / comm_size)][2];
    double complex input[size][2];

    // Create input and output for FFTW
    fftw_complex in[size], out[size];

    // Create FFTW plan
    fftw_plan p = fftw_plan_dft_1d(size, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

    // Create input and output arrays for 2D FFTW
    fftw_complex *in2;
    fftw_complex *out2;

    in2 = (fftw_complex*) fftw_malloc(size*size * sizeof(fftw_complex));
    out2 = (fftw_complex*) fftw_malloc(size*size * sizeof(fftw_complex));

    // Create FFTW plan
    fftw_plan q = fftw_plan_dft_2d(size, size, in2, out2, FFTW_FORWARD, FFTW_ESTIMATE);

    if (dimension==1){
        if(rank == 0){
            // Seed random number with current time
            srand(time(NULL));

            // Generate random input values
            for(int i = 0; i < size/2; i++){
                input[i][0] = i;
                input[i][1] = (((int) rand())%1000)+(((int) rand())%1000)*I;
                input[i+size/2][0] = i+size/2;
                input[i+size/2][1] = 0+0*I;

                // Fill in FFTW input
                in[i][0] = creal(input[i][1]);
                in[i][1] = cimag(input[i][1]);
                in[i+size/2][0] = 0;
                in[i+size/2][1] = 0;
            }
            // Run FFTW
            fftw_execute(p);
        }

        // Size of chunks to be scattered and gathered
        int chunk_size = (size / comm_size) * 2;

        // Scatter the array to sub arrays for each process
        MPI_Scatter(input,chunk_size,MPI_DOUBLE_COMPLEX,sub_array,chunk_size,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD);

        // Main loop, only have to loop over half due to symmetry
        for (int k = 0; k < size / 2; k++){
            double complex sumeven = 0.0 + 0.0*I;
            double complex sumodd = 0.0 + 0.0*I;

            // Calculate even and odd parts
            for(int i = 0; i < (size/comm_size)/2; i++){
                double factoreven , factorodd = 0.0;

                // Shift index numbers within sub arrays
                int evenshift = rank * creal(sub_array[2*i][0]);
                int oddshift = rank * creal(sub_array[2*i + 1][0]);

                double complex even = sub_array[2*i][1];
                double complex odd = sub_array[2*i + 1][1];

                // If master rank, don't shift
                // else, shift index of sub array in order for results to be in correct place
                if(rank == 0){
                    sub_even[i] = (even * (cos((((2*PI)*((2*i)*k))/size)) - (sin((((2*PI)*((2*i)*k))/size))*I)));
                    sub_odd[i] = (odd * (cos((((2*PI)*((2*i+1)*k))/size)) - (sin((((2*PI)*((2*i+1)*k))/size))*I)));
                }
                else{
                    sub_even[i] = (even * (cos((((2*PI)*((evenshift)*k))/size)) - (sin((((2*PI)*((evenshift)*k))/size))*I)));
                    sub_odd[i] = (odd * (cos((((2*PI)*((oddshift)*k))/size)) - (sin((((2*PI)*((oddshift)*k))/size))*I)));
                }
            }

            // Master rank gathers even and odd parts from sub processes to create master even and odd arrays
            MPI_Gather(sub_even,(size / comm_size / 2),MPI_DOUBLE_COMPLEX,master_even,(size / comm_size / 2), MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD);
            MPI_Gather(sub_odd,(size / comm_size / 2),MPI_DOUBLE_COMPLEX,master_odd,(size / comm_size / 2), MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD);

            if(rank == 0){
                // Master rank sums even and odd parts
                for(int i = 0; i < (size / comm_size / 2) * comm_size; i++){
                    sumeven += master_even[i];
                    sumodd += master_odd[i];
                }

                // Add even and odd parts
                results[k] = (sumeven + sumodd);
                // Take advantage of symmetry
                results[k+size/2] = sumeven - sumodd; // Symmetry
            }
        }

        if(rank == 0){
            double error =0;

            for (int i = 0; i < size; i++){
            // Print to check
//            printf("FFTW: %3d %+9.5f %+9.5f I\n", i, out[i][0], out[i][1]);
//            printf("Mine: %3d %+9.5f %+9.5f I\n", i, creal(results[i]), cimag(results[i]));

            // Calculate RMS error SUM[(observed_i-expected_i)*(observed_i-expected_i)/size] for real and complex parts
            error += ((creal(results[i]) - out[i][0])*(creal(results[i]) - out[i][0])/size);
            error += ((cimag(results[i]) - out[i][1])*(cimag(results[i]) - out[i][1])/size);
        }

        printf("\nRMSE: %f\n", error);
        fftw_destroy_plan(p);
        }

        MPI_Barrier(MPI_COMM_WORLD); // Pause until all processes have caught up
        MPI_Finalize(); // End MPI
        return 0;
    }

/*###################################################################################################################*/

    if (dimension==2){
        // Allocate pointer for GSL matrix
        gsl_matrix_complex *rand_matrix = NULL;
        // Allocate memory for GSL matrix
        rand_matrix = gsl_matrix_complex_alloc(size, size);

        if(rank == 0){
            // Seed random number with current time
            srand(time(NULL));

            // Fill matrix with random values
            for(int row = 0; row < size; row++){
                for(int col = 0; col < size; col++){
                    // Create gsl complex number
                    gsl_complex rand_fill = gsl_complex_rect((((int) rand())%1000), (((int) rand())%1000));
                    // Add gsl complex to matrix
                    gsl_matrix_complex_set(rand_matrix, row, col, rand_fill);

                    // Fill in FFTW input array
                    in2[col+size*row][0] = GSL_REAL(gsl_matrix_complex_get(rand_matrix, row, col));
                    in2[col+size*row][1] = GSL_IMAG(gsl_matrix_complex_get(rand_matrix, row, col));

                    // Print to check
//                    printf("Mine: %d%d: %f+%fj\n", row, col, GSL_REAL(gsl_matrix_complex_get(rand_matrix, row, col)),GSL_IMAG(gsl_matrix_complex_get(rand_matrix, row, col)));
//                    printf("FFTW: %d%d: %f+%fj\n", row, col, in2[col+size*row][0],in2[col+size*row][1]);
                }
            }

            // Carry out FFTW after input values filled in
            fftw_execute(q);

        }

        // Loop over each row of matrix and do FFT
        for (int row=0; row < size; row++){
            // Put data into input for FFT along this row
            for(int i = 0; i < size; i++){
                input[i][1] = GSL_REAL(gsl_matrix_complex_get(rand_matrix, row, i))+GSL_IMAG(gsl_matrix_complex_get(rand_matrix, row, i))*I;
                input[i][0] = i;
            }

            // Size of chunks to be scattered and gathered
            int chunk_size = (size / comm_size) * 2;

            // Scatter the array to sub arrays for each process
            MPI_Scatter(input,chunk_size,MPI_DOUBLE_COMPLEX,sub_array,chunk_size,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD);

            // Main loop, only have to loop over half due to symmetry
            for (int k = 0; k < size / 2; k++){
                double complex sumeven = 0.0 + 0.0*I;
                double complex sumodd = 0.0 + 0.0*I;

                // Calculate even and odd parts
                for(int i = 0; i < (size/comm_size)/2; i++){
                    double factoreven , factorodd = 0.0;

                    // Shift index numbers within sub arrays
                    int evenshift = rank * creal(sub_array[2*i][0]);
                    int oddshift = rank * creal(sub_array[2*i + 1][0]);

                    double complex even = sub_array[2*i][1];
                    double complex odd = sub_array[2*i + 1][1];

                    // If master rank, don't shift
                    // else, shift index of sub array in order for results to be in correct place
                    if(rank == 0){
                        sub_even[i] = (even * (cos((((2*PI)*((2*i)*k))/size)) - (sin((((2*PI)*((2*i)*k))/size))*I)));
                        sub_odd[i] = (odd * (cos((((2*PI)*((2*i+1)*k))/size)) - (sin((((2*PI)*((2*i+1)*k))/size))*I)));
                    }
                    else{
                        sub_even[i] = (even * (cos((((2*PI)*((evenshift)*k))/size)) - (sin((((2*PI)*((evenshift)*k))/size))*I)));
                        sub_odd[i] = (odd * (cos((((2*PI)*((oddshift)*k))/size)) - (sin((((2*PI)*((oddshift)*k))/size))*I)));
                    }
                }

                // Master rank gathers even and odd parts from sub processes to create master even and odd arrays
                MPI_Gather(sub_even,(size / comm_size / 2),MPI_DOUBLE_COMPLEX,master_even,(size / comm_size / 2), MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD);
                MPI_Gather(sub_odd,(size / comm_size / 2),MPI_DOUBLE_COMPLEX,master_odd,(size / comm_size / 2), MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD);

                if(rank == 0){
                    // Master rank sums even and odd parts
                    for(int i = 0; i < (size / comm_size / 2) * comm_size; i++){
                        sumeven += master_even[i];
                        sumodd += master_odd[i];
                    }

                    results[k] = (sumeven + sumodd); // Add even and odd parts
                    results[k+size/2] = sumeven - sumodd; // Symmetry

                    gsl_matrix_complex_set(rand_matrix, row, k, gsl_complex_rect(creal(results[k]), cimag(results[k])));
                    gsl_matrix_complex_set(rand_matrix, row, k+size/2, gsl_complex_rect(creal(results[k+size/2]), cimag(results[k+size/2])));
                }
            }
        }

        // Pause until all processes have caught up
        MPI_Barrier(MPI_COMM_WORLD);

        // Transpose matrix
        gsl_matrix_complex_transpose(rand_matrix);

        // Loop over each col of matrix and do FFT (access by row due to transpose)
        for (int col=0; col < size; col++){
            // Put data into input for FFT along this row
            for(int i = 0; i < size; i++){
                input[i][1] = GSL_REAL(gsl_matrix_complex_get(rand_matrix, col, i))+GSL_IMAG(gsl_matrix_complex_get(rand_matrix, col, i))*I;
                input[i][0] = i;
            }

            // Size of chunks to be scattered and gathered
            int chunk_size = (size / comm_size) * 2;

            // Scatter the array to sub arrays for each process
            MPI_Scatter(input,chunk_size,MPI_DOUBLE_COMPLEX,sub_array,chunk_size,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD);

            // Main loop, only have to loop over half due to symmetry
            for (int k = 0; k < size / 2; k++){
                double complex sumeven = 0.0 + 0.0*I;
                double complex sumodd = 0.0 + 0.0*I;

                // Calculate even and odd parts
                for(int i = 0; i < (size/comm_size)/2; i++){
                    double factoreven , factorodd = 0.0;

                    // Shift index numbers within sub arrays
                    int evenshift = rank * creal(sub_array[2*i][0]);
                    int oddshift = rank * creal(sub_array[2*i + 1][0]);

                    double complex even = sub_array[2*i][1];
                    double complex odd = sub_array[2*i + 1][1];

                    // If master rank, don't shift
                    // else, shift index of sub array in order for results to be in correct place
                    if(rank == 0){
                        sub_even[i] = (even * (cos((((2*PI)*((2*i)*k))/size)) - (sin((((2*PI)*((2*i)*k))/size))*I)));
                        sub_odd[i] = (odd * (cos((((2*PI)*((2*i+1)*k))/size)) - (sin((((2*PI)*((2*i+1)*k))/size))*I)));
                    }
                    else{
                        sub_even[i] = (even * (cos((((2*PI)*((evenshift)*k))/size)) - (sin((((2*PI)*((evenshift)*k))/size))*I)));
                        sub_odd[i] = (odd * (cos((((2*PI)*((oddshift)*k))/size)) - (sin((((2*PI)*((oddshift)*k))/size))*I)));
                    }
                }

                // Master rank gathers even and odd parts from sub processes to create master even and odd arrays
                MPI_Gather(sub_even,(size / comm_size / 2),MPI_DOUBLE_COMPLEX,master_even,(size / comm_size / 2), MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD);
                MPI_Gather(sub_odd,(size / comm_size / 2),MPI_DOUBLE_COMPLEX,master_odd,(size / comm_size / 2), MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD);

                if(rank == 0){
                    // Master rank sums even and odd parts
                    for(int i = 0; i < (size / comm_size / 2) * comm_size; i++){
                        sumeven += master_even[i];
                        sumodd += master_odd[i];
                    }

                    results[k] = (sumeven + sumodd); // Add even and odd parts
                    results[k+size/2] = sumeven - sumodd; // Symmetry

                    gsl_matrix_complex_set(rand_matrix, col, k, gsl_complex_rect(creal(results[k]), cimag(results[k])));
                    gsl_matrix_complex_set(rand_matrix, col, k+size/2, gsl_complex_rect(creal(results[k+size/2]), cimag(results[k+size/2])));
                }
            }
        }

        // Pause until all processes have caught up
        MPI_Barrier(MPI_COMM_WORLD);

        // Transpose matrix back for final results, stored in rand_matrix
        gsl_matrix_complex_transpose(rand_matrix);

        if(rank == 0){
            double error = 0;
            // Calculate RMS error SUM[(observed_i-expected_i)*(observed_i-expected_i)/size] for real and complex parts
            for (int row = 0; row < size; row++){
                for (int col =0; col< size; col++){
                    // Print to check
//                    printf("FFTW: %3d %3d: %+9.5f %+9.5f I\n", row,col, out2[col+size*row][0], out2[col+size*row][1]);
//                    printf("Mine: %3d %3d: %+9.5f %+9.5f I\n", row, col, GSL_REAL(gsl_matrix_complex_get(rand_matrix, row, col)), GSL_IMAG(gsl_matrix_complex_get(rand_matrix, row, col)));
                    error += (GSL_REAL(gsl_matrix_complex_get(rand_matrix, row, col)) - (out2[col+size*row][0]))*(GSL_REAL(gsl_matrix_complex_get(rand_matrix, row, col)) - (out2[col+size*row][0]))/size;
                    error += (GSL_IMAG(gsl_matrix_complex_get(rand_matrix, row, col)) - (out2[col+size*row][1]))*(GSL_IMAG(gsl_matrix_complex_get(rand_matrix, row, col)) - (out2[col+size*row][1]))/size;
                }
            }
            printf("\nRMSE: %f\n", error);
            fftw_destroy_plan(p);
        }
        // End MPI
        MPI_Finalize();
        return 0;
    }
}
