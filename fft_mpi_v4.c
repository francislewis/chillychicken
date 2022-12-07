// This one gets correct answers with random values, need to add some extra functionality though
// This is good for 1D, going to try add 2d functionality in v5
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<complex.h>
#include<time.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_complex_math.h>
#include<omp.h>
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
        printf("Size: %d, Dimension: %d, Number of Processes: %d\n", size, dimension, comm_size);
    }

    // Initialise arrays
    double complex sub_even[(size / comm_size / 2)];
    double complex sub_odd[(size / comm_size / 2)];
    double complex master_even[ (size / comm_size / 2) * comm_size];
    double complex master_odd[ (size / comm_size / 2) * comm_size];
    double complex results[size];

    // [0] - index, [1] - value, shared between master and sub array(s) to keep track when gathering results
    double complex sub_array[(size / comm_size)][2];
    double complex master_array[size][2];

    if(rank == 0){
        // Start timer on master rank
        double start = MPI_Wtime();

        // Seed random number with current time
        srand(time(NULL));

        // Generate random input values
        for(int i = 0; i < size/2; i++){
            master_array[i][0] = i;
            master_array[i][1] = (((int) rand())%1000)+(((int) rand())%1000)*I;
//            // Print to test
//            printf("inputdat1[%d] = %f+%fj\n",i,creal(master_array[i][1]),cimag(master_array[i][1]));
            master_array[i+size/2][0] =i;
            master_array[i+size/2][1] =0;
        }
    }

    // Size of chunks to be scattered and gathered
    int chunk_size = (size / comm_size) * 2;

    // Scatter the array to sub arrays for each process
    MPI_Scatter(master_array,chunk_size,MPI_DOUBLE_COMPLEX,sub_array,chunk_size,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD);

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

            results[k] = (sumeven + sumodd); // add even and odd parts
            results[k+size/2] = sumeven - sumodd; //Symmetry

//            // Print first 8 to check
//            if(k <= 8){
//                if(k == 0){
//                    printf(" \n\n TOTAL PROCESSED SAMPLES : %d\n",size);
//                }
//                printf("================================\n");
//                printf("XR[%d]: %.4f XI[%d]: %.4f \n",k,creal(results[k]),k,cimag(results[k]));
//                printf("================================\n");
//            }
        }
    }

	if(rank == 0){
	    // Print time taken
		double final = MPI_Wtime();
		printf("Time: %f\n",final);
	}

	MPI_Barrier(MPI_COMM_WORLD); // Pause until all processes have caught up
	MPI_Finalize(); // End MPI
	return 0;
}
