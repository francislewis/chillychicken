#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<complex.h>
#include<time.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_complex_math.h>
#include<omp.h>

#define PI acos(-1)

// Input pointer
double complex* input = NULL;

// Compute the even and odd parts of the FFT
void even_odd(int size, int k, double complex* Ek, double complex* Ok){
    if(k >= size/2){
        k -= size/2;
    }
    *Ek=0+0*I;
    *Ok=0+0*I;

    for(int m = 0; m <= size/2 - 1; m++)
    {
        // Update odd and even parts
        *Ek = *Ek + (input[2*m]*(cos(-2.0 * PI * m * k / (size/2)) + sin(-2.0 * PI * m * k / (size/2))*I));
        *Ok = *Ok + (input[2*m+1]*(cos(-2.0 * PI * m * k / (size/2)) + sin(-2.0 * PI * m * k / (size/2))*I));
    }
}

// Compute the twiddle factors in Cooley-Turkey FFT
double complex find_twiddle_term(int size, int k){
    double complex twiddle_term = cos(-2.0 * PI * k / size)+sin(-2.0 * PI * k / size)*I;
    return twiddle_term;
}

//Main function
int main(int argc, char* argv[]){
	// Ensure correct number of arguments
	if (argc < 4){
            printf("Usage: %s <SIZE> <DIMENSION> <THREADS>\n",argv[0]);
            return 1;
    }

	// Get command line arguments
	int size = atoi(argv[1]);
    int dimension = atoi(argv[2]);
    int threads = atoi(argv[3]);

    // Check if dimension in [1,2] using DeMorgan
    if (!(dimension>0 && dimension<3)){
        printf("Number of dimensions can only be 1 or 2\n");
        return 1;
    }

    // Set number of omp threads
    omp_set_num_threads(threads);

    printf("Size: %d, Dimension: %d, Threads: %d\n", size, dimension, threads);

    double start_time, end_time;
    start_time = omp_get_wtime();

    // Allocate memory for input values
    // input is the 1-D array passed through FFT
    // Once for 1-D and multiple times with different data for 2D
    input = malloc(size * sizeof(complex));

    double complex* temp_results = NULL;
    temp_results = malloc(size * sizeof(complex));

    double complex result, twiddle_term;
    double complex* Ek = malloc(sizeof(double complex));
    double complex* Ok = malloc(sizeof(double complex));

    // Seed random number with current time
    srand(time(NULL));

	if(dimension==1){
	    // Generate random input values
	    #pragma omp parallel for
	    for(int i = 0; i < size; i++){
            input[i] = (((int) rand())%1000)+(((int) rand())%1000)*I;
            //printf("inputdat1[%d] = %f+%fj\n",i,creal(input[i]),cimag(input[i]));
	    }
        #pragma omp parallel for
        for(int k = 0; k < size; k++){
            // Find the twiddle factor term
            twiddle_term = find_twiddle_term(size, k);

            // Compute FFT for even and odd values
            even_odd(size, k, Ek, Ok);

            // Compute result;
            result = *Ek + (twiddle_term* *Ok);

            // Store result
            temp_results[k] = result;
            //printf("%d: %f + %fj\n",k,creal(temp_results[k]),cimag(temp_results[k]));
        }
        end_time = omp_get_wtime();
        printf("Time: %8.6f s\n", end_time-start_time);
	}

    // dimension 2 using GSL
	if(dimension==2){

        // Allocate pointer for GSL matrix
        gsl_matrix_complex *rand_matrix = NULL;
        // Allocate memory for GSL matrix
        rand_matrix = gsl_matrix_complex_alloc(size, size);

        // Fill matrix with complex values
        #pragma omp parallel for collapse(2)
        for(int row = 0; row < size; row++){
            for(int col = 0; col < size; col++){
                int rand1 =(((int) rand())%1000);
                int rand2 =(((int) rand())%1000);
                // Create gsl complex number
                gsl_complex rand_fill = gsl_complex_rect((((int) rand())%1000), (((int) rand())%1000));
                // Add gsl complex to matrix
                gsl_matrix_complex_set(rand_matrix, row, col, rand_fill);
            }
        }

//        // Print to check
//        printf("start:\n");
//        for(int k = 0; k < size; k++){
//            for(int j = 0; j < size; j++){
//                printf("inputdat[%d][%d]= %f + %fj\n", k, j, GSL_REAL(gsl_matrix_complex_get(rand_matrix, k, j)), GSL_IMAG(gsl_matrix_complex_get(rand_matrix, k, j)));
//            }
//        }
        printf("\n");

        // Loop over each row of matrix and do FFT
        for (int row=0; row < size; row++){
            // Put data into input for FFT along this row
            for(int i = 0; i < size; i++){
                input[i] = GSL_REAL(gsl_matrix_complex_get(rand_matrix, row, i))+GSL_IMAG(gsl_matrix_complex_get(rand_matrix, row, i))*I;
	            // would be good to include this in loop below and do a loop collapse
	            // But getting wrong results when replacing i with col
	        }
            #pragma omp parallel for
            for (int col=0; col < size; col++){
                //input[col] = GSL_REAL(gsl_matrix_complex_get(rand_matrix, row, col))+GSL_IMAG(gsl_matrix_complex_get(rand_matrix, row, col))*I;
                // Do FFT
                twiddle_term = find_twiddle_term(size, col);
                even_odd(size, col, Ek, Ok);
                result = *Ek + (twiddle_term* *Ok);
                // Update matrix with result
                gsl_matrix_complex_set(rand_matrix, row, col, gsl_complex_rect(creal(result), cimag(result)));
            }
        }

        #pragma omp critical
        // Transpose matrix (in place - thanks to GSL)
        gsl_matrix_complex_transpose(rand_matrix);

        // Loop over each column (although accessing by row due to transpose)
        for (int col=0; col < size; col++){
            // Put data into input for FFT along this row
            for(int i = 0; i < size; i++){
                input[i] = GSL_REAL(gsl_matrix_complex_get(rand_matrix, col, i))+GSL_IMAG(gsl_matrix_complex_get(rand_matrix, col, i))*I;
	        }
	        #pragma omp parallel for
            for (int row=0; row < size; row++){
                // Do FFT
                twiddle_term = find_twiddle_term(size, row);
                even_odd(size, row, Ek, Ok);
                result = *Ek + (twiddle_term* *Ok);
                gsl_matrix_complex_set(rand_matrix, col, row, gsl_complex_rect(creal(result), cimag(result)));
            }
        }

        #pragma omp critical
        // Transpose matrix back for final results, stored in rand_matrix
        gsl_matrix_complex_transpose(rand_matrix);
        end_time = omp_get_wtime();
        printf("Time: %8.6f s\n", end_time-start_time);

//        // Print to check
//        for(int k = 0; k < size; k++){
//            for(int j = 0; j < size; j++){
//                printf("%d %d : %f + %f j\n", k, j, GSL_REAL(gsl_matrix_complex_get(rand_matrix, k, j)), GSL_IMAG(gsl_matrix_complex_get(rand_matrix, k, j)));
//            }
//        }
        printf("\n");
	}
    return 0;
}
