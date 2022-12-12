#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<complex.h>
#include<time.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_complex_math.h>
#include<omp.h>

#define PI acos(-1)

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

    // Print info
    printf("Size: %d, Dimension: %d, Threads: %d\n", size, dimension, threads);

    // Start timer
    double start_time, end_time;
    start_time = omp_get_wtime();

    // Setup arrays
    double complex input[size];
    double complex results[size];

    // Seed random number with current time
    srand(time(NULL));

/*--------------------------------------------------------------------------------------------------------------------*/

    if(dimension==1){
        // Generate random input values
        for(int i = 0; i < size; i++){
            input[i] = (((int) rand())%1000)+(((int) rand())%1000)*I;

            /*        Print to check          */
            //printf("inputdat1[%d] = %f+%fj\n",i,creal(input[i]),cimag(input[i]));
        }

//        double setup_time = omp_get_wtime();
//        printf("Array Setup Time: %8.6f s\n", setup_time-start_time);

        // Main FFT loop
        #pragma omp parallel for
        for (int k = 0; k < size / 2; k++){
            double complex even = 0.0 + 0.0*I;
            double complex odd = 0.0 + 0.0*I;

            // Sum even and odd terms
            for(int i = 0; i < size/2; i++){
                even += (input[2*i]) * (cos(((2*PI)*((2*i)*k))/size) - (sin(((2*PI)*((2*i)*k))/size)*I));
                odd += (input[2*i+1]) * (cos(((2*PI)*((2*i+1)*k))/size) - (sin(((2*PI)*((2*i+1)*k))/size)*I));
            }

            // Store results
            results[k] = even + odd;
            results[k+size/2] = even - odd;

            // Print to check
//            printf("%d: %f + %fj \n",k,creal(results[k]),cimag(results[k]));
//            printf("%d: %f + %fj \n",k+size/2,creal(results[k+size/2]),cimag(results[k+size/2]));
        }
        end_time = omp_get_wtime();
        printf("Final Time: %f s\n", end_time-start_time);
    }
/*--------------------------------------------------------------------------------------------------------------------*/
	if(dimension==2){
	    // Allocate pointer for GSL matrix
        gsl_matrix_complex *rand_matrix = NULL;
        // Allocate memory for GSL matrix
        rand_matrix = gsl_matrix_complex_alloc(size, size);
        // gsl_matrix_complex_set_all(rand_matrix, GSL_COMPLEX_ONE);

        // Fill matrix with complex values
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
//        printf("\n");

//        double setup_time = omp_get_wtime();
//        printf("Matrix Setup Time: %8.6f s\n", setup_time-start_time);

        // Loop over each row of matrix and do FFT
        for (int row=0; row < size; row++){
            // Put data into input for FFT along this row
            for(int i = 0; i < size; i++){
                input[i] = GSL_REAL(gsl_matrix_complex_get(rand_matrix, row, i))+GSL_IMAG(gsl_matrix_complex_get(rand_matrix, row, i))*I;
	        }
	        // Perform FFT
            #pragma omp parallel for
	        for (int k = 0; k < size / 2; k++){
                double complex even = 0.0 + 0.0*I;
                double complex odd = 0.0 + 0.0*I;

                // Sum even and odd terms
                for(int i = 0; i < size/2; i++){
                    even += (input[2*i]) * (cos(((2*PI)*((2*i)*k))/size) - (sin(((2*PI)*((2*i)*k))/size)*I));
                    odd += (input[2*i+1]) * (cos(((2*PI)*((2*i+1)*k))/size) - (sin(((2*PI)*((2*i+1)*k))/size)*I));
                }

                results[k] = even + odd;
                results[k+size/2] = even - odd;

                gsl_matrix_complex_set(rand_matrix, row, k, gsl_complex_rect(creal(results[k]), cimag(results[k])));
                gsl_matrix_complex_set(rand_matrix, row, k+size/2, gsl_complex_rect(creal(results[k+size/2]), cimag(results[k+size/2])));
            }
        }

//        double first_fft = omp_get_wtime();
//        printf("First FFT Time: %8.6f s\n", first_fft-setup_time);

        // Transpose matrix (in place - thanks to GSL)
        #pragma omp critical
        gsl_matrix_complex_transpose(rand_matrix);

//        double first_transpose = omp_get_wtime();
//        printf("First Transpose  Time: %8.6f s\n", first_transpose-first_fft);

        // Loop over each col of matrix and do FFT (access via row due to transpose for speed)
        for (int col=0; col < size; col++){
            // Put data into input for FFT along this row
            for(int i = 0; i < size; i++){
                input[i] = GSL_REAL(gsl_matrix_complex_get(rand_matrix, col, i))+GSL_IMAG(gsl_matrix_complex_get(rand_matrix, col, i))*I;
	        }
            // Perform FFT
            #pragma omp parallel for
	        for (int k = 0; k < size / 2; k++){
                double complex even = 0.0 + 0.0*I;
                double complex odd = 0.0 + 0.0*I;

                // Sum even and odd terms
                for(int i = 0; i < size/2; i++){
                    even += (input[2*i]) * (cos(((2*PI)*((2*i)*k))/size) - (sin(((2*PI)*((2*i)*k))/size)*I));
                    odd += (input[2*i+1]) * (cos(((2*PI)*((2*i+1)*k))/size) - (sin(((2*PI)*((2*i+1)*k))/size)*I));
                }

                results[k] = even + odd;
                results[k+size/2] = even - odd;

                gsl_matrix_complex_set(rand_matrix, col, k, gsl_complex_rect(creal(results[k]), cimag(results[k])));
                gsl_matrix_complex_set(rand_matrix, col, k+size/2, gsl_complex_rect(creal(results[k+size/2]), cimag(results[k+size/2])));
            }
        }

//        double second_fft = omp_get_wtime();
//        printf("Second FFT Time: %8.6f s\n", second_fft-first_transpose);

        // Transpose matrix back for final results - stored in rand_matrix
        #pragma omp critical
        gsl_matrix_complex_transpose(rand_matrix);

//        // Print to check
//        for(int k = 0; k < size; k++){
//            for(int j = 0; j < size; j++){
//                printf("%d %d : %f + %f j\n", k, j, GSL_REAL(gsl_matrix_complex_get(rand_matrix, k, j)), GSL_IMAG(gsl_matrix_complex_get(rand_matrix, k, j)));
//            }
//        }
//        printf("\n");

//        double second_transpose = omp_get_wtime();
//        printf("Second Transpose Time: %8.6f s\n", second_transpose-second_fft);
        end_time = omp_get_wtime();
        printf("Total time: %8.6f s\n", end_time-start_time);

	}
	return 0;
}
