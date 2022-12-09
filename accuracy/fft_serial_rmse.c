#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<fftw3.h>
#include<complex.h>
#include<time.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_complex_math.h>

#define PI acos(-1)

int main(int argc, char* argv[]){

    // Ensure correct number of arguments
	if (argc > 3){
            printf("Usage: %s <SIZE> <DIMENSION> \n",argv[0]);
            return 1;
    }

	// Get command line arguments
	int size = atoi(argv[1]);
    int dimension = atoi(argv[2]);

    // Check if dimension in [1,2] using DeMorgan
    if (!(dimension>0 && dimension<3)){
        printf("Number of dimensions can only be 1 or 2\n");
        return 1;
    }

    // Print info
    printf("Size: %d, Dimension: %d\n", size, dimension);

    // Setup arrays
    double complex input[size];
    double complex results[size];

    // Seed random number with current time
    srand(time(NULL));

/*--------------------------------------------------------------------------------------------------------------------*/

    if(dimension==1){
        // Create input and output for FFTW
        fftw_complex in[size], out[size];
        // Create FFTW plan
        fftw_plan p = fftw_plan_dft_1d(size, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

        // Generate random input values
        for(int i = 0; i < size; i++){
            input[i] = (((int) rand())%1000)+(((int) rand())%1000)*I;
            in[i][0] = creal(input[i]);
            in[i][1] = cimag(input[i]);
        }

        // FFTW
        fftw_execute(p);

        // Main FFT loop
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

        }

        double error = 0;
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
/*--------------------------------------------------------------------------------------------------------------------*/
	if(dimension==2){
	    // Create input and output arrays for FFTW
        fftw_complex *in;
        fftw_complex *out;
        in = (fftw_complex*) fftw_malloc(size*size * sizeof(fftw_complex));
        out = (fftw_complex*) fftw_malloc(size*size * sizeof(fftw_complex));

        // Create FFTW plan
	    fftw_plan p = fftw_plan_dft_2d(size, size, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

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

                // Fill in FFTW input array
                in[col+size*row][0] = GSL_REAL(gsl_matrix_complex_get(rand_matrix, row, col));
                in[col+size*row][1] = GSL_IMAG(gsl_matrix_complex_get(rand_matrix, row, col));
            }
        }

        // Carry out FFTW after input values filled in
        fftw_execute(p);

        // Loop over each row of matrix and do FFT
        for (int row=0; row < size; row++){
            // Put data into input for FFT along this row
            for(int i = 0; i < size; i++){
                input[i] = GSL_REAL(gsl_matrix_complex_get(rand_matrix, row, i))+GSL_IMAG(gsl_matrix_complex_get(rand_matrix, row, i))*I;
	        }
	        // Perform FFT
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

        // Transpose matrix (in place - thanks to GSL)
        gsl_matrix_complex_transpose(rand_matrix);

        // Loop over each col of matrix and do FFT (access via row due to transpose for speed)
        for (int col=0; col < size; col++){
            // Put data into input for FFT along this row
            for(int i = 0; i < size; i++){
                input[i] = GSL_REAL(gsl_matrix_complex_get(rand_matrix, col, i))+GSL_IMAG(gsl_matrix_complex_get(rand_matrix, col, i))*I;
	        }

            // Perform FFT
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

        // Transpose matrix back for final results - stored in rand_matrix
        gsl_matrix_complex_transpose(rand_matrix);

        double error = 0;
        for (int row = 0; row < size; row++){
            for (int col =0; col< size; col++){
                // Print to check
//                printf("FFTW: %3d %3d: %+9.5f %+9.5f I\n", row,col, out[col+size*row][0], out[col+size*row][1]);
//                printf("Mine: %3d %3d: %+9.5f %+9.5f I\n", row, col, GSL_REAL(gsl_matrix_complex_get(rand_matrix, row, col)), GSL_IMAG(gsl_matrix_complex_get(rand_matrix, row, col)));

                // Calculate RMS error SUM[(observed_i-expected_i)*(observed_i-expected_i)/size] for real and complex parts
                error += (GSL_REAL(gsl_matrix_complex_get(rand_matrix, row, col)) - (out[col+size*row][0]))*(GSL_REAL(gsl_matrix_complex_get(rand_matrix, row, col)) - (out[col+size*row][0]))/size;
                error += (GSL_IMAG(gsl_matrix_complex_get(rand_matrix, row, col)) - (out[col+size*row][1]))*(GSL_IMAG(gsl_matrix_complex_get(rand_matrix, row, col)) - (out[col+size*row][1]))/size;
            }
        }
        printf("\nRMSE: %f\n", error);
        fftw_destroy_plan(p);

	}
	return 0;
}
