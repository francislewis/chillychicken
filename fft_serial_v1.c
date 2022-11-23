#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<complex.h>

#define PI acos(-1)

// Input pointer
complex* xinput = NULL;

// Compute the even and odd parts of the FFT
// Returns complex pointer to results
void even_odd(int size, int k, float complex* Ek, float complex* Ok)
{
    if(k >= size/2){
        k -= size/2;
    }
    *Ek=0+0*I;
    *Ok=0+0*I;

    float complex exponent;
    double exp;

    for(int m = 0; m <= size/2 - 1; m++)
    {
        // Find exponent term
        exponent = cos(-2.0 * PI * m * k / (size/2)) + sin(-2.0 * PI * m * k / (size/2))*I;

        // Update odd and even parts
        *Ek = *Ek + (xinput[2*m]*exponent);
        *Ok = *Ok + (xinput[2*m+1]*exponent);
    }

}

// Compute the twiddle factors in Cooley-Turkey FFT
float complex find_twiddle_term(int size, int k)
{
    double exp = -2.0 * PI * k / size;
    float complex twiddle_term = cos(exp)+sin(exp)*I;

    return twiddle_term;
}

//Main function
int main(int argc, char* argv[])
{
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
        printf("Number of dimensions can only be 1 or 2");
        return 1;
    }

	// Allocate memory for N complex inputs
	xinput = malloc(size * sizeof(complex));

	// Set some defined values to test
	xinput[0]=2.3+4.5*I;
	xinput[1]=2.9+6.9*I;
	xinput[2]=2+3.5*I;
	xinput[3]=8.8+2.4*I;
	xinput[4]=8.4+0.4*I;
	xinput[5]=5.3+5.3*I;
	xinput[6]=4.9+2.3*I;
	xinput[7]=8.7+8.1*I;

	// Zero padding
	for (int i = 8; i < size; i++)
	{
		xinput[i]=0+0*I;
	}

    float complex result;
    float complex twiddle_term;

    // Allocate memory
    float complex* Ek = malloc(sizeof(complex));
    float complex* Ok = malloc(sizeof(complex));

	printf("Size: %d, Dimension: %d, Threads: %d\n", size, dimension, threads);

    for(int k = 0; k < 8; k++)
    {
		printf("======================\n");

        // Find the twiddle factor term
        twiddle_term = find_twiddle_term(size, k);

        // Compute FFT for even and odd values
        even_odd(size, k, Ek, Ok);

        // Compute result;
        result = *Ek + (twiddle_term* *Ok);

		printf("Real Value at [%d]: %f\n", k, creal(result));
		printf("Complex Value at [%d]: %f\n", k, cimag(result));
    }

	printf("======================\n");

    return 0;
}
