#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include<complex.h>

int main(int argc, char* argv[]){
    int N = atoi(argv[1]);
    int dimension = atoi(argv[2]);
    srand(time(NULL));
    float complex data[N][dimension];
    int i;
    int d;
    for (i = 0; i < N; i++) {
        for (d = 0; d < dimension; d++){
            data[i][d] = ((float) rand())+((float) rand())*I;
        }
    }
    for(int k = 0; k < N; k++){
        for(int j = 0; j < dimension; j++){
            printf("%d %d : real: %f + imag: %f\n", k, j, creal(data[k][j]), cimag(data[k][j]));
        }
    }
}

