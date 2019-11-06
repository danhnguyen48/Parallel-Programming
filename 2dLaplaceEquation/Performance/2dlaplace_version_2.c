#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

float stencil ( float v1, float v2, float v3, float v4 )
{
    return (v1 + v2 + v3 + v4) * 0.25;
}

// Left is mentioning where the current data is
void laplace_step ( float *M, int n, int m, bool left )
{
    int index;
    float tError = 0.0;
    if (left == true) index = 0;
    else index = 1;
    for ( int i=1; i < m-1; i++ )
        for ( int j=2; j < 2*n-2; j+=2 )
            {
                // calculate +1 position because we are calculating the new value of a point
                // v1 top, v2 right, v3 bottom, v4 left
                *(M + i*n*2 + j + 1 - index) = stencil(*(M + (i-1)*n*2 + j + index), 
                                                     *(M + i*n*2 + j + 2 + index), 
                                                     *(M + (i+1)*n*2 + j + index), 
                                                     *(M + i*n*2 + j - 2 + index));
                // Calculate the error 
                tError = fmaxf( tError, fabsf( *(M + i*n*2 + j) - *(M + i*n*2 + j + 1) ));
            }

    return sqrtf(tError);
}


void laplace_init ( float *M, int n, int m )
{
    const float pi  = 2.0f * asinf(1.0f);
    memset(M, 0, n*m*sizeof(float));
    for (int i=0; i<n*2; i++) {
        *(M + i) = 0.f; // top boundary
        *(M + n*(m-1) + i) = 0.f; // bottom boundary
    } 
    for (int j=0; j<m; j++) {
        *(M + j*n)  = sinf(pi*j / (n-1)); // left boundary
        *(M + j*n + n - 1) = sinf(pi*j / (n-1))*expf(-pi); // right boundary
    }
}

int main(int argc, char** argv)
{
    
    // m is number of rows
    // n is number of columns
    int n = 4096, m = 4096;
    int iter_max = 1000;
    float *A;
        
    const float tol = 1.0e-4f;
    float error= 1.0f;    

    // get runtime arguments: n, m and iter_max 
    if (argc>1) {  n        = atoi(argv[1]); }
    if (argc>2) {  m        = atoi(argv[2]); }
    if (argc>3) {  iter_max = atoi(argv[3]); }

    // remove A new for creating a bigger matrix contained A and Anew
    A    = (float*) malloc( n*m*sizeof(float) );

    //  set boundary conditions
    laplace_init (A, n, m);
    *(A + (n/128)*m+m/128) = 1.0f; // set singular point

    printf("Jacobi relaxation Calculation: %d x %d mesh,"
            " maximum of %d iterations\n", 
            n, m, iter_max );

    int iter = 0;
    while ( error > tol && iter < iter_max )
    {
        iter++;
        error = laplace_step (A, n, m, iter % 2 == 1);
        if (iter % (iter_max/10) == 0) printf("%5d, %0.6f\n", iter, error);
    }
    printf("Total Iterations: %5d, ERROR: %0.6f, \n", iter, error);
    printf("A[%d][%d]= %0.6f\n", n/128, m/128, *(A + (n/128)*m+m/128));

    free(A);
}