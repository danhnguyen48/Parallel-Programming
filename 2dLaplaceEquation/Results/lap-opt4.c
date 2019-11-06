#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

float stencil ( float v1, float v2, float v3, float v4 )
{
  return (v1 + v2 + v3 + v4) * .25;
}

float laplace_step_error ( float *in1, float *in2, int iter, int n, int m )
{
  int i, j;
  float error=0.0f;
  for ( i=1; i < m-1; i++ )
    for ( j=1; j < n-1; j++ ){
      if (iter % 2 == 0){
      in2[i*n+j]= stencil(in1[i*n+j+1], in1[i*n+j-1], in1[(i-1)*n+j], in1[(i+1)*n+j]);
      error = fmaxf( error, sqrtf( fabsf( in2[i*n+j] - in1[i*n+j] )));
    } else {
      in1[i*n+j]= stencil(in2[i*n+j+1], in2[i*n+j-1], in2[(i-1)*n+j], in2[(i+1)*n+j]);
      error = fmaxf( error, sqrtf( fabsf( in2[i*n+j] - in1[i*n+j] )));
    }
  }
  return error;

}


void laplace_init ( float *in1, float *in2, int n, int m )
{
  int i, j;
  const float pi  = 2.0f * asinf(1.0f);
  memset(in1, 0, n*m*sizeof(float));
  for (i=0; i<m; i++) {  in1[    i    ] = 0.f;                           in2[    i    ] = 0.f; }
  for (i=0; i<m; i++) {  in1[(n-1)*m+i] = 0.f;                           in2[(n-1)*m+i] = 0.f; }
  for (j=0; j<n; j++) {  in1[   j*m   ] = sinf(pi*j / (n-1));            in2[   j*m   ] = sinf(pi*j / (n-1)); }
  for (j=0; j<n; j++) {  in1[ j*m+m-1 ] = sinf(pi*j / (n-1))*expf(-pi);  in2[ j*m+m-1 ] = sinf(pi*j / (n-1))*expf(-pi); }
    
  
  
}

int main(int argc, char** argv)
{
  int n = 4096, m = 4096;
  int iter_max = 1000;
  float *A, *Anew;
    
  const float tol = 1.0e-4f;
  float error= 1.0f;    

  // get runtime arguments: n, m and iter_max 
  if (argc>1) {  n        = atoi(argv[1]); }
  if (argc>2) {  m        = atoi(argv[2]); }
  if (argc>3) {  iter_max = atoi(argv[3]); }

  A    = (float*) malloc( n*m*sizeof(float) );
  Anew = (float*) malloc( n*m*sizeof(float) );

  //  set boundary conditions
  laplace_init (A, Anew, n, m);
  A[(n/128)*m+m/128] = 1.0f; // set singular point
  Anew[(n/128)*m+m/128] = 1.0f; // set singular point

  printf("Jacobi relaxation Calculation: %d x %d mesh,"
         " maximum of %d iterations\n", 
         n, m, iter_max );

  int iter = 0;
  while ( error > tol && iter < iter_max ) 
  {
    error= laplace_step_error (A, Anew, iter, n, m);
    iter++;
    if (iter % (iter_max/10) == 0) printf("%5d, %0.6f\n", iter, error);
  }   
  printf("Total Iterations: %5d, ERROR: %0.6f, ", iter, error);
  if (iter % 2 == 0) printf("A[%d][%d]= %0.6f\n", n/128, m/128, A[(n/128)*m+m/128]);
  else printf("Anew[%d][%d]= %0.6f\n", n/128, m/128, Anew[(n/128)*m+m/128]);

  free(A); 
  free(Anew);
}