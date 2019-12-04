#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

float stencil ( float v1, float v2, float v3, float v4)
{
  return (v1 + v2 + v3 + v4) * 0.25f;
}

float max_error ( float prev_error, float old, float new )
{
  float t= fabsf( new - old );
  return t>prev_error? t: prev_error;
}

float laplace_step(float *in, float *out, int n)
{
  int i, j;
  float error=0.0f;
  for ( j=1; j < n-1; j++ )
    for ( i=1; i < n-1; i++ )
    {
      out[j*n+i]= stencil(in[j*n+i+1], in[j*n+i-1], in[(j-1)*n+i], in[(j+1)*n+i]);
      error = max_error( error, out[j*n+i], in[j*n+i] );
    }
  return error;
}


void laplace_init(float *in, int n, int nprocs, int me)
{
  int i;
  const float pi  = 2.0f * asinf(1.0f);
  int block_size  = (int)(n/nprocs);
  memset(in, 0, block_size*n*sizeof(float));
  
  int source = me * block_size;
  int destination = me * block_size + block_size;
  
  for (i=source; i<destination; i++) {
    float V = in[(i % block_size)*n] = sinf(pi*i / (n-1));
    in[ (i % block_size)*n+n-1 ] = V*expf(-pi);
  }
}

int main(int argc, char** argv)
{ 
  int n = 4096;
  int iter_max = 1000;
  
  const float tol = 1.0e-5f;
  float error= 1.0f;    
  
  int err, me, nprocs;
  
  

  // get runtime arguments 
  if (argc>1) {  n        = atoi(argv[1]); }
  if (argc>2) {  iter_max = atoi(argv[2]); }

  err=MPI_Init( &argc, &argv );
  //check succesful initialization of MPI
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &me);
  
  printf("size = %d\n", nprocs);
  printf("i am number = %d\n", me);
  
  float *A, *temp;
  
  int block_size 		= (int)(n/nprocs);
  int fr_index   		= 0;
  int lr_index   		= block_size - 1;
  int last_mes_index 	= block_size;
  int first_mes_index 	= block_size + 1;

  A    = (float*) malloc( (block_size+2)*n*sizeof(float) );
  temp = (float*) malloc( (block_size+2)*n*sizeof(float) );

 

  //  set boundary conditions
  laplace_init (A, n, nprocs, me);
  laplace_init (temp, n, nprocs, me);
  //  A[(n/128)*n+n/128] = 1.0f; // set singular point

  printf("Jacobi relaxation Calculation: %d x %d mesh,"
         " maximum of %d iterations\n", 
         n, n, iter_max );
         
  int iter = 0;
  while ( error > tol*tol && iter < iter_max )
  {
    iter++;
    error= laplace_step (A, temp, n);
    float *swap= A; A=temp; temp= swap; // swap pointers A & temp
  }
  error = sqrtf( error );
  printf("Total Iterations: %5d, ERROR: %0.6f, ", iter, error);
  printf("A[%d][%d]= %0.6f\n", n/128, n/128, A[(n/128)*n+n/128]);

  free(A); free(temp);
  
  err=MPI_Finalize();
  //check finalization
  return 0;
}
