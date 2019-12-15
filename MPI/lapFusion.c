#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <stdbool.h>

float stencil ( float v1, float v2, float v3, float v4)
{
  return (v1 + v2 + v3 + v4) * 0.25f;
}

float max_error ( float prev_error, float old, float new )
{
  float t= fabsf( new - old );
  return t>prev_error? t: prev_error;
}

float laplace_step(float *in, float *out, int n, int nprocs, int me)
{
  int i, j;
  float error=0.0f;

  int block_size  = (int)(n/nprocs);  
  int source = 0;
  int destination = block_size;

  // Check if me == fisrt matrix, source = 1
  // me == last matrix, destination = block_size - 1;
  if (me == 0) {
    source = 1;
  } else if (me == nprocs - 1) {
    destination = block_size - 1;
  }

  for ( i=source; i < destination; i++ )
    for ( j=1; j < n-1; j++ )
    {
      float aboveE = i == 0 ? in[block_size*n + j] : in[(i-1)*n + j];
      float belowE = i == block_size - 1 ? in[(block_size+1)*n + j] : in[(i+1)*n + j];
      out[i*n+j]= stencil(in[i*n+j+1], in[i*n+j-1], aboveE, belowE);
      error = max_error( error, out[i*n+j], in[i*n+j] );
    }
  return error;
}

void laplace_init(float *in, int n, int nprocs, int me)
{
  int i;
  const float pi  = 2.0f * asinf(1.0f);
  int block_size  = (int)(n/nprocs);
  memset(in, 0, (block_size+2)*n*sizeof(float));
  
  int source = me * block_size;
  int destination = me * block_size + block_size;
  
  for (i=source; i<destination; i++) {
    float V = in[(i % block_size)*n] = sinf(pi*i / (n-1));
    in[ (i % block_size)*n+n-1 ] = V*expf(-pi);
    if (i == source) {
      in[block_size*n] = V;
      in[block_size*n+n-1] = V*expf(-pi);
    } else if (i == destination - 1) {
      in[(block_size+1)*n] = V;
      in[(block_size+1)*n+n-1] = V*expf(-pi);
    }
  }
}

int main(int argc, char** argv)
{ 
  
  double t_before, t_after;
  int n = 4096;
  int iter_max = 1000;
  
  const float tol = 1.0e-5f;

  float *A, *temp;

  MPI_Request firstRowRequest, lastRowRequest;
  MPI_Status firstRowStatus, lastRowStatus;

  int err, me, nprocs;
  
  // get runtime arguments 
  if (argc>1) {  n        = atoi(argv[1]); }
  if (argc>2) {  iter_max = atoi(argv[2]); }

  MPI_Init( &argc, &argv );

  t_before = MPI_Wtime();

  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &me);

  int block_size                = (int)(n/nprocs);
  
  // we assumed that the row at block_size position will be data receiving from last row of the previous matrix
  // the row at block_size + 1 position will be data receiving from first row of the next matrix
  A    = (float*) malloc( (block_size+2)*n*sizeof(float) );
  temp = (float*) malloc( (block_size+2)*n*sizeof(float) ); 

  //  set boundary conditions
  laplace_init (A, n, nprocs, me);
  laplace_init (temp, n, nprocs, me);
  //  A[(n/128)*n+n/128] = 1.0f; // set singular point

  printf("Processor %d is running Jacobi relaxation Calculation: %d x %d mesh,"
         " maximum of %d iterations\n", 
        me, block_size, n, iter_max );
         
  int iter = 0;
  float my_error = 1.0f;

  while ( my_error > tol*tol && iter < iter_max )
  {
    iter++;
    if (me > 0) {
      // send first row and receive last row from the previous processor
      MPI_Send(A, n, MPI_FLOAT, me - 1, 0, MPI_COMM_WORLD);
      MPI_Irecv(A + block_size*n, n, MPI_FLOAT, me - 1, 0, MPI_COMM_WORLD, &firstRowRequest);
      MPI_Wait(&firstRowRequest, &firstRowStatus);
    }
    if (me < nprocs - 1) {
      // send last row and receive first row from the next processor
      MPI_Irecv(A + (block_size+1)*n, n, MPI_FLOAT, me + 1, 0, MPI_COMM_WORLD, &lastRowRequest);
      MPI_Send(A + (block_size-1)*n, n, MPI_FLOAT, me + 1, 0, MPI_COMM_WORLD);
      MPI_Wait(&lastRowRequest, &lastRowStatus);
    }
    my_error= laplace_step (A, temp, n, nprocs, me);   
    float *swap= A; A=temp; temp= swap; // swap pointers A & temp
    MPI_Barrier(MPI_COMM_WORLD);
  }

  free(A); free(temp);

  if (me > 0) {
    MPI_Send(&my_error, 1, MPI_FLOAT, 0, 0, MPI_COMM_WORLD);
  } else if (me == 0) {
    float receiver;
    for (int iMe=1; iMe<nprocs; iMe++) {
      MPI_Irecv(&receiver, 1, MPI_FLOAT, iMe, 0, MPI_COMM_WORLD, &lastRowRequest);
      MPI_Wait(&lastRowRequest, &lastRowStatus);
      my_error = max_error(my_error, 0, receiver);
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);

  err=MPI_Finalize();
  //check finalization

  if (me == 0) {
    my_error = sqrtf( my_error );
    printf("Total Iterations: %d, ERROR: %0.6f\n", iter, my_error);
    t_after = MPI_Wtime();
    printf("Total seconds the program took: %f\n", t_after - t_before);
  }
  
  return 0;
}
