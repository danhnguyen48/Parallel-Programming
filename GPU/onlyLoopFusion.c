#include <stdlib.h>

void loop_fusion(float * restrict x, float * y, int N, float a1, float a2, float a3, float a4, float a5, float a6)
{
  #pragma acc kernels present(x,y)
  for (int i = 0; i < N; i++)
  {
    x[i] += a1 + i*a3;
    y[i] = (y[i] + a2 + i*a4 + a5*x[i])*a6;
  }
} 

void loop_4(float * x, int N, float a)
{
  #pragma acc kernels present(x)
  for (int i = 0; i < N-1; i++) 
    x[i] = a * x[i+1] ;
} 

void loop_5(float *x, int N, float a)
{
  #pragma acc kernels present(x)
  for (int i = 1; i < N; i++) 
    x[i] = a * x[i-1] ;
} 

float fuse_loop_67(float * restrict x, float *y, int N, float a9)
{
  float S=0.0;
  #pragma acc kernels present(x,y)
  for (int i = 0; i < N; i++)
  {
    if (x[i] > a9)
      x[i] = 0;
    S += x[i] + y[i];
  }
  return S;
}

void initVec ( float *x, int N, int seed ) 
{
  srand48 (seed);  // set initial seed for random number generation
  for (int i=0; i<N; i++)
    x[i]=drand48();
}

int main(int argc, char **argv) 
{
  int N= 1000000;
  float a1=3.0, a2=5.2, a3= 4.0, a4= -3.7, a5= 0.95, a6= 1.12, a7= -0.91, a8= 0.98, a9= 0.8;

  float S;
  float * V = (float*)malloc(N * sizeof(float));
  float * W = (float*)malloc(N * sizeof(float));

  initVec (V, N, 0);

  #pragma acc data copyin(V[:N],W[:N])
  {
    loop_fusion(V, W, N, a1, a2, a3, a4, a5, a6);
    loop_4(V, N, a7);
    loop_5(W, N, a8);
    S = fuse_loop_67(V, W, N, a9);
  }

  printf("Checksum= %12.7f\n", S);
  return 0;
}
