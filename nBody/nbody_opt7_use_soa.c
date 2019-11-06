// BASELINE VERSION
//
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

typedef float real;
#define SOFTENING_SQUARED  0.01
#define N 20000

typedef struct { real x[N], y[N], z[N], w[N]; } population4;
typedef struct { real x[N], y[N], z[N]; } population3;
typedef struct { real x, y, z; } instance;

void recalculatePosition(population3 force, population4 * in, population4 * out, population3 * vel, int i, real dt)
{
    real fx = force.x[i], fy = force.y[i], fz = force.z[i];
    real px = in->x[i],    py = in->y[i],    pz = in->z[i],    invMass = in->w[i];
    real vx = vel->x[i],   vy = vel->y[i],   vz = vel->z[i];

    // acceleration = force / mass; 
    // new velocity = old velocity + acceleration * deltaTime
    vx += (fx * invMass) * dt;
    vy += (fy * invMass) * dt;
    vz += (fz * invMass) * dt;

    // new position = old position + velocity * deltaTime
    px += vx * dt;
    py += vy * dt;
    pz += vz * dt;

    out->x[i] = px;
    out->y[i] = py;
    out->z[i] = pz;
    out->w[i] = invMass;

    vel->x[i] = vx;
    vel->y[i] = vy;
    vel->z[i] = vz;
}


void integrate(population4 * out, population4 * in,
               population3 * vel, population3 * force,
               real    dt,  int n)
{
  int i, j;
  memset(force, 0.0, n*sizeof(population3));
  for (i = 0; i < n; i++)
  {

    for (j = i+1; j < n; j++)
    {

      real rx, ry, rz, distSqr, s;

      rx = in->x[j] - in->x[i];  ry = in->y[j] - in->y[i];  rz = in->z[j] - in->z[i];

      distSqr = rx*rx+ry*ry+rz*rz;

      if (distSqr < SOFTENING_SQUARED) s = 1.0/sqrtf(SOFTENING_SQUARED);
      else s = 1.0/sqrtf(distSqr);
      s = s*s*s;
      s = in->w[j] * s; 

      force->x[i] += rx * s;   force->y[i] += ry * s;   force->z[i] += rz * s;
      force->x[j] += -rx * s;  force->y[j] += -ry * s;  force->z[j] += -rz * s;

    }

    recalculatePosition(*force, in, out, vel, i, dt);
  }
}


real dot(real v0[3], real v1[3])
{
  return v0[0]*v1[0]+v0[1]*v1[1]+v0[2]*v1[2];
}


real normalize(real vector[3])
{
  float dist = sqrt(dot(vector, vector));
  if (dist > 1e-6)
  {
    vector[0] /= dist;
    vector[1] /= dist;
    vector[2] /= dist;
  }
  return dist;
}


void cross(real out[3], real v0[3], real v1[3])
{
  out[0] = v0[1]*v1[2]-v0[2]*v1[1];
  out[1] = v0[2]*v1[0]-v0[0]*v1[2];
  out[2] = v0[0]*v1[1]-v0[1]*v1[0];
}


void randomizeBodies(population4* pos, 
                     population3* vel, 
                     float clusterScale, 
                     float velocityScale, 
                     int   n)
{
  srand(42);
  float scale = clusterScale;
  float vscale = scale * velocityScale;
  float inner = 2.5f * scale;
  float outer = 4.0f * scale;

  int p = 0, v=0;
  int i = 0;
  while (i < n)
  {
    real x, y, z;
    // generate real numbers between -1.0 and +1.0
    x = rand() / (float) RAND_MAX * 2 - 1;
    y = rand() / (float) RAND_MAX * 2 - 1;
    z = rand() / (float) RAND_MAX * 2 - 1;

    real point[3] = {x, y, z};
    real len = normalize(point);
    if (len > 1) // discard position and generate new one
      continue;

    pos->x[i]= point[0] * (inner + (outer - inner)*rand() / (real) RAND_MAX);
    pos->y[i]= point[1] * (inner + (outer - inner)*rand() / (real) RAND_MAX);
    pos->z[i]= point[2] * (inner + (outer - inner)*rand() / (real) RAND_MAX);
    pos->w[i]= 1.0f;

    real axis[3] = {0.0f, 0.0f, 1.0f};

    if (1 - dot(point, axis) < 1e-6)
    {
      axis[0] = point[1];
      axis[1] = point[0];
      normalize(axis);
    }
    real vv[3] = {(real)pos->x[i], (real)pos->y[i], (real)pos->z[i]};
    real vv0[3];

    cross(vv0, vv, axis);
    vel->x[i] = vv0[0] * vscale;
    vel->y[i] = vv0[1] * vscale;
    vel->z[i] = vv0[2] * vscale;

    i++;
  }
}


instance average(population4 * p, int n)
{
  int i;
  instance av= {0.0, 0.0, 0.0};
  for (i = 0; i < n; i++)
  {
    av.x += p->x[i];
    av.y += p->y[i];
    av.z += p->z[i];
  }
  av.x /= n;
  av.y /= n;
  av.z /= n;
  return av;
}


int main(int argc, char** argv)
{
  int i, j, n = N;
  int iterations = 10;
  real dt = 0.01667;

  if (argc >= 2) n = atoi(argv[1]);
  if (argc >= 3) iterations = atoi(argv[2]);

  population4 *pin = (population4*)malloc(sizeof(population4));
  population4 *pout = (population4*)malloc(n * sizeof(population4));
  population3 *v    = (population3*)malloc(n * sizeof(population3));
  population3 *f    = (population3*)malloc(n * sizeof(population3));

  // Third parameter is cluster scale, Fouth parameter is velocity scale
  randomizeBodies(pin, v,  1.54f, 8.0f, n);

  printf("n=%d bodies for %d iterations:\n", n, iterations);

  for (i = 0; i < iterations; i++)
  {
    if (i%2 == 0)
        integrate (pout, pin, v, f, dt, n); 
    else 
        integrate (pin, pout, v, f, dt, n); 
  }

  instance p_av= average( pin, n);
  printf("Average position: (%f,%f,%f)\n", p_av.x, p_av.y, p_av.z);
  printf("Body-0  position: (%f,%f,%f)\n", pin->x[0], pin->y[0], pin->z[0]);

  free(pin);  free(pout);  free(v);  free(f);

  return 0;
}