#include <stdio.h>
#include <math.h>

#define X 3
#define T 4
#define L 0.345678

double U[X+1][T+2];

void initialazation();
void doing(double *sum);
double calculateTheEquation(int x, int t);


int main(int argv, char *argc[]) {

    double checksum;

    initialazation();
    doing(&checksum);
    printf("%lf\n", checksum);

}


void initialazation() {

    // initialize beginning position of rope
    for (int i=0; i<X; i++) {
        U[i][0] = sin(i*M_PI / X);
        U[i][1] = sin(i*M_PI / X) * cos(M_PI / T);
    }
    // initialize the boundary matrix
    for (int i=0; i<T+2; i++) {
        U[0][i] = 0;
        U[X][i] = 0;
    }

}

void doing(double *sum) {

    *sum = 0;
    for (int x=0; x<X+1; x++) {
        for (int t=0; t<T+2; t++) {
            U[x][t+1] = calculateTheEquation(x, t);
        }
        *sum += U[x][T+1];
    }

}

double calculateTheEquation(int x, int t) {

    return 2 * (1-L) * U[x][t] + L * U[x+1][t] + L * U[x-1][t] - U[x][t-1];

}