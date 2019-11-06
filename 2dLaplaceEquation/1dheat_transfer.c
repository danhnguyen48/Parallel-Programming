#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>

#define n 10
#define T 40

double U[T][n+2];
const double L = 0.345678;

void initialazation();
double randomTheTemperature(double *pointOfBar, double minTemp, double maxTemp);
double doing();


int main(int agrv, char *argv[]) {

    initialazation();
    const double result = doing();
    printf("\nResult: %.2lf\n", result);

}

void initialazation() {

    for (int t=0; t<T; t++) {
        U[t][0] = 0;
        U[t][n+1] = 0;
    }
    for (int i=1; i<n+1; i++) {
        U[0][i] = randomTheTemperature(U[0] + i, 0.0, 100.0);
    }

}

double randomTheTemperature(double *pointOfBar, double minTemp, double maxTemp) {

    double range = maxTemp - minTemp;
    double div = RAND_MAX / range;
    double randomValue;
    randomValue = (double) (rand() / div) + minTemp;
    return randomValue;

}

double doing() {

    double result = 0.0;

    for (int t=1; t<T; t++)
        for (int x=1; x<n+1; x++) {
            U[t][x] = (1.0 - 2 * L) * U[t-1][x] + L * (U[t-1][x-1] + U[t-1][x+1]);
            if (t == T-1) {
                printf("%.2lf ", U[t][x]);
                result += U[t][x];
            }
        }

    return result;

}