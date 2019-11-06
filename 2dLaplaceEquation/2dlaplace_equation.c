#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define pi 3.141516
#define e 2.7182818284590423536
// #define n 10
// #define m 15
#define iter_max 100

float error = 0;
float tol = 0.0000001;
int n, m;

void calculateError(float valueOfAElement, float valueOfAnewElement);
void initialization(float A[][m], float Anew[][m]);
void printA();
void printAnew();

int main(int argc, char *argv[]) {

    int iter = 0;

    if (argc < 3 || atoi(argv[1]) <= 1 || atoi(argv[2]) <= 1) {
        // display the error message
        printf("Error! We will assign n to 10 and m to 15.\n");
        n = 10;
        m = 15;
    } else {
        n = atoi(argv[1]);
        m = atoi(argv[2]);
    }
    float A[n][m], Anew[n][m];

    initialization(A, Anew);

    do {

        // First loop that calculates Anew 
        for (int i=1; i<n-1; i++) 
            for (int j=1; j<m-1; j++) {
                Anew[i][j] = (A[i][j-1] + A[i+1][j] + A[i][j+1] + A[i-1][j]) / 4;
                // printf("A[%d][%d] = %f\n", i, j, (A[i][j-1] + A[i+1][j] + A[i][j+1] + A[i-1][j]) / 4);
            }

        // Second loop that calculates error
        error = 0;
        for (int i=1; i<n-1; i++) 
            for (int j=1; j<m-1; j++) {
                calculateError(A[i][j], Anew[i][j]);
            }

        // Third loop that updates A from Anew
        for (int i=1; i<n-1; i++) 
            for (int j=1; j<m-1; j++) {
                A[i][j] = Anew[i][j];
            }
        if (++iter % 10 == 0)
            printf("This is a error at iteration %d: %f\n", iter, error);

    } while (error > tol && iter <= iter_max);

    printf("The interation is: %d\n", iter);
    for (int i=0; i<n; i++) {
        for (int j=0; j<m; j++) {
            printf("%.2f ", A[i][j]);
        }
        printf("\n");
    }

}


void initialization(float A[][m], float Anew[][m]) {
    for (int i=0; i<n; i++)
        for (int j=0; j<m; j++) {
            if (j == 0) {
                Anew[i][j] = sin(i * pi / (n-1));
            } else if (j == m-1) {
                Anew[i][j] = sin(i * pi / (n-1)) * (pow(e, -pi));
            } else {
                Anew[i][j] = 0;
            }
            A[i][j] = Anew[i][j];
            // calculateError(A[i][j], 0);
        }
    return;
}

void calculateError(float valueOfAElement, float valueOfAnewElement) {
    const float tempError = sqrt(fabsf(valueOfAElement - valueOfAnewElement));
    if (tempError > error) error = tempError;
    return;
}