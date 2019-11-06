#include <stdio.h>

#define NumberOfElements 9

float X[NumberOfElements] = {10.0, 4.0, 12.0, 6.0, 3.0, 52.0, 7.0, 33.0, 65.0};
float Y[NumberOfElements];

void mergeSort(float X[], float Y[], int n);
void copy(float X, float *Y);
void merge(float Y[], float T1[], float T2[], int n1, int n2);

int main(int argc, char *argv[]) {

    mergeSort(X, Y, NumberOfElements);
    for (int i=0; i<NumberOfElements; i++) {
        printf("%f ", Y[i]);
    }
    printf("\n");

}

void mergeSort(float X[], float Y[], int n) {

    if (n == 1) {
        copy(*X, &(*Y));
        return;
    }
    float T[n];
    mergeSort(X, T, n/2);
    mergeSort(X + n/2, (T + n/2), n - n/2);
    merge(Y, T, T + n/2, n/2, n - n/2);

}

void copy(float X, float *Y) {

    *Y = X;

}

void merge(float Y[], float T1[], float T2[], int n1, int n2) {

    int i = 0;
    int j = 0;
    int w = 0;
    while (i < n1 || j < n2 ) {

        if (j >= n2 || (i < n1 && *(T1 + i) <= *(T2 + j))) {
            *(Y + w++) = *(T1 + i++);
        } else {
            *(Y + w++) = *(T2 + j++);
        }

    }

}