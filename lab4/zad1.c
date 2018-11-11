#include <stdio.h>
#include <math.h>
#include <stdlib.h>

const double eps = 1e-12;

void set_matrix(double *B, double ** A){
    A[0][0] = 4; A[0][1] = -2; A[0][2] = 4; A[0][3] = -2; B[0] = 8;
    A[1][0] = 3; A[1][1] = 1;  A[1][2] = 4; A[1][3] = 2;  B[1] = 7;
    A[2][0] = 2; A[2][1] = 4;  A[2][2] = 2; A[2][3] = 1;  B[2] = 10;
    A[3][0] = 2; A[3][1] = -2; A[3][2] = 4; A[3][3] = 2;  B[3] = 2;
}

int Doolitle(int n, double** A){
    double s;
    for(int j=0; j<n; j++){
        //0 on diagonal
        if(fabs(A[j][j] < eps)) return 0;
        //set u matrix
        for(int i = 0; i<=j; i++){
            s = 0.0;
            for(int k = 0; k<i; k++) s += A[i][k] * A[k][j];
            A[i][j] -= s;
        }
        //set l matrix
        for(int i=j+1; i<n; i++){
            s = 0.0;
            for(int k = 0; k<j; k++) s += A[i][k] * A[k][j];
            A[i][j] = (A[i][j] - s) / A[j][j];
        }
    }
    return 1;
}

int solveEquation(int n, double** A, double* B, double* X){
    double s;
    // Y(1) = B(1)
    X[0] = B[0];

    // Y(i) = B(i) - ((sum from j=1 to i-1) L(i,j) x Y(j) )
    for(int i=0; i<n; i++){
        s = 0.0;
        for(int j=0; j<i; j++) s += A[i][j] * X[j]; 
        X[i] = B[i] - s;
    }

    // X(n) = Y(n)/U(n,n) 
    if(fabs(A[n-1][n-1]) < eps) return 0;
    X[n-1] /= A[n-1][n-1];

    // X(i) = 1/U(i,i) * (Y(i) - ((sum from j = i+1 to n) U(i,j) x X(j) )
    for(int i=n-2; i>=0; i--){
        s = 0.0;

        for(int j=i+1; j<n; j++) s += A[i][j] * X[j];

        if(fabs(A[i][i]) < eps) return 0;
        X[i] = (X[i] - s)/A[i][i];
    }
    return 1;
}

int  main()
{
    int n,i,j;
    n = 4 ;
    double **A = (double **)malloc(n * sizeof(double *)); 
    for (i=0; i<n; i++) 
        A[i] = (double *)malloc((n) * sizeof(double)); 

    double *B = (double *)malloc((n) * sizeof(double)); 
    double *X = (double *)malloc((n) * sizeof(double)); 

    set_matrix(B,A);

    if(Doolitle(n, A) && solveEquation(n,A,B,X)){
        for(int i=0; i<n; i++){
            printf("%lf\n", X[i]);
        }
    } else{
        printf("Error during calculation.");
    }


    for(i = 0; i < n; i++) free(A[i]);
    free(A);
    free(B);
    free(X);

    return 0;
} 