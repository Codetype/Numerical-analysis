#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

const double eps = 1e-12;

void set_example_matrix(double *B, double ** A){
    A[0][0] = 4; A[0][1] = -2; A[0][2] = 2;  B[0] = -6;
    A[1][0] = -2; A[1][1] = 2; A[1][2] = 2;   B[1] = 4;
    A[2][0] = 2; A[2][1] = 2;  A[2][2] = 14;  B[2] = 0;
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

void Crout(double **A, double **L, double **U, int n) {
	int i, j, k;
	double sum = 0;

	for (i = 0; i < n; i++) {
		U[i][i] = 1;
	}

	for (j = 0; j < n; j++) {
		for (i = j; i < n; i++) {
			sum = 0;
			for (k = 0; k < j; k++) {
				sum = sum + L[i][k] * U[k][j];	
			}
			L[i][j] = A[i][j] - sum;
		}

		for (i = j; i < n; i++) {
			sum = 0;
			for(k = 0; k < j; k++) {
				sum = sum + L[j][k] * U[k][i];
			}
			if (L[j][j] == 0) {
				printf("Divide by 0...\n");
				exit(1);
			}
			U[j][i] = (A[j][i] - sum) / L[j][j];
		}
	}
    for(int i=0; i<n; i++){   
        for(int j=0; j<=i; j++){
            A[i][j] = L[i][j];
        }
        for(int j=i+1; j<n; j++){
            A[i][j] = U[i][j];
        }
    }
}

void cholesky(double **A, int n) {
    double **L = (double **)malloc(n * sizeof(double *)); 
    for (int i=0; i<n; i++) 
        L[i] = (double *)malloc((n) * sizeof(double)); 
 
    double s;
    for (int i = 0; i < n; i++) 
        for (int j = 0; j < (i+1); j++) {
            s = 0.0;
            for (int k = 0; k < j; k++){
                s += L[i][k] * L[j][k];
            }
            if(i == j) L[i][j] = sqrt(fabs(A[i][i] - s));
            else L[i][j] = (1.0 / L[j][j] * (A[i][j] - s));
        } 
        
    for(int i=0; i<n; i++){
        
        for(int j=0; j<i; j++){
            A[i][j] = L[i][j];
        }
        for(int j=i; j<n; j++){
            A[i][j] = L[j][i];
        }
    }
}

int solveEquation1(int n, double** A, double* B, double* X){
    double s;
    // Y(1) = B(1)
    X[0] = B[0];

    // Y(i) = B(i) - ((sum from j=1 to i-1) L(i,j) x Y(j) )
    for(int i=1; i<n; i++){
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

int solveEquation2(int n, double** A, double* B, double* X){
    double s;
    // Y(1) = B(1)/L(1)(1)
    X[0] = B[0]/A[0][0];

    // Y(i) = B(i) - ((sum from j=1 to i-1) L(i,j) x Y(j) )
    for(int i=1; i<n; i++){
        s = 0.0;
        for(int j=0; j<i; j++) s += A[i][j] * X[j]; 
        X[i] = B[i] - s;
    }

    // X(n) = Y(n)/U(n,n) 
    if(fabs(A[n-1][n-1]) < eps) return 0;
    X[n-1] /= A[n-1][n-1];

    // X(i) = (Y(i) - ((sum from j = i+1 to n) U(i,j) x X(j) )
    for(int i=n-2; i>=0; i--){
        s = 0.0;

        for(int j=i+1; j<n; j++) s += A[i][j] * X[j];

        if(fabs(A[i][i]) < eps) return 0;
        X[i] = (X[i] - s);
    }
    return 1;
}

int solveEquation3(int n, double** A, double* B, double* X){
    double s;
    // Y(1) = B(1)/L(1)(1)
    X[0] = B[0]/A[0][0];

    // Y(i) = B(i) - ((sum from j=1 to i-1) L(i,j) x Y(j) )
    for(int i=1; i<n; i++){
        s = 0.0;
        for(int j=0; j<i; j++) s += A[i][j] * X[j]; 
        X[i] = B[i] - s;
    }

    // X(n) = Y(n)/LT(n,n) 
    if(fabs(A[n-1][n-1]) < eps) return 0;
    X[n-1] /= A[n-1][n-1];

    // X(i) = 1/LT(i,i) * (Y(i) - ((sum from j = i+1 to n) LT(i,j) x X(j) )
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
  srand(time(NULL));
  // odczytujemy liczbę niewiadomych

  n = 3;

  // tworzymy macierze A, B i X

    double **A = (double **)malloc(n * sizeof(double *)); 
    for (i=0; i<n; i++) 
        A[i] = (double *)malloc((n) * sizeof(double)); 

    double **U = (double **)malloc(n * sizeof(double *)); 
    for (i=0; i<n; i++) 
        U[i] = (double *)malloc((n) * sizeof(double)); 

    double **L = (double **)malloc(n * sizeof(double *)); 
    for (i=0; i<n; i++) 
        L[i] = (double *)malloc((n) * sizeof(double)); 

    double *B = (double *)malloc((n) * sizeof(double)); 
    double *X = (double *)malloc((n) * sizeof(double)); 

    set_example_matrix(B,A);
    
    Doolitle(n,A);
    
    if(solveEquation1(n,A,B,X)){
        printf("\nDoolitle results:\n");
        for(int i=0; i<n; i++){
            printf("x%d = %lf\n", i+1, X[i]);
        }
    } else{
        printf("Error during calculation.");
    }

    set_example_matrix(B,A);

    Crout(A,L,U,n);

    if(solveEquation2(n,A,B,X)){
        printf("\nCrout results:\n");
        for(int i=0; i<n; i++){
            printf("x%d = %lf\n", i+1, X[i]);
        }
    } else{
        printf("Error during calculation.");
    }

    set_example_matrix(B,A);

    cholesky(A,n);

    if(solveEquation3(n,A,B,X)){
        printf("\nCholesky results:\n");
        for(int i=0; i<n; i++){
            printf("x%d = %lf\n", i+1, X[i]);
        }
    } else{
        printf("Error during calculation.");
    }

  for(i = 0; i < n; i++) free(A[i]);
  free(A);
  for(i = 0; i < n; i++) free(U[i]);
  free(U);
  for(i = 0; i < n; i++) free(L[i]);
  free(L);
  free(B);
  free(X);

  return 0;
} 