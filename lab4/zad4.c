#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

const double eps = 1e-12;
void rand_matrix(int n, double **AB);
void set_matrix(int n, double **AB, double **A, double *B);

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

int  main(int argc, char* argv[])
{
  int n,i,j;
  srand(time(NULL));

    if(argc < 2){
        printf("Wrong number of arguments! Set equation range.\n");
        return 1;
  }

  n = atoi(argv[1]);

    double **AB = (double **)malloc(n * sizeof(double *)); 
    for (i=0; i<n; i++) 
        AB[i] = (double *)malloc((n+1) * sizeof(double)); 

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

    rand_matrix(n,AB);
    
    set_matrix(n,AB,A,B);

    clock_t start = clock();
    Doolitle(n,A);
    
    if(!solveEquation1(n,A,B,X)){
        printf("Error during calculation.");
    }
    clock_t end = clock();
    double time1 = (double)(end - start)/CLOCKS_PER_SEC;

    set_matrix(n,AB,A,B);
    
    start = clock();
    Crout(A,L,U,n);
    
    if(!solveEquation2(n,A,B,X)){
        printf("Error during calculation.");
    }
    end = clock();
    double time2 = (double)(end - start)/CLOCKS_PER_SEC;

    set_matrix(n,AB,A,B);
    
    start = clock();
    cholesky(A,n);
    
    if(!solveEquation3(n,A,B,X)){
        printf("Error during calculation.");
    }
    end = clock();
    double time3 = (double)(end - start)/CLOCKS_PER_SEC;

    printf("%lf\n%lf\n%lf\n", time1, time2, time3);

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

void rand_matrix(int n, double **AB){
    for(int i=0; i<n; i++){
        for(int j=0; j<=n; j++){
            AB[i][j] = rand()%100000/10000.0;
        }
    }
}

void set_matrix(int n, double **AB, double **A, double *B){
    int i,j;
    for(i=0; i<n; i++){  
        for(j=0; j<n; j++){
            A[i][j] = AB[i][j];
        }
        B[i] = AB[i][j+1];
    }
}
