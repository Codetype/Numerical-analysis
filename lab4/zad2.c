#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

const double eps = 1e-12;

void set_example_matrix(double *B, double ** A){
    A[0][0] = 4; A[0][1] = -2; A[0][2] = 4; A[0][3] = -2; B[0] = 8;
    A[1][0] = 3; A[1][1] = 1;  A[1][2] = 4; A[1][3] = 2;  B[1] = 7;
    A[2][0] = 2; A[2][1] = 4;  A[2][2] = 2; A[2][3] = 1;  B[2] = 10;
    A[3][0] = 2; A[3][1] = -2; A[3][2] = 4; A[3][3] = 2;  B[3] = 2;
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


/*
Gauss algorithm
*/

void swap(int *a, int *b){
  int tmp = *a;
  *a = *b;
  *b = tmp;
}

int Gauss(int n, double ** A, double *B, double * X, int * W)
{
  int   i,j,k;
  double m,s;

  for(i = 0; i <= n; i++) W[i] = i;

  for(i = 0; i < n - 1; i++)
  {
    k = i;
    for(j = i + 1; j < n; j++)
      if(fabs(A[i][W[k]]) < fabs(A[i][W[j]])) k = j;
    swap(&W[k], &W[i]);
    for(j = i + 1; j < n; j++)
    {

      if(fabs(A[i][W[i]]) < eps) return 0;
      m = -A[j][W[i]] / A[i][W[i]];
      for(k = i + 1; k < n; k++)
        A[j][W[k]] += m * A[i][W[k]];
      B[j] += m * B[i]; 
    }
  }

  for(i = n - 1; i >= 0; i--)
  {
    if(fabs(A[i][W[i]]) < eps) return 0;
    s = B[i];
    for(j = n - 1; j >= i + 1; j--)
      s -= A[i][W[j]] * X[W[j]];
    X[W[i]] = s / A[i][W[i]];
  }
  return 1;
}


/*
End of Gauss algorithm
*/

int  main(int argc, char *argv[])
{
  srand(time(NULL));
  int n,i,j;
  if(argc < 2){
      printf("Wrong number of arguments! Set equation range.\n");
      return 1;
  }

  n = atoi(argv[1]);

    double **A = (double **)malloc(n * sizeof(double *)); 
    for (i=0; i<n; i++) 
        A[i] = (double *)malloc((n) * sizeof(double)); 

    double *B = (double *)malloc((n) * sizeof(double)); 
    double *X = (double *)malloc((n) * sizeof(double)); 
    int *W = (int *)malloc((n) * sizeof(int));

    double **AB = (double **)malloc(n * sizeof(double *)); 
    for (i=0; i<n; i++) 
        AB[i] = (double *)malloc((n+1) * sizeof(double)); 

    double **Test = (double **)malloc(n * sizeof(double *)); 
    for (i=0; i<n; i++) 
        Test[i] = (double *)malloc(2 * sizeof(double)); 

    rand_matrix(n,AB);
    set_matrix(n,AB,A,B);

    clock_t start = clock();
    if(Doolitle(n, A) && solveEquation(n,A,B,X)){
        //printf("\nDoolitle results:\n");
        for(int i=0; i<n; i++){
            //printf("x%d = %lf\n", i+1, X[i]);
            Test[i][0] = X[i];
        }
    } else{
        printf("Error during calculation.");
    }
    clock_t end = clock();
    double time1 = (double)(end - start)/CLOCKS_PER_SEC;

    set_matrix(n,AB,A,B);
    
    start = clock();  
    if(Gauss(n,A,B,X,W)){
        //printf("\nGauss results:\n");
        for(int i=0; i<n; i++){
            //printf("x%d = %lf\n", i+1, X[i]);
            Test[i][1] = X[i];
        }
    } else{
        printf("Error during calculation.");
    }
    end = clock();
    double time2 = (double)(end - start)/CLOCKS_PER_SEC;

    //printf("Time1: %lf s\nTime2: %lf s\n", time1, time2);

    int flag = 1;
    for(int i=0; i<n;i++){
        if(fabs(Test[i][0] - Test[i][1]) <= eps) continue;
        else flag = 0;
    }
    if(flag == 1) printf("Test passed\n");
    else printf("Test failed\n");
    
  for(i = 0; i < n; i++) free(A[i]);
  free(A);
  free(B);
  free(X);

  return 0;
} 