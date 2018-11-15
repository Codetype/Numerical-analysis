#include <stdio.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <stdlib.h>
#include <time.h>

#define MAX_ITER 10000
const double eps = 1e-12;

void gauss_seidel();
void init();
int convergence();
void print_vector();
void print_equation();
void print_gsl_equation_set();
void print_gsl_vector();
void fill_gsl_equation_set();

int main(int argc, const char * argv[])
{
    srand(time(NULL));

    int n = atoi(argv[1]);
    double error = atof(argv[2]);

    int s;

      //gsl matrix and vectors declaration and allocation
      gsl_vector *b, *x;
      gsl_matrix *a;
      gsl_permutation *p;

      b = gsl_vector_alloc (n);
      a = gsl_matrix_alloc(n, n);

      x = gsl_vector_alloc(n);
      p = gsl_permutation_alloc(n);

      //double matrix declaration and allocation
      double **A = (double **)malloc(n * sizeof(double *));
      for (int i=0; i<n; i++)
          A[i] = (double *)malloc((n) * sizeof(double));

      double *B = (double *)malloc((n) * sizeof(double));
      double *X = (double *)malloc((n) * sizeof(double));

      //fill gsl and double matrix with the same data
      init(A, B, X, n);
      fill_gsl_equation_set(a,b,n,A,B);
      print_gsl_equation_set(a,b,n);
      
      gsl_linalg_LU_decomp(a, p, &s);
      gsl_linalg_LU_solve(a, p, b, x);      
      print_gsl_vector(x,n);


      print_equation(A,B,n);
      gauss_seidel(A,B,X,n,error);
      print_vector(X,n);
      

      free(X);
      free(B);
      for(int i=0; i<n; i++)
          free(A[i]);
      free(A);
      gsl_vector_free(b);
      gsl_vector_free(x);
      gsl_permutation_free(p);
      gsl_matrix_free(a);

    return 0;
}

void gauss_seidel(double **A, double* b, double* x, int n, double error){
  int i,j,k = 0;

      double **L = (double **)malloc(n * sizeof(double *));
      for (int i=0; i<n; i++)
          L[i] = (double *)malloc((n) * sizeof(double));
      double **U = (double **)malloc(n * sizeof(double *));
      for (int i=0; i<n; i++)
          U[i] = (double *)malloc((n) * sizeof(double));
      double **D = (double **)malloc(n * sizeof(double *));
      for (int i=0; i<n; i++)
          D[i] = (double *)malloc((n) * sizeof(double));
    
    // Divide A into L + D + U
    for (i=0; i<n; i++)
        for (j=0; j<n; j++) {
            if (i < j) {
                U[i][j] = A[i][j];
            }
            else if (i > j) {
                L[i][j] = A[i][j];
            }
            else {
                D[i][j] = A[i][j];
            }
        }
    
    // Calculate D^-1
    for (i=0; i<n; i++)
        D[i][i] = 1/D[i][i];
    
    // Calculate D^-1 * b
    for (i=0; i<n; i++)
        b[i] *= D[i][i];
    
    //Calculate D^-1 * L
    for (i=0; i<n; i++)
        for (j=0; j<i; j++)
            L[i][j] *= D[i][i];
    
    //Calculate D^-1 * U
    for (i=0; i<n; i++)
        for (j=i+1; j<n; j++)
            U[i][j] *= D[i][i];
    
    //Initialize x
    for (i=0; i<n; i++)
        x[i] = 0;
    
    while (!convergence(k,A,b,x,error) && (k < MAX_ITER)) {
        for (i=0; i<n; i++) {
            x[i] = b[i];                       // x = D^-1*b -
            for (j=0; j<i; j++)
                x[i] -= L[i][j]*x[j];    // D^-1*L * x -
            for (j=i+1; j<n; j++)
                x[i] -= U[i][j]*x[j];    // D^-1*U * x
        }
        k++;
    }
}

// returns 1 if converged else 0
int convergence(int iter, double **A, double* B, double* X, int n, double error){
  int i,j,flag=1;
  float sum;
  for (i = 0; i < n; i++) {
    sum = 0.0;
    for(j = 0; j < n; j++){
      sum = sum + A[i][j] * X[j];
      if ((sum - B[i]) < error) {
      	flag = 0;
      	return flag;
      }
    }
  }

  return flag;
}

void init(double **A, double* B, double* X, int n){
  int i,j;
  for (i=0;i<n;i++) {
    for (j=0;j<n;j++) {
      if (i == j) A[i][j] = ((rand()%100000)/10000.0);
      //else if(i == j+1) A[i][j] = (rand()%100000)/10000.0;
      else A[i][j] = ((rand()%1000)/100000.0);
    }
  }

  for (i=0;i<n;i++) X[i] = 1.0;


  for (i=0;i<n;i++){
    B[i] = (rand()%100000)/10000.0;
  }

}

void print_equation(double **A, double* B, int n){
  int i,j;
    printf("\nGauss Seidel method\n");
  for (i=0;i<n;i++) {
    printf("|");
    for (j=0;j<n;j++) printf(" %lf ",(A[i][j]));
    printf("|| x%d | | %lf |\n",i+1,(B[i]));
  }
  printf("\n");
}

void print_vector(double *l, int n){
  int i;
  for (i=0; i<n; i++) printf("x%d = %lf\n",i+1,l[i]);
  printf("\n");
}

//gsl
void fill_gsl_equation_set(gsl_matrix* a, gsl_vector* b, int n, double** A, double* B){
    for(int i=0; i<n; i++){
        gsl_vector_set(b, i, B[i]);
        for(int j=0; j<n; j++){
            gsl_matrix_set(a, i, j, A[i][j]);
        }
    }
}

void print_gsl_equation_set(gsl_matrix* a, gsl_vector* b, int n){
      printf("\nLibrary method\n");
    for(int i=0; i<n; i++){ printf("|");
        for(int j=0; j<n; j++){
            printf(" %lf ", gsl_matrix_get(a, i, j));
        }    printf("|| x%d | | %lf |\n",i+1, gsl_vector_get(b, i));
    }
}

void print_gsl_vector(gsl_vector* x, int n){
    printf("\n");
    for(int i=0; i<n; i++){
      printf("x%d = %lf\n",i+1, gsl_vector_get(x, i));
    }
}