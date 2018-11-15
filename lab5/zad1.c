#include <stdio.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <stdlib.h>
#include <time.h>

#define MAX_ITER 10000
const double eps = 1e-12;

int jacobi();
int convergence();
void print_vector();
void print_equation();
void rand_gsl_equation_set();
void print_gsl_equation_set();
void print_gsl_vector();

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
      rand_gsl_equation_set(a,b,n,A,B,X);
      print_gsl_equation_set(a,b,n);
      
      gsl_linalg_LU_decomp(a, p, &s);
      gsl_linalg_LU_solve(a, p, b, x);      
      print_gsl_vector(x,n);


      print_equation(A,B,n);
      jacobi(A,B,X,n,error);



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

int jacobi(double **A, double* B, double* X, int n, double error){
  int i,j,k = 0;
  double sum;
  double *buff = (double *)malloc((n) * sizeof(double));

  while (!convergence(k,A,B,X,error) && (k < MAX_ITER)) {
        for (i = 0; i < n; i++) {
            sum = 0.0;
            for (j = 0; j < n; j++) {
                if (i != j)
                    sum = sum + A[i][j] * X[j];
            }
            buff[i] = (B[i] - sum) / A[i][i];
        }
        k++;
  }

    print_vector(buff,n);
    return k;
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

void print_equation(double **A, double* B, int n){
  int i,j;
    printf("\nJacobi method\n");
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
void rand_gsl_equation_set(gsl_matrix* a, gsl_vector* b, int n, double** A, double* B, double* X){
    for(int i=0; i<n; i++){
        double v = (rand()%100000)/10000.0;
        gsl_vector_set(b, i, v);
        B[i] = v;
        for(int j=0; j<n; j++){
            v = (rand()%100000)/10000.0;
            A[i][j] = v;
            gsl_matrix_set(a, i, j, v);
        }
    }

    for (int i=0;i<n;i++) X[i]=1;
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