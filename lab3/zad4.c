#include <stdio.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <stdlib.h>
#include <time.h>

const double eps = 1e-12;

//correct solution
int gauss(int n, double ** AB, double * X)
{

  int i,j,k;
  double m,s;

  for(i = 0; i < n - 1; i++)
  {
    for(j = i + 1; j < n; j++)
    {
      if(fabs(AB[i][i]) < eps) return 0;
      m = -AB[j][i] / AB[i][i];
      for(k = i + 1; k <= n; k++)
        AB[j][k] += m * AB[i][k];
    }
  }

  for(i = n - 1; i >= 0; i--)
  {
    s = AB[i][n];
    for(j = n - 1; j >= i + 1; j--)
      s -= AB[i][j] * X[j];
    if(fabs(AB[i][i]) < eps) return 0;
    X[i] = s / AB[i][i];
  }
  return 1;
}

//gsl
void rand_gsl_equation_set(gsl_matrix* a, gsl_vector* b, int n, double** AB){
    for(int i=0; i<n; i++){
        double v = (rand()%100000)/10000.0;
        gsl_vector_set(b, i, v);
        AB[i][n] = v;
        for(int j=0; j<n; j++){
            v = (rand()%100000)/10000.0;
            AB[i][j] = v;
            gsl_matrix_set(a, i, j, v);
        }
    }
}

int main(int argc, const char * argv[])
{
    srand(time(NULL));
    int s, n = 100;
    //for(int i=0; i<10; i++){

      //gsl matrix and vectors declaration and allocation
      gsl_vector *b, *x;
      gsl_matrix *a;
      gsl_permutation *p;

      b = gsl_vector_alloc (n);
      a = gsl_matrix_alloc(n, n);

      x = gsl_vector_alloc(n);
      p = gsl_permutation_alloc(n);

      //double matrix declaration and allocation
      double **AB = (double **)malloc(n * sizeof(double *));
      for (int i=0; i<n; i++)
          AB[i] = (double *)malloc((n+1) * sizeof(double));

      double *X = (double *)malloc((n) * sizeof(double));

      //fill gsl and double matrix with the same data
      rand_gsl_equation_set(a,b,n,AB);

      clock_t start = clock();
      gsl_linalg_LU_decomp(a, p, &s);
      gsl_linalg_LU_solve(a, p, b, x);
      clock_t end = clock();
      double time1  = (double)(end-start)/CLOCKS_PER_SEC;

      start = clock();
      if(!gauss(n,AB,X))
          return 1;
      end = clock();
      double time2  = (double)(end-start)/CLOCKS_PER_SEC;

      printf("%lf %lf\n", time1, time2);

      free(X);
      for(int i=0; i<n; i++)
          free(AB[i]);
      free(AB);
      gsl_vector_free(b);
      gsl_vector_free(x);
      gsl_permutation_free(p);
      gsl_matrix_free(a);
    //  n += 100;
    //}

    return 0;
}
