#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

const float eps = 1e-6;

//naive solution
void naive_gauss(int n, float ** AB, float * X)
{

  int i,j,k;
  float m,s;

  for(i = 0; i < n - 1; i++)
  {
    for(j = i + 1; j < n; j++)
    {
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

    X[i] = s / AB[i][i];
  }
}


void fill_matrix1(float** AB, float* R);
void fill_matrix2(float** AB, float* R);
void fill_matrix3(float** AB, float* R);

//testing
void test(float** AB, float* X, float* R, int n){
    naive_gauss(n, AB, X);

    int test_flag = 1;
    for(int i=0; i<n; i++){
        if(fabs(X[i]-R[i]) < eps) continue;
        else test_flag = 0;
    }

    if(test_flag) printf("Test passed\n");
    else printf("Test failured\n");
}

int main(int argc, const char * argv[])
{
    srand(time(NULL));
    int n=4;
    //float matrix declaration and allocation
    float **AB = (float **)malloc(n * sizeof(float *));
    for (int i=0; i<n; i++)
        AB[i] = (float *)malloc((n+1) * sizeof(float));

    float *X = (float *)malloc((n) * sizeof(float));
    float *R = (float *)malloc((n) * sizeof(float));

    /*
    Input data:
    0.00*x0 + 2.00*x1 + 3.00*x2 + 4.00*x3 = 49.00
    1.00*x0 + 0.00*x1 + 3.00*x2 + 4.00*x3 = 45.00
    1.00*x0 + 2.00*x1 + 0.00*x2 + 4.00*x3 = 36.00
    1.00*x0 + 2.00*x1 + 3.00*x2 + 0.00*x3 = 23.00

    Expected result
    { x1 = 2
    { x2 = 3
    { x3 = 5
    { x4 = 7
    */

    fill_matrix1(AB, R);
    test(AB, X, R, n);

    /*
    Input data:
    10.000000*x0 + -7.000000*x1 + 0.000000*x2 = 7.000000
    0.000000*x0 + -0.001000*x1 + 6.000000*x2 = 6.001000
    0.000000*x0 + 2.500000*x1 + 5.000000*x2 = 2.500000

    Expected result
    {x0 = 0.0
    {x1 = -1.0
    {x2 = 1.0
    */
    fill_matrix2(AB, R);
    test(AB, X, R, n-1);

    /*
    Input data:
    3.000000*x0 + -13.000000*x1 + 9.000000*x2 + 3.000000*x3 = -19.000000
    -6.000000*x0 + 4.000000*x1 + 1.000000*x2 + -18.000000*x3 = -34.000000
    6.000000*x0 + -2.000000*x1 + 2.000000*x2 + 4.000000*x3 = 16.000000
    12.000000*x0 + -8.000000*x1 + 6.000000*x2 + 10.000000*x3 = 26.000000

    Expected result
    { x1 = 3
    { x2 = 1
    { x3 = -2
    { x4 = 1
    */
    fill_matrix3(AB, R);
    test(AB, X, R, n);

    //free memory
    free(X);
    free(R);
    for(int i=0; i<n; i++)
        free(AB[i]);
    free(AB);

    return 0;
}

//fill matrices
void fill_matrix1(float** AB, float* R){
    AB[0][0] = 0; AB[0][1] = 2; AB[0][2] = 3; AB[0][3] = 4; AB[0][4] = 49;
    AB[1][0] = 1; AB[1][1] = 0; AB[1][2] = 3; AB[1][3] = 4; AB[1][4] = 45;
    AB[2][0] = 1; AB[2][1] = 2; AB[2][2] = 0; AB[2][3] = 4; AB[2][4] = 36;
    AB[3][0] = 1; AB[3][1] = 2; AB[3][2] = 3; AB[3][3] = 0; AB[3][4] = 23;

    R[0] = 2; R[1] = 3; R[2] = 5; R[3] = 7;
}
void fill_matrix2(float** AB, float* R){
    AB[0][0] = 10; AB[0][1] = -7;     AB[0][2] = 0; AB[0][3] = 7;
    AB[1][0] = 0;  AB[1][1] = -0.001; AB[1][2] = 6; AB[1][3] = 6.001;
    AB[2][0] = 0;  AB[2][1] = 2.5;    AB[2][2] = 5; AB[2][3] = 2.5;

    R[0] = 0; R[1] = -1; R[2] = 1;
}
void fill_matrix3(float** AB, float* R){
    AB[2][0] = 6;  AB[2][1] = -2;  AB[2][2] = 2; AB[2][3] = 4;   AB[2][4] = 16;
    AB[3][0] = 12; AB[3][1] = -8;  AB[3][2] = 6; AB[3][3] = 10;  AB[3][4] = 26;
    AB[0][0] = 3;  AB[0][1] = -13; AB[0][2] = 9; AB[0][3] = 3;   AB[0][4] = -19;
    AB[1][0] = -6; AB[1][1] = 4;   AB[1][2] = 1; AB[1][3] = -18; AB[1][4] = -34;

    R[0] = 3; R[1] = 1; R[2] = -2; R[3] = 1;
}
