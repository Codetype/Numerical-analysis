#include <stdio.h>
#include <math.h>
#include  <stdlib.h>

#define N 5
double lroots[N];
double weight[N];
double lcoef[N + 1][N + 1] = {{0}};

double f0(double x){
	return exp(x);
}

double f1(double x){
    return 3*pow(x, 3)-1;
}
double f2(double x){
    return 2 * pow(x, 2);
}
double f3(double x){
    return 4*sin(x);
}

void lege_coef()
{
	int n, i;
	lcoef[0][0] = lcoef[1][1] = 1;
	for (n = 2; n <= N; n++) {
		lcoef[n][0] = -(n - 1) * lcoef[n - 2][0] / n;
		for (i = 1; i <= n; i++)
			lcoef[n][i] = ((2 * n - 1) * lcoef[n - 1][i - 1] - (n - 1) * lcoef[n - 2][i] ) / n;
	}
}
 
double lege_eval(int n, double x)
{
	int i;
	double s = lcoef[n][n];
	for (i = n; i; i--)
		s = s * x + lcoef[n][i - 1];
	return s;
}
 
double lege_diff(int n, double x)
{
	return n * (x * lege_eval(n, x) - lege_eval(n - 1, x)) / (x * x - 1);
}
 
void lege_roots()
{
	int i;
	double x, x1;
	for (i = 1; i <= N; i++) {
		x = cos(M_PI * (i - .25) / (N + .5));
		do {
			x1 = x;
			x -= lege_eval(N, x) / lege_diff(N, x);
		} while ( fdim( x, x1) > 2e-16 );

		lroots[i - 1] = x;
 
		x1 = lege_diff(N, x);
		weight[i - 1] = 2 / ((1 - x * x) * x1 * x1);
	}
}

double lege_inte(double (*f)(double), double a, double b)
{
	double c1 = (b - a) / 2, c2 = (b + a) / 2, sum = 0;
	int i;
	for (i = 0; i < N; i++)
		sum += weight[i] * f(c1 * lroots[i] + c2);
	return c1 * sum;
}


double rectangle_method_integral(double (*f)(double), double xp, double xk, int Nt){
  double s, dx;
  int i;

  s  = 0.0;
  dx = (xk - xp) / Nt;
  for(i = 1; i <= Nt; i++) s += f(xp + i * dx);
  s *= dx;

  return s;
}

double trapeze_method_integral(double (*f)(double), double xp, double xk, int Nt){
  double s, dx, a, b;
  int i;

  s  = 0.0;
  dx = (xk - xp) / Nt;
  a = f(xp);

  for(i = 1; i <= Nt; i++){ 
    b = f(xp + i*dx);
    s += (a + b);
    a = b;
  }

  return s*dx*0.5;
}

double simpson_method_integral(double (*f)(double), double xp, double xk, int Nt){
  double s, st, x, dx;
  int i;

  s  = 0; st = 0;
  dx = (xk - xp) / Nt;
  for(i = 1; i <= Nt; i++)
  {
    x = xp + i * dx;
    st += f(x - dx / 2);
    if(i < Nt) s += f(x);
  }
  s = dx / 6 * (f(xp) + f(xk) + 2 * s + 4 * st);

  return s;
}

  

int main(int argc, char* argv[])
{
	int i;
    double xp = atof(argv[1]);
    double xk = atof(argv[2]);
    int Nt  = atoi(argv[3]);

	lege_coef();
	lege_roots();

    printf("\nMetoda Gaussa-Legendre'a\nWartosc calki wynosi: %lf\n", lege_inte(&f0, xp, xk));
    printf("\nMetoda Prostokątów\nWartosc calki wynosi: %lf\n", rectangle_method_integral(&f0, xp, xk, Nt));
    printf("\nMetoda Trapezów\nWartosc calki wynosi: %lf\n", trapeze_method_integral(&f0, xp, xk, Nt)); 
    printf("\nMetoda Simpsona\nWartosc calki wynosi: %lf\n", simpson_method_integral(&f0, xp, xk, Nt));

	return 0;
}