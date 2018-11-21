#include <stdio.h>
#include <math.h>
#include <gsl/gsl_integration.h>

double f (double x, void * params) {
  double alpha = *(double *) params;
  double f = x*x + alpha*x + 1;
  return f;
}

double f1(double x){
    double f1 = x;
    return f1;
}
double f2(double x){
    return (2*(x*x));
}
double f3(double x){
    return (4*sin(x));
}
double f4(double x){
    return (exp(x));
}
double f5(double x){
    return (x * sin(x)*sin(x) + 2*cos(x));
}
double f6(double x){
    return cos((x+1)/(x*x + 0.04))*exp(x);
}

double rectangle_method_integral(double xp, double xk, int N){
  double s, dx;
  int i;

  s  = 0.0;
  dx = (xk - xp) / N;
  for(i = 1; i <= N; i++) s += f1(xp + i * dx);
  s *= dx;

  return s;
}

double trapeze_method_integral(double xp, double xk, int N){
  double s, dx, a, b;
  int i;

  s  = 0.0;
  dx = (xk - xp) / N;
  a = f1(xp);

  for(i = 1; i <= N; i++){ 
    b = f1(xp + i*dx);
    s += (a + b);
    a = b;
  }

  return s*dx*0.5;
}

double simpson_method_integral(double xp, double xk, int N){
  double s, st, x, dx;
  int i;

  s  = 0; st = 0;
  dx = (xk - xp) / N;
  for(i = 1; i <= N; i++)
  {
    x = xp + i * dx;
    st += f1(x - dx / 2);
    if(i < N) s += f1(x);
  }
  s = dx / 6 * (f1(xp) + f1(xk) + 2 * s + 4 * st);

  return s;
}

int main (int argc, char *argv[])
{
  double xp = atof(argv[3]);
  double xk = atof(argv[4]);
  int N  = atoi(argv[5]);
  gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);

  double result, error;
  double alpha = 2.0;

  gsl_function F;
  F.function = &f1;
  F.params = &alpha;

    gsl_integration_qags (&F, xp, xk, 0, 1e-7, 1000, w, &result, &error);
    printf("%lf\n", fabs(rectangle_method_integral(xp,xk,N) - result));
    printf("%lf\n", fabs(trapeze_method_integral(xp,xk,N) - result));
    printf("%lf\n", fabs(simpson_method_integral(xp,xk,N) - result));

  gsl_integration_workspace_free (w);

  return 0;
}