#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double f(double x)
{
  return(x*x + 2*x + 1);
}

double f1(double x){
    return x;
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
  double xp = atof(argv[1]);
  double xk = atof(argv[2]);
  int N  = atoi(argv[3]);

  printf("\nMetoda Prostokątów\nWartosc calki wynosi: %lf\n", rectangle_method_integral(xp, xk, N));
  printf("\nMetoda Trapezów\nWartosc calki wynosi: %lf\n", trapeze_method_integral(xp, xk, N)); 
  printf("\nMetoda Simpsona\nWartosc calki wynosi: %lf\n", simpson_method_integral(xp, xk, N));
  
  //printf("\n%.5lf", rectangle_method_integral(xp, xk, N));
  //printf("\n%.5lf", trapeze_method_integral(xp, xk, N)); 
  //printf("\n%.5lf\n", simpson_method_integral(xp, xk, N));

  return 0;
} 