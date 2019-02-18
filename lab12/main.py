import numpy as np
import matplotlib.pyplot as plt
import approx
from numpy import pi, cos, sin, sqrt

def f1(x):
    return x/(2+x**2)
first_range = (0, 4.5)

def f2(x):
    return 10 + x**2 * 0.5 - 10 * cos(2*x)
second_range = (-2*pi, 2*pi)

def f3(x):
    return -x * sin(sqrt(3*abs(x-1)))
third_range = (-50, 50)

data = [(f1, first_range), (f2, second_range), (f3, third_range)]

def count_error(f, samples):
    return sum((y - f(x))**2 for x,y in samples)

for f, d in data:
    for n in [10]:
        x = np.linspace(d[0], d[1], 100)
        y = list(map(f, x))
        sp = approx.CubicSpline(x, y, smooth=0.85)

        xs = np.linspace(x[0], x[-1], 100)
        ys = sp(xs)

        plt.plot(x, y, 'r--', xs, ys, '-')

        plt.xlabel('x values')
        plt.ylabel('f(x) values')
        plt.title('Spline approximmation with\n' + str(n) + ' nodes.')

        plt.show()
        #plt.savefig('pol_aprox' + str(f.__name__) + str(n) + 'smps' + str(m) + 'basef_s')

