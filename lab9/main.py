import numpy as np
from math import sqrt, sin
from functools import reduce
from operator import mul
import matplotlib.pyplot as plt

def F1(x): return sqrt(x)
def F2(x): return sin(x)
def F3(x): return x**3 + 2*x

FUNCTIONS = [F1,F2,F3]
RANGE = (0,10)
SAMPLES = [3,4,5,8]

def tabularize(f, range, samples):
    return [(x, f(x)) for x in np.linspace(range[0], range[1], samples)]

def lagrange(samples):
    return lambda x: sum(
        samples[i][1] * reduce(mul,(
            (x - samples[j][0]) / (samples[i][0] - samples[j][0])
            for j in range(len(samples))
            if j != i
        ), 1)
        for i in range(len(samples))
    )

def main():
    tables = dict()

    for f in FUNCTIONS:
        tables[f.__name__] = dict()
        for s in SAMPLES:
            tables[f.__name__][s] = tabularize(f, RANGE, s)


    x = np.linspace(RANGE[0], RANGE[1], 100)

    for f in FUNCTIONS:
        for s in SAMPLES:
            plt.figure()
            plt.plot(*zip(*tables[f.__name__][s]), 'ro')
            plt.plot(x, list(map(f,x)), 'r--')
            plt.plot(x,list(map(lagrange(tables[f.__name__][s]), x)), 'c')

main()