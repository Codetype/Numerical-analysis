from scipy import random
import numpy as np
from math import sqrt

def monte_carlo_basic(f, x0, xN, N):
    xrand = np.zeros(N)

    for i in range(len(xrand)):
        xrand[i] = random.uniform(x0, xN)

    integral = 0.0

    for i in range(N):
        integral += f(xrand[i])

    return (xN-x0)/float(N)*integral


def monte_carlo_integral(f, x1, x2, y1, y2, N):
    counter = 0
    for i in range(N):
        x = random.uniform(x1, x2)
        y = random.uniform(y1, y2)

        if 0 < y < f(x):
            counter += 1
        if f(x) < y < 0:
            counter -= 1

    result = (x2 - x1) * (y2 - y1) * (counter / N)
    if y1 > 0:
        result += y1 * (x2 - x1)
    if y2 < 0:
        result += y2 * (x2 - x1)

    return result

def monte_carlo_volume(f, x1, x2, y1, y2, z1, z2, N):
    counter = 0
    for i in range(N):
        x = random.uniform(x1, x2)
        y = random.uniform(y1, y2)
        z = random.uniform(z1, z2)

        if f(x, y, z):
            counter += 1

    return (x2 - x1) * (y2 - y1) * (z2 - z1) * (counter / N)


samples = [
        (lambda x: 1 / x**2,           1, 2 , 1/4, 1, 1/2),
        (lambda x: 1 / sqrt(x**5 + 8), 1, 6 , 0,   1/3, 0.43506),
        (lambda x: 1 / sqrt(x + 8),    1, 22, 0,   1/3, 2 * sqrt(30) - 6),
        (lambda x: x**2 + 2*x,         1, 22, 3,   528, 4032),
        (lambda x: sqrt(x**5),         1, 3,  1,   16, 54 / 7 * sqrt(3) - 2 / 7)
]

def v1(x, y, z):
    return x**2 + y**2 + z**2 <= 5**2

def v2(x, y, z):
    return (x**2 + y**2 <= 10**2) and (0 <= z <= 10 - sqrt(x**2 + y**2))

def v3(x, y, z):
    return (x**2 + y**2 + z**2 <= 10**2) and not ((x**2 + y**2 <= 3**2) and abs(z) <= 6/2)

volumes = [
    (v1,  -5,  5,  -5,  5,  -5,  5),
    (v2, -10, 10, -10, 10,   0, 10),
    (v3, -10, 10, -10, 10, -10, 10)
]

volumes_expected = [
    500 / 3 * np.pi,
    1000 / 3 * np.pi,
    3838 / 3 * np.pi
]

def make_integral_test():
    for i in range(0, 4):
        for n in range(200, 1200, 200):
            print(monte_carlo_basic(*samples[i][:3], n))
            print(monte_carlo_integral(*samples[i][:5], n))

def make_volume_test():
    for i in range(0, 2):
        for n in range(200, 1200, 200):
            print(monte_carlo_volume(*volumes[i], n))
            print(volumes_expected[i])

def main():
    make_integral_test()
    #make_volume_test()

main()
