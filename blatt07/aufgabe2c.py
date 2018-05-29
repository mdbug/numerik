from pylab import *


def newton_polynom(x, xx, d):
    n = len(xx) - 1
    q = d[n, 0]
    for k in range(1, n + 1):
        q = d[n - k, 0] + (x - xx[n - k]) * q

    return q


def divdiff(x, y):
    size = len(x)
    d = zeros([size, size])

    for n in range(0, size):
        d[n, n] = y[n]
        for m in reversed(range(0, n)):
            d[n, m] = (d[n, m + 1] - d[n - 1, m]) / (x[n] - x[m])

    return d


def f(x):
    return 1 / (1 + x**2)


for m in [7, 9, 11]:
    x = zeros(m)
    for i in range(m):
        x[i] = -5 * cos(pi * (2*i+1)/(2*m))

    xx = linspace(-5, 5, 100)
    y = f(x)
    d = divdiff(x, y)
    plot(x, y, 'ro')
    plot(xx, f(xx))
    plot(xx, newton_polynom(xx, x, d))
    show()

