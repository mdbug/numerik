from pylab import *


def newton_polynom(x, xx, yy):
    n = len(xx) - 1
    d = divdiff(xx, yy)
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


x = array([0, 1, 3], dtype=float)
y = array([3, 2, 6], dtype=float)
print(newton_polynom(2.0, x, y))
