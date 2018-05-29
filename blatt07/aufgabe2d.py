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


def x(t):
    return 1/(1+t) * cos(3*pi*t)


def y(t):
    return 1/(1+t) * sin(3*pi*t)


for m in [6, 7, 8]:
    t = linspace(0, 1, m)
    x_t = x(t)
    y_t = y(t)

    tt = linspace(0, 1, 100)

    d1 = divdiff(t, x_t)
    d2 = divdiff(t, y_t)
    p_x = newton_polynom(tt, t, d1)
    p_y = newton_polynom(tt, t, d2)
    plot(x_t, y_t, 'ro')
    plot(x(tt), y(tt))
    plot(p_x, p_y)
    show()
