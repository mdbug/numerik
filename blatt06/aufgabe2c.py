from pylab import *


def newton(f, ff, x0, alpha):
    x = x0
    while True:
        fx = f(x)
        x_neu = x - fx/ff(x)
        if alpha/(1-alpha)*abs(x_neu - x) < 1.0e-6:
            return x_neu

        x = x_neu


def f(x):
    return x + log(x) - 2


def ff(x):
    return 1 + 1/x


print(newton(f, ff, 1.0, 0.25))
