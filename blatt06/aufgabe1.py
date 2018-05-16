from pylab import *


def newton(f, ff, x0):
    x = x0
    fx = f(x)
    while abs(fx) > 1.1e-14:
        x -= fx/ff(x)
        fx = f(x)

    return x


def sekanten(f, x0, x1):
    fx1 = f(x0)
    x = x1
    while abs(fx1) > 1.1e-14:
        fx0 = fx1
        fx1 = f(x1)
        x = x1 - ((x1 - x0) / (fx1 - fx0)) * fx1
        x0 = x1
        x1 = x

    return x


def b(x):
    return 9.8606 / (1 + 1.1085e25*exp(-0.029*x)) - 9


def bb(x):
    return 9.8606*0.029*1.1085e25*exp(0.029*x)/((-1.1085e25 - exp(0.029*x))**2)


print(sekanten(b, 1961.0, 2000.0))
print(newton(b, bb, 1961))
