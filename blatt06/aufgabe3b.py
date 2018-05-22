from pylab import *


def newton(f, ff, x0):
    x = x0
    for i in range(4):
        x -= f(x)/ff(x)
        print("x_%d = %f" % (i+1, x))

    return x


def newton_mod1(f, ff, x0, q):
    x = x0
    for i in range(4):
        x -= q*f(x)/ff(x)
        print("x_%d = %f" % (i+1, x))

    return x


def newton_mod2(f, ff, x0):
    x = x0
    for i in range(4):
        fx = f(x)
        ffx = ff(x)
        x -= (fx*ffx)/(ffx**2 - fx*ffx)
        print("x_%d = %f" % (i+1, x))

    return x


def f(x):
    return arctan(x) - x


def ff(x):
    return 1/(x**2 + 1) - 1


print("Newton:")
newton(f, ff, 1.0)
print("Newton Modifikation 1:")
newton_mod1(f, ff, 1.0, 3)
print("Newton Modifikation 2:")
newton_mod2(f, ff, 1.0)
