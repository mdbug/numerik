from numpy import *


def sqrt_heron(n, x, y_0=1):
    y = y_0
    for i in range(n):
        y = 0.5 * (y + x / y)

    return y


x = 2
print("x   n   y_0   y_n      Abweichung")
for n in range(10):
    for y_0 in [1, 2]:
        y_n = sqrt_heron(n, x, y_0)
        y = sqrt(x)
        d = abs(y_n - y)
        print("%d   %d   %d     %f %f" % (x, n, y_0, y_n, d))
