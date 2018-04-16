from numpy import *


def sqrt_taylor(n, x, y_0=1):
    x_0 = y_0**2
    y = y_0
    a = y
    factor = x/x_0 - 1
    for k in range(1, n + 1):
        a *= (3 / (2 * k) - 1) * factor
        y += a

    return y


def sqrt_heron(n, x, y_0=1):
    y = y_0
    for i in range(n + 1):
        y = 0.5 * (y + x / y)

    return y


x = 2
print("sqrt_taylor:")
print("x   n    y_0   y_n      Abweichung")
for y_0 in [1, 2]:
    for n in range(1, 11):
        y_n = sqrt_taylor(n, x, y_0)
        y = sqrt(x)
        d = abs(y_n - y)
        print("%d   %2d   %d     %f %f" % (x, n, y_0, y_n, d))


# Man benoetigt 9 Schritte fuer y_0 = 1 und 4 Schritte fuer y_0 = 2 bis die Abweichung kleiner als 0.005 ist


print("")
print("sqrt_heron:")
print("x   n    y_0   y_n      Abweichung")
for y_0 in [1, 2]:
    for n in range(1, 11):
        y_n = sqrt_heron(n, x, y_0)
        y = sqrt(x)
        d = abs(y_n - y)
        print("%d   %2d   %d     %f %f" % (x, n, y_0, y_n, d))


# Das Verfahren konvergiert in deutlich weniger Schritten
