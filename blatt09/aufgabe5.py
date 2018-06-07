from scipy import integrate

from pylab import *


def trapez(f, a, b, m):
    l = (b-a)/m
    ret = f(a) + f(b)
    for i in range(1, m):
        ret += 2 * f(a + i*l)

    return l/2*ret


def romberg(f, a, b, q):
    p = empty((q + 1, q + 1), float)
    for i in range(q + 1):
        p[i, 0] = trapez(f, a, b, 2**i)
        for j in range(1, i + 1):
            p[i, j] = p[i, j - 1] + (p[i, j - 1] - p[i - 1, j - 1])/((2**(2*j))-1)

    return p[q, q]


def bulirsch(f, a, b, q):
    p = empty((q + 1, q + 1), float)
    for i in range(q + 1):
        if i < 3:
            m = i + 1
        else:
            m = 2*i - 2
        p[i, 0] = trapez(f, a, b, m)
        for j in range(1, i + 1):
            p[i, j] = p[i, j - 1] + (p[i, j - 1] - p[i - 1, j - 1])/((2**(2*j))-1)

    return p[q, q]


def f(x):
    return sin(pi*x**2)


print("Romberg")
for q in range(8):
    print("Ordnung %i" % (2*q + 2))
    print(romberg(f, -1, 1, q))

print("")
print("Bulirsch")
for q in range(8):
    print("Ordnung %i" % (2*q + 2))
    print(bulirsch(f, -1, 1, q))

print("")
print("scipy.integrate.quad:")
print(integrate.quad(f, -1, 1))

