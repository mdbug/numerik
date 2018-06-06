from pylab import *


def trapez(f, a, b, m):
    l = (b-a)/m
    ret = f(a) + f(b)
    for i in range(1, m):
        ret += 2 * f(a + i*l)

    return l/2*ret


def simpson(f, a, b, m):
    l = (b-a)/m
    ret = f(a) + f(b)
    for i in range(1, m):
        if i%2 == 0:
            ret += 2 * f(a + i*l)
        else:
            ret += 4 * f(a + i*l)

    return l/3*ret


def g_2(f, a, b):
    x = array([-sqrt(3/5), 0, sqrt(3/5)])
    beta = array([5/9, 8/9, 5/9])
    x_tilde = (b - a)/2 * x + (b - a)/2
    beta_tilde = (b - a)/2 * beta
    return beta_tilde.dot(f(x_tilde))


def f(x):
    return 1/(1+x**2)


a = 0
b = 1
print("Trapez:")
t = trapez(f, a, b, 8)
print(t)
print("Absoluter Fehler:")
print(fabs(t - pi/4))
print("====================")
print("Simpson:")
s = simpson(f, a, b, 4)
print(s)
print("Absoluter Fehler:")
print(fabs(s - pi/4))
print("====================")
print("Gauß G_2:")
g = g_2(f, a, b)
print(g)
print("Absoluter Fehler:")
print(fabs(g - pi/4))

