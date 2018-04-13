import numpy as np


def integral_trapez(n, a, b, f):
    h = (b - a)/n
    v = np.arange(1, n)
    v = a + h*v
    v = f(v)
    return h/2 * (f(a) + 2*np.sum(v) + f(b))


def f(x):
    return 1/(x**2)


def g(x):
    return np.log(x)


n = 1000

a = 0.1
b = 10
int_f_trapez = integral_trapez(n, a, b, f)
int_f_exakt = 9.9
print("n = %d; int_a^b f:" % n)
print("trapez = %f, exakt = %f" % (int_f_trapez, int_f_exakt))

a = 1
b = 2
int_g_trapez = integral_trapez(n, a, b, g)
int_g_exakt = np.log(4) - 1
print("n = %d; int_a^b g:" % n)
print("trapez = %f, exakt = %f" % (int_g_trapez, int_g_exakt))
