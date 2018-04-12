from math import log


def integral_untersumme(n, a, b, f):
    summe = 0
    h = (b - a)/n
    for i in range(0, n):
        summe += f(a + i*h)

    return h*summe


def integral_obersumme(n, a, b, f):
    summe = 0
    h = (b - a)/n
    for i in range(1, n + 1):
        summe += f(a + i*h)

    return h*summe


def integral_trapez(n, a, b, f):
    summe = 0
    h = (b - a)/n
    for i in range(1, n):
        summe += f(a + i*h)

    return h/2 * (f(a) + 2*summe + f(b))


def f(x):
    return 1/(x**2)


def g(x):
    return log(x)


n = 1000
a = 0.1
b = 10
int_f_untersumme = integral_untersumme(n, a, b, f)
int_f_obersumme = integral_obersumme(n, a, b, f)
int_f_trapez = integral_trapez(n, a, b, f)
int_f_exakt = 9.9
print("n = %d; int_a^b f:" % (n))
print("untersumme = %f, obersumme = %f, trapez = %f, exakt = %f" % (int_f_untersumme, int_f_obersumme, int_f_trapez, int_f_exakt))

a = 1
b = 2
int_g_untersumme = integral_untersumme(n, a, b, g)
int_g_obersumme = integral_obersumme(n, a, b, g)
int_g_trapez = integral_trapez(n, a, b, g)
int_g_exakt = log(4) - 1
print("n = %d; int_a^b g:" % (n))
print("untersumme = %f, obersumme = %f, trapez = %f, exakt = %f" % (int_g_untersumme, int_g_obersumme, int_g_trapez, int_g_exakt))
