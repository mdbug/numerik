from pylab import *


def newton(f, jf, v0):
    v = copy(v0)
    fv = f(v)
    while norm(fv) > 1.1e-14:
        v -= jf(v).dot(fv)
        fv = f(v)

    return v


def f(v: array) -> array:
    return array([sin(v[0]) - v[1], exp(-v[1]) - v[0]])


def jf(v: array) -> array:
    """ Jocobi Matrix zur Funktion f """
    return array([[cos(v[0]), -1],
                 [-1, -exp(-v[1])]])


print(newton(f, jf, array([0., 0.])))