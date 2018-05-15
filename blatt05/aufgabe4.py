from pylab import *


def jacobi(a):
    def signum(t):
        if t == 0:
            return 1
        else:
            return sign(t)

    def max_off_diag(a):
        l, k = unravel_index(argmax(abs(triu(a))), a.shape)
        return a[l, k], l, k

    def rotate(a, p, k, l):
        n = len(a)

    n = len(a)
    p = identity(n)

    while True:
        a_max, k, l = max_off_diag(a)
        rotate(a, p, k, l)
