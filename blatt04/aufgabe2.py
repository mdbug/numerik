from numpy import array, ones, sign, inf, argmax, absolute, zeros, copy
from numpy.linalg import norm, cond


def swap_rows(matrix: array, row_1: int, row_2: int):
    tmp = copy(matrix[row_1, :])
    matrix[row_1, :] = matrix[row_2, :]
    matrix[row_2, :] = tmp


def zerlegung_pivot(a: array) -> (array, array):
    lu = array(a, copy=True)
    n = a.shape[0]
    p = zeros(n - 1, dtype=int)
    for i in range(n - 1):
        max_zeile = i
        for k in range(i + 1, n):
            if abs(lu[k, i]) > abs(lu[max_zeile, i]):
                max_zeile = k

        p[i] = max_zeile + 1
        if max_zeile != i:
            swap_rows(lu, i, max_zeile)

        for j in range(i + 1, n):
            assert lu[i, i] != 0
            lu[j, i] /= lu[i, i]
            for k in range(i + 1, n):
                lu[j, k] -= (lu[j, i] * lu[i, k])

    return lu, p


def permutation(p: array, x: array) -> array:
    y = copy(x)
    for i in range(y.size - 1):
        tmp = y[i]
        y[i] = y[p[i] - 1]
        y[p[i] - 1] = tmp

    return y


def vorwaerts(lu: array, x: array) -> array:
    y = copy(x)
    for i in range(x.size):
        for j in range(i):
            y[i] -= lu[i, j] * y[j]

    return y


def vorwaerts2(l: array, x: array) -> array:
    y = copy(x)
    for i in range(x.size):
        for j in range(i):
            y[i] -= l[i, j] * y[j]

        y[i] /= l[i, i]

    return y


def rueckwaerts(lu: array, x: array) -> array:
    y = copy(x)
    for i in reversed(range(x.size)):
        for j in range(i + 1, x.size):
            y[i] -= lu[i, j] * y[j]

        y[i] /= lu[i, i]

    return y


def permutation_rueckwaerts(p: array, w: array) -> array:
    z = copy(w)
    for i in reversed(range(z.size - 1)):
        tmp = z[i]
        z[i] = z[p[i] - 1]
        z[p[i] - 1] = tmp

    return z


def hager(lu: array, p: array) -> array:
    n = lu.shape[0]
    x = ones(n, dtype=float) / n

    for i in range(5):
        v = vorwaerts2(lu.T, x)

        l = copy(lu)
        for i in range(l.shape[0]):
            l[i, i] = 1.0

        w = rueckwaerts(l.T, v)
        y = permutation_rueckwaerts(p, w)
        xi = sign(y)

        v = vorwaerts(lu, permutation(p, xi))
        z = rueckwaerts(lu, v)
        if norm(z, inf) <= z.T.dot(x):
            break
        j = argmax(absolute(z))
        x = zeros(n)
        x[j] = 1

    return norm(y, 1)


for delta in [10 ** -8, 10 ** -10, 10 ** -12]:
    a = array([[3, 2, 1],
               [2, 2 * delta, 2 * delta],
               [1, 2 * delta, -delta]], dtype=float)

    lu, p = zerlegung_pivot(a)
    hag = hager(lu, p)
    norm_a = norm(a, inf)
    print("delta: %e" % delta)
    print("Hager: %f " % hag)
    print("Zeilensummennorm A: %f" % norm_a)
    print("Kondition A: %f " % (norm_a * hag))
    print("numpy.linalg.cond(A): %f " % cond(a, inf))
    print("")
