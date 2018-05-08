from numpy import copy, array, zeros
from numpy.linalg import solve


def swap_rows(matrix: array, row_1: int, row_2: int):
    tmp = copy(matrix[row_1, :])
    matrix[row_1, :] = matrix[row_2, :]
    matrix[row_2, :] = tmp


def zerlegung_pivot(a: array) -> (array, array):
    lu = array(a, copy=True)
    n = a.shape[0]
    p = zeros(n-1, dtype=int)
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
                lu[j, k] -= (lu[j, i]*lu[i, k])

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
            y[i] -= lu[i, j]*y[j]

    return y


def vorwaerts2(lu: array, x: array) -> array:
    y = copy(x)
    for i in range(x.size):
        for j in range(i):
            y[i] -= lu[i, j]*y[j]

        y[i] /= lu[i,i]

    return y


def rueckwaerts(lu: array, x: array) -> array:
    y = copy(x)
    for i in reversed(range(x.size)):
        for j in range(i + 1, x.size):
            y[i] -= lu[i, j] * y[j]

        y[i] /= lu[i, i]

    return y


def berechne_v(lu: array, c: array) -> array:
    return vorwaerts2(lu.T, c)


def berechne_w(lu: array, v: array) -> array:
    l = copy(lu)
    for i in range(l.shape[0]):
        l[i, i] = 1.0

    return rueckwaerts(l.T, v)


def berechne_z(p: array, w: array) -> array:
    z = copy(w)
    for i in reversed(range(z.size - 1)):
        tmp = z[i]
        z[i] = z[p[i] - 1]
        z[p[i] - 1] = tmp

    return z


def loese(a: array, c: array) -> array:
    lu, p = zerlegung_pivot(a)
    v = berechne_v(lu, c)
    w = berechne_w(lu, v)
    return berechne_z(p, w)


a = array([[0, 0, 0, 1],
           [2, 1, 2, 0],
           [4, 4, 0, 0],
           [2, 3, 1, 0]], dtype=float)

c = array([152, 154, 56, 17], dtype=float)

print(loese(a, c))
