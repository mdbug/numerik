from numpy import dot, abs, array, zeros, copy, inf
from numpy.linalg import cond


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


def rueckwaerts(lu: array, x: array) -> array:
    y = copy(x)
    for i in reversed(range(x.size)):
        for j in range(i + 1, x.size):
            y[i] -= lu[i, j] * y[j]

        y[i] /= lu[i, i]

    return y


def loesung(a: array, b: array) -> array:
    lu, p = zerlegung_pivot(a)
    pb = permutation(p, b)
    y = vorwaerts(lu, pb)
    return rueckwaerts(lu, y)


def prager_oettli(a: array, b: array, x: array) -> float:
    r = abs(b - dot(a, x))
    s = dot(abs(a), abs(x)) + abs(b)
    return max(r/s)


for delta in [10**-8, 10**-10, 10**-12]:
    a = array([[3, 2, 1], [2, 2*delta, 2*delta], [1, 2*delta, -delta]])
    b = array([3 + 3*delta, 6*delta, 2*delta])
    x = loesung(a, b)
    epsilon = prager_oettli(a, b, x)
    kondition = cond(a, inf)
    maschinengenauigkeit = 2**-53
    print("delta = %e" % delta)
    print("x = ", end='')
    print(x)
    print("epsilon = %e" % epsilon)
    print("kondition = %f" % kondition)
    print("kondition * maschinengenauigkeit = %e" % (kondition * maschinengenauigkeit))
    print("-"*25)

