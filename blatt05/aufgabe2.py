from numpy import array, copy, zeros, ones
from numpy.linalg import solve


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


def nachiteration(a: array, b: array, lu: array, x: array, n: int) -> array:
    x_k = copy(x)
    for k in range(n):
        r = b - a.dot(x_k)
        p = rueckwaerts(lu, vorwaerts(lu, r))
        x_k += p

    return x_k


# Aufgabe 2 a)
print("Aufgabe 2 a)")
for n in [40, 50, 60]:
    a = zeros([n, n], dtype=float)
    b = zeros(n, dtype=float)

    # a initialisieren
    for i in range(n):
        for j in range(n):
            if i == j or j == n-1:
                a[i, j] = 1
            elif i > j:
                a[i, j] = -1

    # b initialisieren
    for i in range(n-1):
        b[i] = 2 - i

    b[n-1] = 2 - n

    lu, p = zerlegung_pivot(a)
    pb = permutation(p, b)
    x = rueckwaerts(lu, vorwaerts(lu, pb))
    x_nach = nachiteration(a, pb, lu, x, 1)
    print("n = %d" % n)
    print("x")
    print(x)
    print("x mit Nachiteration")
    print(x_nach)
    print(" ")

print("Aufgabe 2 b)")
for n in range(10, 15):
    a = zeros([n, n], dtype=float)
    b = zeros(n, dtype=float)
    b[0] = 1

    # a initialisieren
    for i in range(n):
        for j in range(n):
            if j == i:
                a[i, j] = 1
            elif j < i:
                a[i, j] = i + j

    lu, p = zerlegung_pivot(a)
    pb = permutation(p, b)
    x = rueckwaerts(lu, vorwaerts(lu, pb))
    x_nach = nachiteration(a, pb, lu, x, 10)
    print("n = %d" % n)
    print("x")
    print(x)
    print("x mit Nachiterationen")
    print(x_nach)
    print(" ")
    print(solve(a, b))
