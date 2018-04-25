from numpy import array, zeros, copy


def swap_rows(matrix: array, row_1: int, row_2: int):
    tmp = copy(matrix[row_1, :])
    matrix[row_1, :] = matrix[row_2, :]
    matrix[row_2, :] = tmp


def zerlegung(a: array) -> (array, array):
    lu = array(a, copy=True)
    n = a.shape[0]
    p = zeros(n-1, dtype=int)
    for i in range(n):
        if lu[i, i] == 0:
            k = i
            while lu[k, i] == 0:
                k += 1

            p[i] = k+1
            swap_rows(lu, i, k)

        for j in range(i + 1, n):
            lu[j, i] /= lu[i, i]
            for k in range(i + 1, n):
                lu[j, k] -= (lu[j, i]*lu[i, k])

    return lu, p


def zerlegung_pivot(a: array) -> (array, array):
    lu = array(a, copy=True)
    n = a.shape[0]
    p = zeros(n-1, dtype=int)
    for i in range(n - 1):
        mx = i
        # find pivot
        for k in range(i + 1, n):
            if abs(lu[k, i]) > abs(lu[mx, i]):
                mx = k

        p[i] = mx + 1
        if mx != i:
            swap_rows(lu, i, mx)

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


# a = array([
#     [0., 0, 0, 1],
#     [2, 1, 2, 0],
#     [4, 4, 0, 0],
#     [2, 3, 1, 0]])
#
# b1 = array([3., 5, 4, 5])
# b2 = array([4., 10, 12, 11])
#
# lu, p = zerlegung_pivot(a)
#
# pb1 = permutation(p, b1)
# pb2 = permutation(p, b2)
#
# y1 = vorwaerts(lu, pb1)
# y2 = vorwaerts(lu, pb2)
#
# x1 = rueckwaerts(lu, y1)
# x2 = rueckwaerts(lu, y2)
#
# print(x1)
# print(x2)
# print("")
#

beta = 10.0
for n in [10, 15, 20]:
    print("n = %d" % n)

    # Initialisiere a
    a = zeros([n, n])
    for i in range(n - 1):
        a[i, i] = 1.0
        a[i + 1, i] = -beta

    a[0, n - 1] = beta

    # Initialisiere b
    b = array([1 - beta] * n)
    b[0] = 1 + beta
    b[n - 1] = -beta

    lu, p = zerlegung(a)
    pb = permutation(p, b)
    y = vorwaerts(lu, pb)
    x = rueckwaerts(lu, y)

    print("Loesung ohne Spaltenpivot")
    print(x)

    lu, p = zerlegung_pivot(a)
    pb = permutation(p, b)
    y = vorwaerts(lu, pb)
    x = rueckwaerts(lu, y)

    print("Loesung mit Spaltenpivot")
    print(x)


