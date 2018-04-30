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

    print("-"*50)


