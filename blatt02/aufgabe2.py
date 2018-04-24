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


a = array([
    [0, 0, 0, 1],
    [2, 1, 2, 0],
    [4, 4, 0, 0],
    [2, 3, 1, 0]], dtype=float)

u = array([0, 1, 2, 3], dtype=float)
v = array([0, 0, 0, 1], dtype=float)

a_dach = a + u.dot(v)
b_dach = array([3, 5, 4, 5], dtype=float)

lu, p = zerlegung(a)
pu = permutation(p, u)
y = vorwaerts(lu, pu)
z = rueckwaerts(lu, y)

# Ueberpruefung ob A_dach regulaer ist
assert 1 + v.dot(z) != 0

alpha = 1 / (1 + v.dot(z))

b_dach = permutation(p, b_dach)
y_dach = vorwaerts(lu, b_dach)
z_dach = rueckwaerts(lu, y_dach)

x_dach = z_dach - alpha*(v.dot(z_dach))*z

print(x_dach)
