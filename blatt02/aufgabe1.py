from numpy import array, zeros, copy


def swap_rows(matrix, row_1, row_2):
    tmp = copy(matrix[row_1, :])
    matrix[row_1, :] = matrix[row_2, :]
    matrix[row_2, :] = tmp


def zerlegung(a):
    lu = array(a, copy=True)
    n = a.shape[0]
    p = zeros(n-1)
    for i in range(n):
        if lu[i, i] == 0:
            k = i
            while lu[k, i] == 0:
                k += 1

            p[i] = k+1
            # Tausche Zeile i mit k
            swap_rows(lu, i, k)

        for j in range(i + 1, n):
            lu[j, i] /= lu[i, i]
            for k in range(i + 1, n):
                lu[j, k] -= (lu[j, i]*lu[i, k])

    return lu, p


def permutation(p, x):
    for i in range(x.size):
        tmp = x[i]
        x[i] = x[p[i] - 1]
        x[p[i] - 1] = tmp


def vorwaerts(lu, x):
    # TODO


def rueckwaerts(lu, x):
    # TODO


a = array([[0, 0, 0, 1], [2, 1, 2, 0], [4, 4, 0, 0], [2, 3, 1, 0]])
lu, p = zerlegung(a)

print(lu)
print(p)
