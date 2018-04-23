from numpy import array, zeros, sqrt, copy


def cholesky(d1: array, d2: array) -> (array, array):
    l_d1 = copy(d1)
    l_d2 = copy(d2)

    for i in range(l_d1.size):
        l_d1[i] = sqrt(l_d1[i])
        if i < l_d1.size - 1:
            l_d2[i] /= l_d1[i]
            l_d1[i + 1] -= l_d2[i]**2

    return l_d1, l_d2


def diagonalmatrix(d1: array, d2: array) -> array:
    mat = zeros([d1.shape[0], d1.shape[0]])
    for i in range(d1.shape[0]):
        mat[i, i] = d1[i]

    for i in range(d2.shape[0]):
        mat[i +1, i] = d2[i]

    return mat


for n in [100, 1000, 10000]:
    print("n = %d" % n)

    # Initialisiere Diagonalen
    d1 = array([2.0] * n)
    d2 = array([-1.0] * (n - 1))

    l_d1, l_d2 = cholesky(d1, d2)
    print("Hauptdiagonale")
    print(l_d1)
    print("Nebendiagonale 1")
    print(l_d2)
    print("-" * 20)

