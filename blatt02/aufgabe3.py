from numpy import array, zeros, sqrt, copy, transpose


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
    """
    Erstellt eine Matrix aus zwei Vektoren,
    welche die Haupt und eine Nebendiagonale repraesentieren
    """
    mat = zeros([d1.shape[0], d1.shape[0]])
    for i in range(d1.shape[0]):
        mat[i, i] = d1[i]

    for i in range(d2.shape[0]):
        mat[i +1, i] = d2[i]

    return mat


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


for n in [100, 1000, 10000]:
    print("n = %d" % n)

    # Initialisiere b
    b = array([-1.0] * n)
    b *= 1 / (n + 1)**2

    # Initialisiere Diagonalen
    d1 = array([2.0] * n)
    d2 = array([-1.0] * (n - 1))

    l_d1, l_d2 = cholesky(d1, d2)
    l = diagonalmatrix(l_d1, l_d2)
    y = vorwaerts(l, b)
    x = rueckwaerts(transpose(l), y)

    print(x)
    print("-" * 20)

