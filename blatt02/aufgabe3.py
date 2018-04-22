from numpy import array, zeros, sqrt


def cholesky(a: array) -> array:
    d1 = zeros(a.size)
    d2 = zeros(a.size - 1)

    for i in range(a.size):
        d1 = a[i, i]

    for i in range(a.size-1):
        d2 = a[i+1, i]

    for i in range(d1):
        d1[i] = sqrt(i)
