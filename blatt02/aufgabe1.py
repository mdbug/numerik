from numpy import array, zeros, copy


def swap_rows(a, i, j):
    tmp = copy(a[i, :])
    a[i, :] = a[j, :]
    a[j, :] = tmp


def zerlegung(a):
    lu = array(a, copy=True)
    n = a.shape[0]
    p = zeros(n)
    for i in range(n-1):
        if lu[i, i] == 0:
            k = i
            while lu[k, i] == 0:
                k += 1

            # Tausche Zeile i mit k
            swap_rows(a, i, k)

            p[i] = lu[k, i]

        for j in range(i + 1, n):
            lu[j, i] /= lu[i, i]

    return lu, p


a = array([[0, 0, 0, 1], [2, 1, 2, 0], [4, 4, 0, 0], [2, 3, 1, 0]])
lu, p = zerlegung(a)

print(lu)
print(p)




