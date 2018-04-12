from math import sqrt


def pq_formel(p, q):
    return -p + sqrt(p**2 + q)


def alt_formel(p, q):
    x1 = -p - sqrt(p**2 + q)
    return -q / x1


q = 1
for p in [10**2, 10**4, 10**6, 10**7, 10**8]:
    x2_pq = pq_formel(p, q)
    x2_alt = alt_formel(p, q)
    print("p = %e, q = %2f, x2_pq = %e, x2_alt = %e" % (p, q, x2_pq, x2_alt))

# alt_formel ist fuer grosse p genauer
