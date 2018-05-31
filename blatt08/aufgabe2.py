from pylab import *


def polynom_normal(x, y, n):
    a = empty([len(x), n+1])
    for i in range(n+1):
        a[:,i] = x**i

    return linalg.solve(a.T.dot(a), a.T.dot(y))


def polynom_qr(x, y, n):
    a = empty([len(x), n+1])
    for i in range(n+1):
        a[:,i] = x**i

    q, r = qr(a)
    return inv(r).dot(q.T.dot(y))


epsilon = 1.0e-4
x = array([0, epsilon, 2*epsilon, 3*epsilon, 1-3*epsilon, 1-2*epsilon, 1-epsilon, 1])
y = x - x**3
n = 5
print("Normalgleichung: ")
print(polynom_normal(x, y, n))
print("QR-Zerlegung: ")
print(polynom_qr(x, y, n))