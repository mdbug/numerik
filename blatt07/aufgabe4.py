from pylab import *


def nat_kub_spline(x, y):
    n = len(x)
    h = empty(n - 1)
    for i in range(n-1):
        h[i] = x[i + 1] - x[i]

    a = zeros([n - 2, n - 2])
    for i in range(n-2):
        a[i][i] = 2 * (h[i] + h[i + 1])
        if i != 0:
            a[i][i - 1] = h[i]
        if i != n-3:
            a[i][i + 1] = h[i + 1]

    gamma = empty(n - 2)
    for i in range(n - 2):
        gamma[i] = 6 * ((y[i + 1] - y[i])/h[i] - (y[i] - y[i - 1])/h[i - 1])

    beta = zeros(n)
    beta[1:n-1] = linalg.solve(a, gamma)

    alpha = empty(n-1)
    for i in range(n-1):
        alpha[i] = (y[i + 1] - y[i])/h[i] - 1/3*beta[i]*h[i] - 1/6*beta[i + 1]*h[i]

    return h, beta, alpha


def nks_auswerten(xx, h, beta, alpha, x, y):
    n = len(x)
    i = 0
    while xx > x[i + 1] and i < n-2:
        i += 1

    return y[i] + alpha[i]*(xx - x[i]) + (beta[i]/2)*((xx-x[i])**2) + (beta[i + 1] - beta[i])/(6*h[i])*((xx - x[i])**3)


def f(x):
    return 1 / (1 + x**2)


for m in [7, 9, 11]:
    x = linspace(-5, 5, m)
    y = f(x)
    xx = linspace(-5, 5, 500)
    yy = empty(500)
    h, beta, alpha = nat_kub_spline(x, y)
    print("h:")
    print(h)
    print("beta:")
    print(beta)
    print("alpha:")
    print(alpha)
    for i in range(500):
        yy[i] = nks_auswerten(xx[i], h, beta, alpha, x, y)
    plot(x, y, 'ro')
    plot(xx, f(xx))
    plot(xx, yy)
    show()
