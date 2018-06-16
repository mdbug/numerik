from sympy import *
from numpy import array, inf, linspace, meshgrid, log, copy
from numpy.linalg import norm
import matplotlib.pyplot as plt


def f(x):
    return (sin(x[0]) - x[1]) ** 2 + (exp(-x[1]) - x[0]) ** 2


def gradient(func, n):
    g = [None] * n
    x = symarray('x', n)

    for i in range(n):
        g[i] = diff(func(x), x[i])

    grad = lambdify(x, g, modules="numpy")

    def gg(v):
        return array(grad(*v))

    return gg


def plot_sd(f, p):
    # Argumentwerte als 1D Arrays erzeugen
    x_1d = linspace(min(p[:,0:1]) - 1, max(p[:,0:1]) + 1, 400)
    y_1d = linspace(min(p[:,1:]) - 1, max(p[:,1:]) + 1, 400)
    x_2d, y_2d = meshgrid(x_1d, y_1d)
    z_2d = f(x_2d, y_2d)
    # plt.figure()
    # fig, ax = plt.subplots()
    plt.pcolormesh(x_2d, y_2d, log(z_2d))
    plt.contour(x_2d, y_2d, log(z_2d), 10, colors="k", linewidths=0.5)
    plt.gca().set_aspect("equal") # x- und y-Skala im gleichen Maßstaab
    plt.plot(p[:,0:1], p[:,1:], 'rx-')

    plt.show()


def sd_armijo(func, x0, alpha, rho, tau, plot=true):
    n = len(x0)
    df = gradient(func, n)
    v = symarray('v', n)
    f = lambdify(v, func(v), "numpy")

    x = x0
    points = [copy(x0)]
    count = 0
    while True:
        p = -df(x)
        if norm(p, inf) < 1e-4:
            break

        if f(*(x + alpha * p)) > f(*x) + rho * p.dot(-p) * alpha:
            alpha *= tau

        x += alpha * p
        count += 1
        if plot:
            points.append(copy(x))

    if plot:
        plot_sd(f, array(points))

    print("Iterationen: %i" % count)
    return x


for startpukt in [(5,2), (6,2), (-1,-1), (-2,-2)]:
    x0 = array(startpukt, dtype=float)
    print("Startpunkt: (%f, %f)" % tuple(x0))
    xn = sd_armijo(f, x0, 1.0, 0.5, 0.5)
    print("Endpunkt: (%f, %f)" % tuple(xn))
    print("")
