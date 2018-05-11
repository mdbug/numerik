from pylab import ones, linspace, meshgrid, sqrt, figure, cm, colorbar, zeros, norm, copy, plt
from mpl_toolkits.mplot3d import Axes3D

from scipy.sparse import spdiags


def system(m):
    n = m*m

    e = ones(n)
    l = ones(n)
    l[m-1::m] = 0.0
    r = ones(n)
    r[::m] = 0.0

    A = spdiags([-e, -l, 4.0*e, -r, -e], [-m, -1, 0, 1, m], n, n, format='csr')

    b = -e / float(n)

    return A, b


def plotxk(xk):
    n = len(xk)
    m = int(sqrt(n))

    print(m)

    h = linspace(0, 1, m)
    yy,xx = meshgrid(h,h)

    fig = figure('xk, m = {0}'.format(m))
    ax = fig.gca(projection='3d')

    surf = ax.plot_surface(xx, yy, xk.reshape(m,m), cmap = cm.jet, rstride = 5, cstride = 5)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("Hoehe")
    colorbar(surf)
    plt.show()


def cg(A, b, x0):
    p = b - A.dot(x0)
    r = copy(p)
    norm_r0 = norm(r)
    x = copy(x0)

    while True:
        r_k_dot_r_k = r.dot(r)
        alpha = r_k_dot_r_k / p.dot(A.dot(p))
        x += alpha * p
        r -= alpha * A.dot(p)
        beta = r.dot(r) / p.dot(A.dot(p))
        p = r + beta * p
        if norm(r)/norm_r0 <= 1e-6:
            return x


for m in [50, 100, 200]:
    n = m**2
    A, b = system(m)
    x = cg(A, b, zeros(n))
    print(x)
    plotxk(x)

