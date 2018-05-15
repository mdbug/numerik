from pylab import *
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import axes3d

from scipy.sparse import spdiags


def Ablock(m):
    # Blockmatrix fuerr 2d-Laplace
    n = m*m

    e = ones(n)
    l = ones(n)
    l[m-1::m] = 0.0
    r = ones(n)
    r[::m] = 0.0

    A = spdiags([-e, -l, 4.0*e, -r, -e], [-m, -1, 0, 1, m], n, n, format='csr')

    return A.toarray()


def ew_exakt(m):
    # exakte Eigenwerte fuer 2d-Laplace Blockmatrix, aufsteigend sortiert
    ew1d = 2.0 * (1.0 - cos( (arange(m) + 1.0) * pi / (m + 1.0) ))

    ew = (c_[ew1d] + ew1d).flatten()
    ew.sort()

    return ew


def plotev(xk):
    # Eigenvektoren fuer 2d-Laplace Blockmatrix graphisch darstellen
    n = len(xk)
    m = int(sqrt(n))

    h = linspace(0, 1, m)
    yy,xx = meshgrid(h,h)

    fig = figure()
    ax  = fig.gca(projection='3d')

    surf = ax.plot_surface(xx, yy, xk.reshape(m,m), cmap = cm.jet, rstride = 1, cstride = 1)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("Hoehe")
    colorbar(surf)
    plt.show()


def animev(xk):
    # Eigenvektoren fuer 2d-Laplace Blockmatrix animieren

    n = len(xk)
    m = int(sqrt(n))

    h = linspace(0, 1, m)
    yy,xx = meshgrid(h,h)
    zz    = xk.reshape(m,m)

    zmax = 1.1 * abs(zz).max()

    def a(nf = 100, inter = 100, rep = False):
        fig = figure()
        #ax  = fig.gca(projection='3d')
        ax = axes3d.Axes3D(fig)
        ax.set_axis_off()
        ax.grid(False)

        surf = ax.plot_surface(xx, yy, zz, cmap = cm.jet, rstride = 1, cstride = 1)
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.set_zlabel("Hoehe")
        #colorbar(surf)

        ax.set_zlim(-zmax, zmax)

        def update(i, ax, fig):
            ax.cla()
            phi = i * 180.0 / 3.14159 / nf
            zzi = cos(phi) * zz
            wframe = ax.plot_surface(xx, yy, zzi, cmap = cm.jet, rstride = 1, cstride = 1)
            ax.set_zlim(-zmax, zmax)
            ax.set_axis_off()
            ax.grid(False)
            return wframe,

        return animation.FuncAnimation(fig, update,
                                       frames = range(nf),
                                       fargs  = (ax, fig), interval = inter, repeat = rep)

    return a


def jacobi(aa):
    a = copy(aa)

    def signum(t):
        if t == 0:
            return 1
        else:
            return sign(t)

    def max_off_diag(a):
        i, j = unravel_index(argmax(abs(triu(a, 1))), a.shape)
        return a[i, j], i, j

    n = len(a)
    q_k = eye(n)
    while sum(square(triu(a, 1)))*2 > 10**-3:
        a_max, i, j = max_off_diag(a)
        alpha = (a[j, j] - a[i, i]) / (2 * a[i, j])
        c = sqrt(0.5 + 0.5 * sqrt((alpha**2) / (1 + alpha**2)))
        s = signum(alpha)/(2*c*sqrt(1+alpha**2))
        q = eye(n)
        q[i, i] = c
        q[j, j] = c
        q[i, j] = s
        q[j, i] = -s
        q_k = q_k.dot(q)
        a = q.T.dot(a).dot(q)

    ew = copy(a.diagonal())
    permutation = argsort(ew)

    return ew[permutation], q_k[:,permutation]


m = 10
a = Ablock(m)
ew_j, q_k = jacobi(a)
ew_e = ew_exakt(m)

print("Eigenwerte Jacobi")
print(ew_j)
print("Eigenwerte Exakt")
print(ew_e)

# for j in range(4):
#     v = q_k[:,j]
#     plotev(v)
#     animev(v)()
