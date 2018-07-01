from pylab import *
from numpy.ma import empty, cos
from scipy.constants import pi


def poly(x, y, grad):
    a = empty([len(x), grad+1])
    for k in range(grad+1):
        a[:, k] = x**k

    return linalg.solve(a.T.dot(a), a.T.dot(y))


def g(p, x):
    y = 0.0
    for k in range(len(p)):
        y += p[k]*x**k

    return y


x = empty(6)
y = empty(6)
for i in range(6):
    x[i] = 2*i/5 - 1
    y[i] = cos(pi*x[i])

p = poly(x, y, 2)
print(p)

xx = linspace(min(x), max(x))
plot(xx, g(p, xx), label = 'g(x)')
plot(x, y, 'o')
title('Ausgleichspolynom')
xlabel('x')
ylabel('y')
grid(True, linestyle='--')
legend()
show()

