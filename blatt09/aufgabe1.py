from math import fabs, pi


def trapez(f, a, b, m):
    l = (b-a)/m
    ret = f(a) + f(b)
    for i in range(1, m):
        ret += 2 * f(a + i*l)

    return l/2*ret


def simpson(f, a, b, m):
    l = (b-a)/m
    ret = f(a) + f(b)
    for i in range(1, m):
        if i%2 == 0:
            ret += 2 * f(a + i*l)
        else:
            ret += 4 * f(a + i*l)

    return l/3*ret


def f(x):
    return 1/(1+x**2)


print("Trapez:")
t = trapez(f, 0, 1, 8)
print(t)
print("Absoluter Fehler:")
print(fabs(t - pi/4))
print("====================")
print("Simpson")
s = simpson(f, 0, 1, 8)
print(s)
print("Absoluter Fehler:")
print(fabs(s - pi/4))


