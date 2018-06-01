from pylab import *
from scipy import randn


n = 100
gamma = 0.05

# Initialisiere x
x = zeros(n)
x[45:56] = 1
x[60:66] = 0.5

# Initialisiere a
c = 1/(gamma*sqrt(2*pi))
a = empty([n, n])
for i in range(n):
    for j in range(n):
        a[i, j] = c/n * exp(-((i-j)/(sqrt(2)*n*gamma))**2)

db = 1.0e-6 * randn(n)
b = a.dot(x)

y = linalg.pinv(a).dot(b+db)
plot(x, label="Original-Signal")
plot(b + db, label="Bild-Signal")
plot(y, label="Rekonstruktion über Pseudoinverse")
legend()
show()

u, s, vh = linalg.svd(a)
for k in range(9):
    alpha = 10.0**-k
    sp = copy(s)
    for i in range(len(s)):
        if s[0]/s[i] >= 1/alpha:
            sp[0:i] = 1 / s[0:i]
            sp[i:len(sp)] = 0
            break

    z = vh.T.dot(diag(sp)).dot(u.T).dot(b+db)
    plot(x, label="Original-Signal")
    plot(b + db, label="Bild-Signal")
    plot(z, label='TSVD alpha=10^%d' % k)
    legend()
    show()

