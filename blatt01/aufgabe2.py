from math import factorial


def exp1(n, x):
    result = 0
    for k in range(0, n+1):
        result += x**k / factorial(k)

    return result


def exp2(n, x):
    result = 0
    for k in range(0, n+1):
        result += x**(n-k) / factorial(n-k)

    return result


for n in [1, 10, 100]:
    for x in [1, 10, 100, 1000]:
        print("exp1(%d, %5d) = %e, exp2(%d, %5d) = %e" % (n, x, exp1(n, x), n, x, exp2(n, x)))
