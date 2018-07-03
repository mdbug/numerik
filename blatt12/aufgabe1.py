import numpy as np
from numpy.linalg import norm

def f(x):
    return (1 / 2) * norm(F(x))**2

def q_s(x, x_neu):
    return (1 / 2) * norm(F(x) + dF(x).dot(x_neu - x))**2
    
def F(x):
    return np.array([np.sin(x[0]) - x[1], np.exp(-x[1]) - x[0]], dtype=float)

def dF(x):
    return np.array([[np.cos(x[0]), -1], [-1, -np.exp(-x[1])]], dtype=float)
    
def gauss_newton(x0, delta):
    x = x0
    x0 = np.array(x, dtype=float)
    print("Startpunkt: (%f, %f)" % tuple(x0))
    while True:
        x_last = x
        x = x - (np.linalg.inv(dF(x))).dot(F(x))
        if norm(F(x_last) - F(x)) <= delta:
            break;
            
    x = np.array(x, dtype=float)
    return x
    
def levenberg_marquardt(x0, v, my, delta):
    i = np.identity(2)
    x = x0
    x0 = np.array(x, dtype=float)
    print("Startpunkt: (%f, %f)" % tuple(x0))
    while True:        
        p = np.linalg.solve((dF(x).T.dot(dF(x))+ my * i), (-dF(x).T.dot(F(x))))
        x_neu = x + p
        rho = (f(x) - f(x_neu)) / (q_s(x, x) - q_s(x, x_neu))        
        if rho > 0:
            x = x_neu
            my = my * max(1/3, 1 - (2 * rho - 1)**3)
            v = 2
            
        elif rho <= 0:
            my = v * my
            v = 2 * v
            
        if norm(p) <= delta:
            break;
            
    x = np.array(x, dtype=float)
    return x
    
startpunkte = np.array([[5, 2], [6, 2], [-1, -1]], dtype=float)
print("Gauß-Newton-Verfahren:")
delta = 10**-6
for x0 in startpunkte:
    x = gauss_newton(x0, delta)
    print("Ergebnis: (%f, %f)" % tuple(x) + "\n")
    
print("Levenberg-Marquardt-Verfahren:")
my = 1
v = 1
for x0 in startpunkte:
    x = levenberg_marquardt(x0, v, my, delta)
    print("Ergebnis: (%f, %f)" % tuple(x) + "\n")
