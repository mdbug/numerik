# -*- coding: utf-8 -*-
"""
Created on Wed May 23 13:41:14 2018

@author: pharnisch
"""

import math
import numpy as np
import matplotlib.pyplot as plt

# Bib (eigene Methoden aus vorherigen Aufgaben zum Lösen eines LGS)

def zerlegung(A):
    n = A[0].size
    LU = np.copy(A)
    p = [None] * (n-1)
    for i in range(n):
        # Permutationen festhalten (Indizes sind jeweils eines kleiner, damit der Umgang mit Arrays erleichtert wird)
        if(i < n - 1):
            # Prüfen, ob getauscht werden muss -> in p wird Manipulation festgehalten
            if LU[i][i] == 0:
                # save swap
                s = i
                while s < n - 1:
                    s += 1
                    if LU[s][s] != 0:
                        break
                # swap
                LU[[i, s]] = LU[[s, i]]
                p[i] = s
            else:
                p[i] = i
        # Eliminiert die i-te Spalte (Zeile für Zeile)
        for j in range(i + 1, n): # Start bei i+1, damit unter den Spalten angefangen wird, die nicht mehr verändert werden  
            # Werte der L-Matrix (in der Spalte von oben nach unten)
            l = (float(LU[j][i]) / LU[i][i])
            LU[j][i] = l
            # Verändert in der jeweiligen Zeile direkt die restlichen Werte von A mit
            for k in range(i + 1, n): # Manipulieren der A Spalten mit negativem l -> wird dann iterativ zu U im oberen Dreieck
                LU[j][k] = float(LU[j][k]) - float(l) * LU[i][k]
    return LU, p



def vorwaerts(LU, x, diagonale_ist_eins):
    n = len(x)
    y = np.copy(x) 
    # jedes y_i berechnen vom Gesamtergebnis-Vektor y
    for i in range(len(y)):
        linke_seite_ohne_y_i = 0
        rechte_seite = x[i]
        # aufkommende y_0, ..., y_i-1 -Variablen zum berechnen von y_i
        for j in range(i):
            if j == i: # y[i] ist y_i, also die Zahl, die wir berechnen wollen
                break
            linke_seite_ohne_y_i += float(LU[i][j]) * float(y[j])
        if diagonale_ist_eins:
            y_i = rechte_seite - float(linke_seite_ohne_y_i)
            y[i] = y_i
        else:
            y_i_mal_faktor = float(rechte_seite) - linke_seite_ohne_y_i
            faktor = LU[i][i]
            y_i = y_i_mal_faktor / float(faktor)
            y[i] = y_i
    return y


def rueckwaerts(LU, x, diagonale_ist_eins):
    n = len(x)
    y = np.copy(x)
    i = n-1    
    # jedes y_i berechnen vom Gesamtergebnis-Vektor y, mit der letzten Zeile anfangen
    for v in range(n):
        linke_seite_ohne_y_i = 0.
        rechte_seite = x[i] # von unten 
        # aufkommende y_n, ..., y_i+1 -Variablen zum berechnen von y_i, von rechts anfangen in der Matrix
        j = n-1
        for w in range(n-i):
            if j == i: # y[i] ist y_i, also die Zahl, die wir berechnen wollen
                break
            linke_seite_ohne_y_i += LU[i][j] * float(y[j])
            j = j - 1           
        if diagonale_ist_eins:
            y_i = rechte_seite - float(linke_seite_ohne_y_i)
            y[i] = y_i
        else:
            y_i_mal_faktor = float(rechte_seite) - linke_seite_ohne_y_i
            faktor = LU[i][i]
            y_i = y_i_mal_faktor / float(faktor)
            y[i] = y_i
        i -= 1       
    return y

def permutation(x, p, rueckwaerts):
    n = len(x)
    y = np.copy(x)
    for l in range(n - 1):
        if rueckwaerts:
            i = n - l - 2
        else:
            i = l
        swap_index = int(p[i])
        if(i == swap_index): # nichts zu tauschen
            continue
        tmp = y[swap_index]
        y[swap_index] = y[i]
        y[i] = tmp
        #y.swapaxes(i, swap_index)
    return y

def nachiteration(LU, x, A, b, epsilon):
    x_k = x
    criteria = 1000
    
    iterationen = 0
    while criteria >= epsilon:
        iterationen += 1
        print("Dies ist die " + str(iterationen) +"-te Iteration:")
        r_k = b - A.dot(x_k)
        h1 = permutation(r_k, p, 0) # LU = Pb => Permutationen werden durchgeführt zuerst
        h2 = vorwaerts(LU, h1, 1) # L * h = r_k (h = U * p_k), dabei ist L eine untere Dreiecksmatrix
        p_k = rueckwaerts(LU, h2, 0) # U * p_k = h, dabei ist U eine obere Dreiecksmatrix
        x_k = x_k + p_k
        criteria = 1.0 * np.linalg.norm(p_k, 2) / np.linalg.norm(x_k, 2)
        print("norm(p_k) / norm(x_k+1): " + str(criteria))
        print(x_k)
    return x_k

def solve(A, b):
    LU, p = zerlegung(A)
    b_permutiert = permutation(b, p, 0)
    y = vorwaerts(LU, b_permutiert, 1)
    x = rueckwaerts(LU, y, 0)
    x_nachiteriert = nachiteration(LU, x, A, b, math.pow(10, -6))
    return x_nachiteriert

# Aufgabe 2
# a)

def dividierte_differenzen(punkte):
    n = len(punkte)
    div_diff = np.empty((n, n))
    for i in range(n):
        if i == 0: # die d_ii auf die erste Spalte setzen
            for j in range(n):
                div_diff[j][0] = punkte[j][1] # bei punkt ist zweites element y
                print("d_" + str(j) + str(j-i) + ": " + str(div_diff[j][i]))
        else:
            for j in range(i, n):
                div_diff[j][i] = ( div_diff[j][i-1] - div_diff[j-1][i-1] ) / ( punkte[j][0] - punkte[j-i][0] ) # bei punkten ist erstes element x
                print("d_" + str(j) + str(j-i) + ": " + str(div_diff[j][i]))
    return np.diag(div_diff);
    
def newton_polynom_auswerten(punkte, dividierte_differenzen, x):
    n = len(punkte)
    q_k = dividierte_differenzen[n-1]
    for i in range(n-1)[::-1]:
        q_k = dividierte_differenzen[i] + (x - punkte[i][0]) * q_k
    return q_k
    
punkte = np.array([[0,3], [1, 2], [3, 6]])
dd = dividierte_differenzen(punkte)
print(dd)
x = [3, 2]
y = newton_polynom_auswerten(punkte, dd, x)
print("y: " + str(y))

# b)

m = np.array([7, 9, 11])
f = lambda x: 1 / (1 + x**2)
plt.figure(0)
for n in m:
    punkte = []
    for x in np.linspace(-5, 5, num=n, endpoint=True):
        punkte.append([x, f(x)])
    dd = dividierte_differenzen(punkte)
    print(dd)
    x = np.linspace(-5, 5, 100)
    plt.plot(x, newton_polynom_auswerten(punkte, dd, x), label = 'p_n(x) bei ' + str(n) + ' Punkten')
plt.plot(x, f(x), label = 'f(x)')
plt.legend()
plt.title('Interpolationspolynom (mit verschieden vielen Punkten)')

# c)

def get_points(n):
    points = []
    for i in range(n):
        x_i = - 5 * np.cos(math.pi * (2 * i + 1) / (2 * n))
        print('x_i: ' + str(x_i))
        print('y_i: ' + str(f(x_i)))
        points.append([x_i, f(x_i)])
    return points;
    
m = np.array([7, 9, 11])
#f = lambda x: 1 / (1 + x**2) # in b bereits definiert
plt.figure(1)
for n in m:
    punkte = get_points(n)
    dd = dividierte_differenzen(punkte)
    print(dd)
    x = np.linspace(-5, 5, 100)
    plt.plot(x, newton_polynom_auswerten(punkte, dd, x), label = 'p_n(x) bei ' + str(n) + ' Punkten')
plt.plot(x, f(x), label = 'f(x)')
plt.legend()
plt.title('Interpolationspolynom (mit speziellen Punkten)')

# d)
    
f_x = lambda t: 1 / (1 + t) * np.cos(3 * math.pi * t)
f_y = lambda t: 1 / (1 + t) * np.sin(3 * math.pi * t)
plt.figure(2)
m = np.array([6, 7, 8])
for n in m:
    punkte_x = []
    punkte_y = []
    for t in np.linspace(0, 1, n):
        punkte_x.append([t, f_x(t)])
        punkte_y.append([t, f_y(t)])
    dd_x = dividierte_differenzen(punkte_x)
    dd_y = dividierte_differenzen(punkte_y)
    print(dd_x)
    print(dd_y)    
    t = np.linspace(0, 1, 101)    
    plt.plot(newton_polynom_auswerten(punkte_x, dd_x, t), newton_polynom_auswerten(punkte_y, dd_y, t), label= 'p_n(t) bei ' + str(n) + ' Punkten')
plt.legend()
plt.title('Kurve (mit je einem Interpolationspolynom pro Dimension)')
    
# Aufgabe 4 (natürliche kubische Splines errechnen und deren Auswertung implementieren)

def natuerlichen_kubischen_spline_berechnen(punkte):
    print("punkte:")
    print(punkte)
    n = len(punkte)
    print("n:")
    print(n)
    # h erstellen (Vektor mit h_i zur Hilfe der folgenden Rechnungen)
    h = np.empty(n-1)
    for i in range(n-1):
        h[i] = punkte[i+1][0] - punkte[i][0] # auf index 0 ist das x! (Breite der n-1 Teilpolynome wird berechnet.)
    print("h:")
    print(h)
    # A erstellen
    A = np.zeros((n-2, n-2)) # es ist n-2 x n-2 
    for i in range(0, n-2):
        if i != 0:
            A[i][i-1] = h[i]
        A[i][i] = 2 * (h[i] + h[i+1])
        if i != n-3:
            A[i][i+1] = h[i+1]
    print("A")
    print(A)
    # y erstellen
    y = np.empty(n-2) # hiervon gibt es nur n-2!
    for i in range(0, n-2):
        y[i] = 6 * ( ((punkte[i+1][1] - punkte[i][1]) / h[i]) - ((punkte[i][1] - punkte[i-1][1]) / h[i-1]) )
    print("y:")
    print(y)
    # Aß = y lösen
    ß_mid = np.linalg.solve(A, y) # solve(A, y) # hat Länge n-2 (es fehlt noch eine 0 darueber und darunter fuer natuerliche Bed.)
    ß = np.zeros(n)
    ß[1:n-1] = ß_mid # ß hat Länge n!
    print("ß:")
    print(ß)
    # alpha bestimmen
    alpha = np.empty(n-1)
    for i in range(n-1):
        alpha[i] = ((punkte[i+1][1] - punkte[i][1]) / h[i]) - 1/3 * ß[i] * h[i] - 1/6 * ß[i+1] * h[i]
    print("alpha:")
    print(alpha)
    return h, ß, alpha

def natuerlichen_kubischen_spline_auswerten(punkte, h, beta, alpha, x):
    # find the right spline zone
    n = len(punkte)
    i = 0
    while x > punkte[i+1][0]:
        if i+1 == n-1:
            break
        i += 1
    return punkte[i][1] + alpha[i] * (x - punkte[i][0]) + (beta[i] / 2) * ((x - punkte[i][0])**2) + ((beta[i+1] - beta[i]) / (6 * h[i])) * ((x - punkte[i][0])**3)
    #return 1 * x

def n_k_s_auswerten(punkte, h, beta, alpha, x):
    n = len(x)
    y = np.empty(n)
    for i in range(n):
        y[i] = natuerlichen_kubischen_spline_auswerten(punkte, h, beta, alpha, x[i])
    return y

m = np.array([7, 9, 11]) # Wie viele Teilpolynome sollen verwendet werden?
f = lambda x: 1 / (1 + x**2)
plt.figure(3)
for n in m:
    punkte = []
    for x in np.linspace(-5, 5, num=n, endpoint=True):
        punkte.append([x, f(x)])
    h, beta, alpha = natuerlichen_kubischen_spline_berechnen(punkte)
    
    x = np.linspace(-5, 5, 100)
    plt.plot(x, n_k_s_auswerten(punkte, h, beta, alpha, x), label = 'p_n(x) bei ' + str(n) + ' Teilpolynomen')
plt.plot(x, f(x), label = 'f(x)')
plt.legend()
plt.title('Natürliche kubische Interpolationspolynome')












        