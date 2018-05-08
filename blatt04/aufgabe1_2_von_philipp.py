# -*- coding: utf-8 -*-
"""
Created on Tue May  8 20:26:17 2018

@author: tschauer
"""

from pylab import *

import numpy as np


#####################################################################################################################
################################################ Aufgabe 1 ##########################################################
#####################################################################################################################
print("########################################## Aufgabe 1")

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
        swap_index = p[i]
        if(i == swap_index): # nichts zu tauschen
            continue
        tmp = y[swap_index]
        y[swap_index] = y[i]
        y[i] = tmp
        #y.swapaxes(i, swap_index)
    return y


# 1. LU-Zerlegung von A erhalten
A = array([[0,0,0,1],[2,1,2,0],[4,4,0,0],[2,3,1,0]], dtype=float)
c = array([152,154,56,17], dtype=float)
LU, p = zerlegung(A)
# 2. U.T * v = c (U.T ist untere Dreiecksmatrix -> vorw einsetzen)
v = vorwaerts(LU.T, c, 0)
#v = np.linalg.solve(triu(LU).T, c)
# 3. L.T * w = v (L.T ist obere Dreiecksmatrix -> rückw einsetzen)
w = rueckwaerts(LU.T, v, 1)
#w =  np.linalg.solve(LU.T, v)
# 4. P * z = w bzw. z = P.T * w (w in umgekehrter Reihenfolge von p vertauschen)
z = permutation(w, p, 1)

print("Die Lösung von A.T z = c ist: z = " + str(z))


# Gegentest

A = array([[0,0,0,1],[2,1,2,0],[4,4,0,0],[2,3,1,0]], dtype=float)
c = array([152,154,56,17], dtype=float)
print("np.linalg.solve(A.T, c): " + str(np.linalg.solve(A.T, c)))


#####################################################################################################################
################################################ Aufgabe 2 ##########################################################
#####################################################################################################################
print("########################################## Aufgabe 2")

def hager(LU, p):
    n = len(p) + 1
    # 0. (Initialisierung) x = (1, ..., 1).T / n
    x = ones(n) / n
    
    for i in range(5):

        # 1. y = B.T * x (für B = A^-1: U.T * v = x, L.T * w = v, y = P.T * w, also wie in Aufgabe 1)
        tmpx = np.copy(x)
        v = vorwaerts(LU.T, tmpx, 0)
        w = rueckwaerts(LU.T, v, 1)
        y = permutation(w, p, 1)
        
        # 2. e = sign(y) (komponentenweise)
        xi = np.sign(y)
        
        # 3. z = B * e (für B = A^-1: L * v = P * e, U * z = v)
        v = vorwaerts(LU, permutation(xi, p, 0), 1)
        z = rueckwaerts(LU, v, 0)
        
        # 4. falls ||z||unendlich <= z.T * x:
        if max(np.abs(z)) <= z.T.dot(x):
            # ||B||unendlich ist ungefähr ||y||eins
            return sum(np.abs(y))
        
        # 5. bestimme j für welches |z_j| maximal ist
        j = np.argmax(np.abs(z))
        
        # 6. x = j-ter Einheitsvektor
        x = np.eye(n)[j]
    
def zeilenSummenNorm(A):
    n = len(A[0])
    ret = 0
    for i in range(n): # jede Zeile
        betragsmaessige_zeilensumme = 0
        for j in range(n): # jedes Element der Zeile
            betragsmaessige_zeilensumme += np.abs(A[i][j])
        if betragsmaessige_zeilensumme > ret:
            ret = betragsmaessige_zeilensumme
    return ret
    

testgroessen = array([math.pow(10, -8), math.pow(10, -10), math.pow(10, -12)])
for i in range(len(testgroessen)):
    testgroesse = testgroessen[i]
    print("Testgröße: " + str(testgroesse))
    A = array([[3,2,1],[2, 2*testgroesse, 2*testgroesse],[1, 2*testgroesse, -testgroesse]], dtype=float)
    # eigene Implementierung
    LU, p = zerlegung(A)
    hag = hager(LU, p)
    norma = zeilenSummenNorm(A)
    print("hager - Zeilensummennorm von A Inverse: " + str(hag))
    print("Zeilensummennorm von A: " + str(norma))
    print("=> Kondition von A: " + str(norma * hag)) # np.linalg.norm(A, np.inf)
    # numpy.linalg
    print("numpy.linalg.cond(A): " + str(np.linalg.cond(A, np.inf)))
    print("")
