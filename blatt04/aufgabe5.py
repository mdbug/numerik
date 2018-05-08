# -*- coding: utf-8 -*-
"""
Created on Mon May  7 13:11:10 2018

@author: tschauer
"""

from numpy import copy, array, zeros
import numpy as np

#Das Householder-Verfahren wird auf die Matrix A angewendet
def householder(A: array) -> (array, array):
    n = len(A[0])
    diagonale = np.empty(n)
    QR = copy(A)
    R = copy(A)
    for i in range(n-1):
        a = copy(R[i:,i])
        norm_a = np.linalg.norm(a)
        d = - np.sign(R[i][i]) * norm_a
        v_1 = R[i][i] - d
        v = a
        v[0] = v_1        
        Q_i = getQ_i(n, i, d, v)
        R = Q_i.dot(R)  
        QR[:, i:] = R[:, i:]
        diagonale[i] = d
        QR[i:, i] = v
        
    diagonale[n - 1] = QR[n-1][n-1]
    return QR, diagonale
        


#Q_i wird mit n als Dimenion der Matrix, dem aktuellen Iterationsschritt i,
#dem Diagonalelement d und dem Vektor v bestimmt
def getQ_i(n: int, i: int, d: int, v: array) -> array:
    v_1 = v[0]
    Q_i = np.eye(n - i) + ((1. /  (v_1 * d)) * (np.multiply(v, v.reshape(n - i, 1))))
    if i > 0:
        Q_tmp = np.eye(n)
        Q_tmp[i:, i:] = Q_i  
        Q_i = Q_tmp
    
    return Q_i
 


#Das y der Gleichung Qy = b wird über die durch das Householder-Verfahren bestimmte
#Matrix QR, Diagonale diagonale und dem gegebenen Vektor b bestimmt
def getY(QR: array, diagonale: array, b: array) -> array:
    n = len(diagonale)
    y = copy(b)
    for i in range(n-1):
        Q_i = getQ_i(n, i, diagonale[i], QR[i:,i])
        y = Q_i.dot(y)
    
    return y



#Das R des Householder-Verfahrens wird über die durch das Householder-Verfahren bestimmte
#Matrix QR und die Diagonale diagonale bestimmt
def getR(QR: array, diagonale: array) -> array:
    n = len(diagonale)
    R = triu(QR, 1)
    for i in range(n):
        R[i][i] = diagonale[i]
        
    return R



#Das LGS Rx = y wird mit der Matrix R und dem Vektor y nach x gelöst
def getX(R: array, y: array) -> array:
    return rueckwaerts(R, y)




def rueckwaerts(lu: array, x: array) -> array:
    y = copy(x)
    for i in reversed(range(x.size)):
        for j in range(i + 1, x.size):
            y[i] -= lu[i, j] * y[j]

        y[i] /= lu[i, i]

    return y

    
    

#Die Matrix A zum testen wird mit der Dimension n erstellt
def getA(n: int) -> array:
    A = np.eye(n)
    for i in range(n):
        for j in range(n):
            if i == j or j == n-1:
                A[i][j] = 1
            elif j < i:
                A[i][j] = -1
            else:
                A[i][j] = 0

    
    return A



#Der Vektor b zum testen wird mit der Dimension n erstellt
def getB(n: int) -> array:
    b = np.empty([n])
    for i in range(n):
        if i != n-1:
            b[i] = 2-i
        
        else:
            b[i] = 2-n
            
    return b

    
######################Tests###################################

#A = array([[20, 18, 44],[0, 40, 45],[-15, 24, -108]], dtype=float)
#b = array([-4, -45, 78], dtype=float)


test = array([40, 50, 60])

for i in test:
    A = getA(i)
    b = getB(i)
    QR, diagonale = householder(A)
    y = getY(QR, diagonale, b)
    R = getR(QR, diagonale)
    x = getX(R, y) 
    print(x)
    
