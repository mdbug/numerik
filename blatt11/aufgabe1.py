# -*- coding: utf-8 -*-
"""
Created on Tue Jun 26 09:26:26 2018

@author: tschauer
"""

#from sympy import *
import numpy as np
from numpy import array, inf, linspace, meshgrid, log, copy
from numpy.linalg import norm
from scipy.linalg import eigh
import matplotlib.pyplot as plt


def f(x):
    return (np.sin(x[0]) - x[1]) ** 2 + (np.exp(-x[1]) - x[0]) ** 2

def df(x):
    df = [2*np.cos(x[0]) * (np.sin(x[0]) - x[1]) - 2*(np.exp(-x[1]) - x[0]), (-2*np.exp(-x[1])*(np.exp(-x[1]) - x[0]) - 2 * (np.sin(x[0]) - x[1]))]
    return array(df, dtype=float)

def d2f(x):
    d2f11 = -2 * np.sin(x[0]) * (np.sin(x[0])-x[1])+2*np.cos(x[0]) ** 2 + 2
    d2f22 = 2 * np.exp(-2 * x[1]) * (np.exp(2 * x[1]) - x[0] *np.exp(x[1]) + 2)
    d2f12 = 2*np.exp(-x[1]) - 2 * np.cos(x[0])
    return array([[d2f11, d2f12], [d2f12, d2f22]], dtype=float)


def levenberg_marquardt(x0, my0):
    x = x0
    my = my0
    i = np.identity(2) 
    while not is_spd(d2f(x) + my * i):
        my = my * 2
        
    p = np.linalg.solve(d2f(x) + my * i, -df(x))
    x = x + p
    print(x)
    return 0


def is_spd(x):
    eigenwerte, eigenvektoren = eigh(x)
    for i in range(len(eigenwerte)):
        if eigenwerte[i] <= 0:
            return False
    
    return True        
    

for startpukt in [(5,2), (6,2), (-1,-1), (-2,-2)]:
    x0 = array(startpukt, dtype=float)
    print("Startpunkt: (%f, %f)" % tuple(x0))
    print(d2f(x0))
    print("")
    levenberg_marquardt(x0, 1)