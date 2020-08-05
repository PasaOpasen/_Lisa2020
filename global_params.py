# -*- coding: utf-8 -*-
"""
Created on Tue Aug  4 16:44:08 2020

@author: qtckp
"""

import math as m
import numpy as np


T = 50 # temperatura
C = 0  # concentraciy PAV
fi = 0.264 # 26.4%
tau = 0.1 # step by time

h = 1
H = 2


def mW(T, C):
    return 3.03 * m.pow(10, -4) * m.exp(31.9283/T) * m.exp(0.16828 * C)

def mN(T, C):
    return 5.10746 * m.pow(10, -3) * m.exp(70.2362/T) * m.exp(-0.25714 * C)

def k0(T, C):
    return 1.47 * m.pow(10, -10) * m.pow(T, -2.1475) * m.exp(0.13563 * C)

def gammaW(T, C):
    return 0

def gammaN(T, C):
    return 7857 * m.pow(10,5) * m.pow(T, -2.27694) * m.exp(-0.16194 * C)

def sigmaConstW():
    return 0.1

def sigmaConstN():
    return 0.8

def kN(sigmaI, T, C):
    return m.pow((sigmaConstN() - sigmaI)/(sigmaConstN() - sigmaConstW()), 3)

def kW(sigmaI, T, C):
    return m.pow((sigmaI - sigmaConstW())/(sigmaConstN() - sigmaConstW()), 3)

def gradP(mi, k, ki, ui):
    return -mi/k*ki*ui



k_param = k0(T, C)
s_1 = sigmaConstW()
s_2 = sigmaConstN()
mu_1 = mW(T, C)
mu_2 = mN(T, C)
gamma = gammaN(T, C)




def k1(sigma, T, C):
    return 1 #kW(sigma, T, C)

def k2(sigma, T, C):
    return 1 #kN(sigma, T, C)


