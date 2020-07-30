import math as m
import numpy as np

k = 2

T = 50 # temperatura
C = 0  # concentraciy PAV
fi = 0.264 # 26.4%
tau = 0.1 # step by time

h = 1#w1(xi1, yj, tk) - w1(xi, yj, tk)

def mW(T, C):
    return 3.03 * m.pow(10, -4) * m.exp(31.9283/T) * m.exp(0.16828 * C)

def mN(T, C):
    return 5.10746 * m.pow(10, -3) * m.exp(70.2362/T) * m.exp(-0.25714 * C)

def k0(T, C):
    return 1.47 * m.pow(10, -10) * m.pow(T, -2.1475) * m.exp(0.13563 * C)

def gammaW(T, C):
    return 0

def gammaN(T, C):
    return 7857 * 105 * m.pow(T, -2.27694) * m.exp(-0.16194 * C)

def sigmaConstW():
    return 0.1

def sigmaConstN():
    return 0.8

def kN(sigmaI, T, C):
    return m.pow((sigmaN() - sigmaI)/(sigmaN() - sigmaW()), 3)

def kW(sigmaI, T, C):
    return m.pow((sigmaI - sigmaW())/(sigmaN() - sigmaW()), 3)

def gradP(mi, k, ki, ui):
    return -mi/k*ki*ui

#sigmaW(x(i), y(j), t(k+1))
def sigmaW(xi, yj, tk, xi1, yj1):
    return tau/(fi)*((w1(xi1, yj, tk) - w1(xi, yj, tk))/h + (w2(xi, yj1, tk) - w2(xi, yj, tk))/H) + sigmaW(xi, yj, tk)

#P(x(i+1), y(j), t(k))
def Px(xi, yj, tk):
    return -h*(mN(T, C)*w1(xi, yj, tk)/(k*kN(s2, T, C)) + gammaN(T, C) * w1(xi, yj, tk)/(w1(xi, yj, tk)*m.sqrt(1 + (H*H * (P(xi1, yj, tk) - P(xi, yj, tk)) * (P(xi1, yj, tk) - P(xi, yj, tk)))/(h*h * (P(xi, yj1, tk) - P(xi, yj, tk)) * (P(xi1, yj1, tk) - P(xi, yj, tk)))))) + Px(xi, yj, tk)#P(x(i), y(j), t(k))

#P(x(i), y(j+1), t(k))
def Py(xi, yj, tk):
    return -H*(mN(T, C)*w2(xi, yj, tk)/(k*k2(s2, T, C)) + gammaN(T,C) * w2(xi, tj, tk)/(w1(xi, yj, tk)*m.sqrt(1 + (H*H * (P(xi1, yj, tk) - P(xi, yj, tk)) * (P(xi1, yj, tk) - P(xi, yj, tk)))/(h*h * (P(xi, yj1, tk) - P(xi, yj, tk)) * (P(xi1, yj1, tk) - P(xi, yj, tk)))))) + Px(xi, yj, tk)#P(x(i), y(j), t(k))

#w1(x(i), y(j), t(k+1))
def w1(xi, yj, tk, xi1, yj1):
    return k*kN(sigmaW(xi, yj, tk, xi1, yj1), T, C)/(h*mN(T,C))*(P(xi, yj, tk)-Px(xi1, yj, tk) - h*gammaN(T,C)/m.sqrt(1 + (H*H*(P(xi1, yj, tk)*P(xi1, yj, tk) - 2*P(xi1, yj, tk)*P(xi, yj, tk) + P(xi, yj, tk)*P(xi, yj, tk)))/(h*h*(P(xi, yj1, tk)*P(xi, yj1, tk) - 2*P(xi, yj1, tk)*P(xi, yj, tk) + P(xi, yj, tk)*P(xi, yj, tk)))))

#w2(x(i), y(j), t(k+1))
def w2(xi, yj, tk, yj1):
    return (Py(xi, yj1, tk)-Px(xi, yj, tk))/(H*(mN(T,C)/(k*kN(sigmaW(xi, yj, tk, xi1, yj1), T, C)) + gammaN(T,C)/(m.sqrt(1 + (H*H*(P(xi1, yj, tk)*P(xi1, yj, tk) - 2*P(xi1, yj, tk)*P(xi, yj, tk) + P(xi, yj, tk)*P(xi, yj, tk)))/(h*h*(P(xi, yj1, tk)*P(xi, yj1, tk) - 2*P(xi, yj1, tk)*P(xi, yj, tk) + P(xi, yj, tk)*P(xi, yj, tk)))))))

#u1(x(i), y(j), t(k+1))
def u1(xi, yj, tk, xi1, yj1):
    return (k*kW(sigmaW(xi, yj, tk, xi1, yj1), T, C)*(P(xi, yj, tk) - P(xi1, yj, tk)))/(mW(T, C)*H)

#u2(x(i), y(j), t(k+1))
def u2(xi, yj, tk, xi1, yj1):
    return (k*kW(sigmaW(xi, yj, tk, xi1, yj1), T, C)*(P(xi, yj, tk) - P(xi1, yj, tk)))/(mW(T, C)*H)

def decision(f, u, X = 1, T = 1, xCount = 10, tCount = 100):
    hh = X/xCount
    t = T/tCount
    coef=t/hh**2

    ix = np.zeros(xCount+1)
    iy = np.zeros(xCount+1)

    for i in range(xCount + 1):
        ix[i] = u(0, i * hh)

    for i in range(tCount):
        iy[0] = u(i * t, 0)

        for j in range(1, xCount):
            iy[j] = ix[j] + coef * (ix[j + 1] - 2 * ix[j] + ix[j - 1]) + t * f(i * t, j * hh)

        iy[xCount] = u(i * t, X)
        
        for j in range(xCount + 1):
            ix[j] = iy[j]
            
    eps = m.fabs(iy[xCount] - u(T, X))
    print("pred = {}   actual = {}  eps = {}".format(iy[xCount], u(T, X), eps))


decision(
    lambda t, x: m.sin(k * x) * (1 + k * k * (-1 + t))+1, 
    lambda t, x: (t - 1) * m.sin(x * k)+t, 
    2, 3, 10, 100)