# -*- coding: utf-8 -*-
"""
Created on Tue Aug  4 16:48:09 2020

@author: qtckp
"""

from global_params import *
import numpy as np
import math


I = 10

J = 10

K = 20




s = np.zeros((I, J, K))
w1 = np.zeros((I, J, K))
w2 = np.zeros((I, J, K))
u1 = np.zeros((I, J, K))
u2 = np.zeros((I, J, K))
P = np.zeros((I, J, K))



def up_sigma(next_level):
    k = next_level - 1
    
    tf = -tau/fi
    
    for i in range(I-1):
        for j in range(J-1):
            s[i, j, k+1] = tf * (s[i,j,k] - (w1[i+1, j, k]-w1[i,j,k])/h + (w2[i,j+1,k]-w2[i,j,k])/H)


def fill_P(which_level):
    
    k = which_level
    
    kkm = -mu_2 / (k_param * k2(s_2, T, C))
    
    for i in range(I):
        for j in range(J):
            hyp = math.hypot(w1[i,j,k], w2[i,j,k])
            P[i+1, j, k] = h * (kkm * w1[i, j, k] - gamma * w1[i, j, k]/ hyp + P[i,j,k])
            P[i, j+1, k] = H * (kkm * w2[i, j, k] - gamma * w2[i, j, k]/ hyp + P[i,j,k])



def up_w(next_level):
    k = next_level - 1
    
    kk = k_param * k2(s_2, T, C)
    
    for i in range(I-1):
        for j in range(J-1):
            coef = math.sqrt(1 + (H*(P[i+1,j,k]-P[i,j,k])/h/(P[i, j+1,k]- P[i,j,k]))**2)
            w1[i, j, k+1] = kk/mu_2/h * (h * P[i, j, k]-P[i+1, j, k] - gamma*h/coef)
            w2[i, j, k+1] = w1[i, j, k+1]/gamma/H * coef * (H * P[i,j,k]-P[i,j+1,k])/ (1+ mu_2*H*w1[i,j,k+1]*coef/kk/gamma/H)


def up_u(next_level):
    k = next_level - 1
    
    kk = k_param * k1(s_1, T, C)/mu_1
    
    for i in range(I-1):
        for j in range(J-1):
            u1[i, j,k+1] = kk*(P[i,j,k]-P[i+1, j, k])/h
            u2[i, j, k+1] = kk/H * (P[i, j, k] - P[i, j+1, k])


def get_det(a11, a12, a13, a21, a22, a23, a31, a32, a33):
    
    return np.linalg.det( np.array([[a11,a12,a13],[a21,a22,a23],[a31,a32,a33]]))


def up_P(next_level):
    k = next_level - 1
    
    l1 = -k_param * k1(s_1, T, C)/mu_1/h
    l2 = -k_param * k1(s_1, T, C)/mu_1/H
    hh = H*H/h
    det = get_det(-1, 1, 0, -1, 0, 1, 1-hh, -1,hh)
    
    for i in range(I-1):
        for j in range(J-1):
            l3 = h * (s[i, j+1, k+1] + fi * s[i, j, k+2]/tau)
            l4 = h/H * u2[i, j+1, k+1]
            l5 = u1[i, j, k+1]/l1
            l6 = u2[i, j, k+1]/l2
            l7 = (l3-u1[i+1, j, k+1]-l4*u2[i, j+1, k+1])/l1
            
            d1 = get_det(l5,1,0,l6,0,1,l7,-1,hh)
            d2 = get_det(-1, l5, 0, -1, l6, 1, 1-hh, l7, hh)
            d3 = get_det(-1, 1, l5, -1, 0, l6, 1-hh, -1, l7)
            
            P[i,j,k+1] = d1/det
            P[i+1,j,k+1] = d2/det
            P[i,j+1,k+1] = d3/det
            
            

P[:,:,0] = np.random.normal(3,2,(I,J))
            

for k in range(1, K-1):
    print(k)
    up_sigma(k)
    up_w(k)
    up_sigma(k+1)
    up_u(k)
    up_P(k)

















