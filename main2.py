# -*- coding: utf-8 -*-
"""
Created on Wed Aug  5 14:47:49 2020

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



# s = np.empty((I, J, K))
# w1 = np.empty((I, J, K))
# w2 = np.empty((I, J, K))
# u1 = np.empty((I, J, K))
# u2 = np.empty((I, J, K))
# P = np.empty((I, J, K))

# s.fill(1488)
# w1.fill(1488)
# w2.fill(1488)
# u1.fill(1488)
# u2.fill(1488)
# P.fill(1488)


def up_sigma(i, j, k):
    
    tf = -tau/fi
    
    s[i, j, k+1] = tf * ( s[i,j,k] - (w1[i+1, j, k]-w1[i,j,k])/h + (w2[i,j+1,k]-w2[i,j,k])/H )


def fill_P(i, j, k):
    
    
    kkm = -mu_2 / (k_param * k2(s_2, T, C))
    
    hyp = math.hypot(w1[i,j,k], w2[i,j,k])
    P[i+1, j, k] = h * (kkm * w1[i, j, k] - gamma * w1[i, j, k]/ hyp + P[i,j,k])
    P[i, j+1, k] = H * (kkm * w2[i, j, k] - gamma * w2[i, j, k]/ hyp + P[i,j,k])



def up_w(i, j, k):
    
    kk = k_param * k2(s_2, T, C)
    
    kk2 = kk/mu_2/h
    
    coef = math.sqrt(1 + (H*(P[i+1,j,k]-P[i,j,k])/h/(P[i, j+1,k]- P[i,j,k]))**2)
    #print(f'{P[i, j+1,k]- P[i,j,k]}')
    w1[i, j, k+1] = kk2 * (h * P[i, j, k]-P[i+1, j, k] - gamma*h/coef)
    w2[i, j, k+1] = w1[i, j, k+1]/gamma/H * coef * (H * P[i,j,k]-P[i,j+1,k])/ (1 + mu_2*H*w1[i,j,k+1]*coef/kk/gamma/H)


def up_u(i, j, k):
    
    kk = k_param * k1(s_1, T, C)/mu_1

    u1[i, j,k+1] = kk*(P[i,j,k]-P[i+1, j, k])/h
    u2[i, j, k+1] = kk/H * (P[i, j, k] - P[i, j+1, k])


def get_det(a11, a12, a13, a21, a22, a23, a31, a32, a33):
    
    return np.linalg.det( np.array([[a11,a12,a13],[a21,a22,a23],[a31,a32,a33]]))


def up_P(i, j, k):
    k = next_level - 1
    
    l1 = -k_param * k1(s_1, T, C)/mu_1/h
    l2 = -k_param * k1(s_1, T, C)/mu_1/H
    hh = H*H/h
    det = get_det(-1, 1, 0, -1, 0, 1, 1-hh, -1,hh)
    

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
            





P[:,:,0] = 0.5

           
up_sigma(0,0,0)
fill_P(0,0,0)
up_w(0,0,0)
up_u(0,0,0)
up_sigma(1, 1, 1)
fill_P(0,0,1)
















