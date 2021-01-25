from scipy.special import j0, j1, jn_zeros
import numpy as np
import matplotlib.pyplot as plt

def integrar(a, b, n, g):
    h = (b-a)/n
    s = 0.0
    for k in range (0, n):
        s+= g(a+k*h)*h
    return s

alfa = jn_zeros(0, 50)
def f(x):
    return 0.00025*x

def g(x):
    return 0

def CfBesselA(f, L, n, m):
    g = lambda x:x*f(x)*j0(alfa[n]*np.sqrt(x/L))
    A = integrar(0, L, m, g)
    return A/(L*j1(alfa[n])*j1(alfa[n]))

def CFBesselB(f, L, n, m):
    g = lambda x:x*f(x)*j0(alfa[n]*np.sqrt(x/L))
    A = integrar(0, L, m, g)
    return 2*A/(L*j1(alfa[n])*j1(alfa[n])*np.sqrt(9.81*L))

def ListaCfBesselA(f, L, n, m):
    An = []
    for i in range(0, n):
        An.append(CfBesselA(f, L, i, m))
    return An

def ListaCfBesselB(f, L, n, m):
    Bn = []
    for i in range(0, n):
        Bn.append(CFBesselB(f, L, i, m))
    return Bn

def BSeries(f,g, L, n, m, k):
    cont = 0
    A = ListaCfBesselA(f, L, n, m)
    B = ListaCfBesselB(g, L, n, m)
    F = np.zeros(shape=(k, ), dtype = float)
    X = np.linspace(0, L, k)
    T = np.linspace(0, 20, 10000)
    for t in T:
        for i in range(0, len(A)):
            F = F + A[i]*j0(alfa[i]*np.sqrt(X/L))*np.cos(np.sqrt(9.8/L)*(alfa[i]/2)*t)+ B[i]*j0(alfa[i]*np.sqrt(X/L))*np.sin(np.sqrt(9.8/L)*alfa[i]*t)
        cont = cont+1
        if cont == 10:
            plt.close('all')
            plt.xlim([-5, 5])
            plt.plot(F, X)
            plt.show(block = False)
            plt.pause(0.1)
            cont = 0

BSeries(f,g, 3, 50, 131072, 1000)
