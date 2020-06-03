# distShape.py
import numpy as np
import math

def distShape(x,ind,mu,sigma,s):
    y = np.zeros(len(x))
    for i in ind:
        if s == 1:
            y[i] = math.exp((-(x[i]-mu)**2)/(2*sigma**2)) / (sigma*math.sqrt(2*math.pi))
        elif s == 2:
            y[i] = (1+np.tanh(sigma*(x[i]-mu)))/2
    return y
