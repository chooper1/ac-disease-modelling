# distShape.py
import numpy as np
import math

def distShape(x,ind,mu,sigma,s):
    if s == 1:
        exp_arg = (-np.square((x-mu)))/(2*sigma**2) #check size of mu
        y = np.exp(exp_arg) / (sigma*math.sqrt(2*math.pi))
    elif s == 2:
        y = (1+np.tanh(sigma*(x-mu)))/2
    return y
