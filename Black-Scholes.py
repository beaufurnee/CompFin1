#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 17 11:22:32 2020

@author: Lotte Felius, Beau Furnée, Philippe Nicolau
"""
import math
import matplotlib.pyplot as plt
from scipy.integrate import quad

def integrand(z):
     return math.exp(-(1/2)*z**2)

def std_nc_pdf(d):
    
    N = (1/math.sqrt(2*math.pi))*quad(integrand, -math.inf, d)[0]
    
    return N

def black_scholes(t,st,k,T,sigma,r):
    
    d1_factor = (math.log(st/k)+(r+(sigma**2/2))*(T-t))
    
    d1 = (1/sigma * math.sqrt(T-t)) * d1_factor
    d2 = d1 - sigma*math.sqrt(T-t)
    
    
    opt_price = std_nc_pdf(d1) * st
    opt_price = opt_price - std_nc_pdf(d2) * k * math.exp(-r*(T-t))
        
    return opt_price