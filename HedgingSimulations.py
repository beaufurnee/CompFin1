# -*- coding: utf-8 -*-
"""
Created on Mon Feb 17 16:20:33 2020

@author: RobbyV

QUESTIONS:
    
    - Are we replicating the process correctly?
    - Does: 'Frequency of hedge adjustments' == delta t?
    - 
"""

import math
import matplotlib.pyplot as plt
#import BinomialTree
import random
from statistics import stdev
import numpy as np

def stock_price_process(s0, sigma, time, steps, r):
    dt = time / steps
    u = math.exp(sigma*math.sqrt(dt))
    d = math.exp(-sigma*math.sqrt(dt))
    p = (math.exp(r*dt) - d) / (u - d)
    t = 0
    
    list_t = [t]
    list_s = [s0]
    
    for i in range(steps):
        t += dt
        
        if random.uniform(0, 1) < p:
            s0 = s0 * u
        else:
            s0 = s0 * d
            
        list_t.append(t)
        list_s.append(s0)
    
    return s0
    
def plot_histogram(samples):
    outs = []
    for i in range(samples):
        outs.append(stock_price_process(100, 0.2, 1, 50, 0.06))
        
    mu = round(sum(outs) / len(outs), 4)
    std = round(stdev(outs), 4)
    
    n, bins, patches = plt.hist(x=outs, bins=np.linspace(0, 202, 101), range=(0, 200), color='#0504aa')
    plt.axvline(mu, color='k', linestyle='dashed', linewidth=1)
    plt.grid(axis='y', alpha=0.75)
    plt.xlabel('Stock Price', size=15)
    plt.ylabel('Frequency', size=15)
    plt.title('Distribution of Stock Price Process runs', size=15)
    plt.text(0, 81000, r'$\mu=$ ' + str(mu) + '\n$ \sigma=$' + str(std), size=15)
    
#x, y = stock_price_process(100, 0.2, 1, 50, 0.06)
#plt.plot(x, y)
#plt.show()