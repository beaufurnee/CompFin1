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

def stock_price_process(s0, sigma, time, steps, r, runs, plot):
    dt = time / steps
    u = math.exp(sigma*math.sqrt(dt))
    d = math.exp(-sigma*math.sqrt(dt))
    p = (math.exp(r*dt) - d) / (u - d)
    
    # if runs > 1, running multiple simulations to plot multiple paths
    for i in range(runs + 1):
        t = 0
        s = s0
    
        list_t = [t]
        list_s = [s]
        
        # Path plotting per simulation
        for i in range(steps):
            t += dt
            
            if random.uniform(0, 1) < p:
                s = s * u
            else:
                s = s * d
            
            if plot == True:
                list_t.append(t)
                list_s.append(s)
        
        if plot == True:
            x, y = list_t, list_s
            plt.plot(x, y)

    if plot == True:
        plt.show()
    else:
        return s
    
def plot_histogram(samples):
    outs = []
    for i in range(samples):
        outs.append(stock_price_process(100, 0.2, 1, 50, 0.06, 1, False))
        
    mu = round(sum(outs) / len(outs), 4)
    std = round(stdev(outs), 4)
    
    n, bins, patches = plt.hist(x=outs, bins=np.linspace(0, 202, 101), range=(0, 200), color='#0504aa')
    plt.axvline(mu, color='k', linestyle='dashed', linewidth=1)
    plt.grid(axis='y', alpha=0.75)
    plt.xlabel('Stock Price', size=15)
    plt.ylabel('Frequency', size=15)
    plt.title('Distribution of Stock Price Process runs', size=15)
    plt.text(0, 81000, r'$\mu=$ ' + str(mu) + '\n$ \sigma=$' + str(std), size=15)