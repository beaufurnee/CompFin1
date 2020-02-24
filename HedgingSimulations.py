# -*- coding: utf-8 -*-
"""
Created on Mon Feb 17 16:20:33 2020

@author: RobbyV

"""

import math
import matplotlib.pyplot as plt
from BinomialTree import black_scholes
import random
from statistics import stdev
import numpy as np

def stock_price_process(s0, sigma, time, steps, r, runs, sigma_stock, k):
    """
    To get a histogram, execute:
        outs = []
        for i in range(100):
            outs.append(stock_price_process(100, 0.2, 1, 365, 0.06, 1, False, 99))
        n, bins, patches = plt.hist(x=outs, bins='auto', range=(-2, 2))
        mu = round(sum(outs) / len(outs), 4)
        std = round(stdev(outs), 4)
        plt.axvline(mu, linestyle='dashed', linewidth=1)
    """
    
    dt = time / steps
    u = math.exp(sigma_stock*math.sqrt(dt))
    d = math.exp(-sigma_stock*math.sqrt(dt))
    p = (math.exp(r*dt) - d) / (u - d)
      

    t = 0
    s = s0

    list_t = [t]
    list_s = [s]
    
    # Portfolio
    portfolio = 0
    debt = 0
    
    delta = black_scholes(t,s,k,time,sigma,r)[0]
    
    # Selling initial option
    debt += black_scholes(t,s,k,time,sigma,r)[1]
    
    # Hedging
    portfolio += delta
    debt -= delta * s
    
    port_list = [portfolio]
    
    # Storing loans to calculate total amount due at maturity
    debt_list = [debt]
    
    # Path plotting per simulation UP UNTILL MATURITY
    for i in range(steps - 1):           
        if random.uniform(0, 1) < p:
            s = s * u
        else:
            s = s * d
            
        t += dt
            
        # Recalculating delta
        delta = black_scholes(t,s,k,time,sigma,r)[0]
        
        # Buying or selling stocks based on delta, but first loaning the necessary money
        debt_list.append(-(delta - portfolio) * s)                
        portfolio += (delta - portfolio)
        port_list.append(portfolio)
        
    # Last price change before settling        
    if random.uniform(0, 1) < p:
        s = s * u
    else:
        s = s * d
        
    t += dt
        
    end_profit = 0
    total_debt = 0
    
    # Selling stock to caller or selling portfolio to market            
    if s > k:
        # Buying remaining stocks to be able to sell to caller
        end_profit -= (1 - portfolio) * s
        end_profit += k
    else:
        end_profit += portfolio * s
    
    # Paying debt including interest
    for loan in debt_list:
        total_debt += loan * math.exp(r * ((steps - debt_list.index(loan)) * dt))
    
    end_profit += total_debt      
        
    return end_profit
    
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
    
def plot_means(list_of_frequencies, list_of_reps):
    
    for freq in list_of_frequencies:
        mus = []
        for reps in list_of_reps:
            outs = []
            for j in range(reps):
                outs.append(stock_price_process(100, 0.2, 1, freq, 0.06, 1, False, 99))
            #n, bins, patches = plt.hist(x=outs, bins='auto', range=(-2, 2))
            mu = round(sum(outs) / len(outs), 4)
            mus.append(mu)
            std = round(stdev(outs), 4)
        plt.plot(list_of_reps, mus)
        plt.text(list_of_reps[-1]*1.01, mu, str(freq), size=15)
            
    plt.show()
    
    
    
def diff_volas(stock_sigmas):
    mus = []
    stds = []
    
    for s_sig in stock_sigmas:
        outs = []
        print(s_sig)
        for i in range(2000):
            outs.append(stock_price_process(100, 0.2, 1, 365, 0.06, 1, s_sig, 99))
        mu = round(sum(outs) / len(outs), 4)
        std = round(stdev(outs), 4)
        
        mus.append(mu)
        stds.append(std)
            
    return(mus, stds)
    """
    outs = []
plt.figure(figsize=(9,6))
plt.grid()
plt.ylabel('Frequency (#)',fontsize=20)
plt.xlabel('Profit (€)',fontsize=20)
plt.tick_params(axis='both', which='major', labelsize=15)
for i in range(20000):
    outs.append(stock_price_process(100, 0.2, 1, 365, 0.06, 1, False, 99))
n, bins, patches = plt.hist(x=outs, bins=50, range=(-2, 2))
mu = round(sum(outs) / len(outs), 4)
std = round(stdev(outs), 4)
plt.axvline(mu, linestyle='dashed', linewidth=1)
plt.text(0.5, 4000, r'$\mu=$ ' + str(mu) + '\n$ \sigma=$' + str(std), size=15)
plt.title('Daily', fontsize=25)
plt.savefig('Images/II-days')
"""