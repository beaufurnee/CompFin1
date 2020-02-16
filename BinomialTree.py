# -*- coding: utf-8 -*-
"""
Created on Tue Feb  4 17:38:09 2020

Code by Philippe Nicolau, Lotte Felius & Beau Furnée
"""

# Imports
import math
import matplotlib.pyplot as plt


def b_tree(K, r, s0, sigma, time, steps, c = 0):

    """
    Estimate the price of an option

    Arguments

        K: strike price
        r: interest rate (risk-free rate)
        s0: stock price at t=0
        sigma: volatility
        time: maturity
        steps: steps in the binomial tree
        c: option type; call option if c == 0, else put option
    """

    dt = time / steps

    u = math.exp(sigma*math.sqrt(dt))
    d = math.exp(-sigma*math.sqrt(dt))
    p = (math.exp(r*dt) - d) / (u - d)

    # Lists to store the stock prices and option prices
    l_stocks = []
    l_option = []

    # Calculate the stock values at expiration
    for i in range(steps+1):
        l_stocks.append(s0 * u**((steps-i) - i))
    

    # Turn stock value into option value
    for el in l_stocks:

        if c == 0:
            l_option.append(max(0, el - K))
        else:
            l_option.append(max(0, K - el))


    # Move backwards through the tree calculating option prices
    while len(l_option) > 1:
        for j in range(len(l_option) - 1):
            l_option[j] = (p * l_option[j] + (1-p) * l_option[j+1])
            l_option[j] = l_option[j] * math.exp(-r * dt)
        l_option.pop()

    return l_option[0], l_stocks


def options(init,steps_range,increment,K, r, s0, sigma, time):

    '''
    Calculate option value for certain lengths of a binomial tree with plot to
    check convergence.

    Arguments

        init: initial number of steps taken in the binomial tree
        steps_range: final number of steps in the binomial tree
        increment: increment size taken from 'init' towards 'steps_range'
        K: strike price
        r: interest rate
        s0: stock price at t=0
        sigma: volatility
        time: maturity
    '''
    
    optionprices = []
    steprange = range(init,steps_range,increment)
    
    for steps in steprange:
        optionprices.append(b_tree(K, r, s0, sigma, time, steps)[0])

    plt.figure(figsize=(8,6))
    plt.grid()
    plt.ylabel('Option Price(€)',fontsize=18)
    plt.xlabel('Step',fontsize=18)
    plt.plot(steprange,optionprices)

    return optionprices


def hedge_per_volatility(initvol,vol_range,vol_increment,K, r, s0, time, steps):
    
    '''
    Loop over a range of values for sigma (the volatility) and check delta 
    (the hedging parameter)
    
    Arguments
    
        initvol: initial volatility
        vol_range: final value for the volatility
        vol_increment: increment in volatility after each computation of a delta
        K: strike price
        r: interest rate (risk-free rate)
        s0: stock price at t=0
        sigma: volatility
        time: maturity
        
    '''
    
    initvol = initvol*100
    vol_range = vol_range*100
    vol_increment = vol_increment*100
    hedge_range = range(int(initvol),int(vol_range),int(vol_increment))
    hedge_parameters = []
    
    for sigma in hedge_range:
    
        
        #Stock Prices at Maturity
        S_TU = b_tree(K, r, s0, sigma/100, time, steps)[1][0]
        S_TD = b_tree(K, r, s0, sigma/100, time, steps)[1][steps]
        #Option Prices at Maturity
        C_U = S_TU-K
        C_D = 0
        
        delta = (C_U - C_D)/(S_TU-S_TD)
        hedge_parameters.append(delta)
        
    plt.figure(figsize=(8,6))
    plt.grid()
    plt.ylabel('Delta',fontsize=18)
    plt.xlabel('Volatility',fontsize=18)
    plt.plot(hedge_range,hedge_parameters)


#options(1,10,1,99,0.06,100,0.2,1)

steps = range(1,2000,2)
topdogs = []
    
for stepstaken in steps:
    
    s0 = 100
    sigma = 0.2
    time = 1
    dt = time / stepstaken
    u = math.exp(sigma*math.sqrt(dt))
    
    topdogs.append(s0 * u**(stepstaken))
    
    
