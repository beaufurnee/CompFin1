# -*- coding: utf-8 -*-
"""
Created on Tue Feb  4 17:38:09 2020

Code by Philippe Nicolau, Lotte Felius & Beau Furnée
"""

# Imports
import math
import matplotlib.pyplot as plt
from scipy.integrate import quad
import numpy as np


#BLACK-SCHOLES####################
##################################

def integrand(z):
     return math.exp(-(1/2)*z**2)

def N(d):
    
    N = (1/math.sqrt(2*math.pi))*quad(integrand, -math.inf, d)[0]
    
    return N

def black_scholes(t,st,k,T,sigma,r):
#black_scholes(0,100,99,1,0.2,0.06)
    
    d1_factor = (math.log(st/k)+(r+(sigma**2/2))*(T-t))

    d1 = (1/(sigma * math.sqrt(T-t))) * d1_factor
    d2 = d1 - sigma*math.sqrt(T-t)
    
    delta = N(d1)
    
    opt_price = delta * st
    opt_price = opt_price - N(d2) * k * math.exp(-r*(T-t))
        
    return delta, opt_price




def b_tree(K, r, s0, sigma, time, steps, c = 0):
#b_tree(99,0.06,100,0.2,1,50,c = 0)

    """
    Estimate the price of an option.
    (Question I.1)

    Arguments

        K: strike price
        r: interest rate (risk-free rate)
        s0: stock price at t=0
        sigma: volatility
        time: maturity
        steps: steps in the binomial tree
        c: option type; call option if c == 0, else put option
    """

    dt = time / (steps)

    u = math.exp(sigma*math.sqrt(dt))
    d = math.exp(-sigma*math.sqrt(dt))
    p = (math.exp(r*dt) - d) / (u - d)

    # Lists to store the stock prices and option prices
    l_stocks = []
    l_option = []
    dt_options = []

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
        
        if len(l_option) == 2:
            dt_options.append(l_option[j-1])
            dt_options.append(l_option[j]) 
                 

    return dt_options, l_option[0], l_stocks


def an_opt_vs_vol(initvol,vol_max,vol_increment,t,st,k,T,r):
    
    options = []    

    initvol = initvol*100
    vol_max = vol_max*100
    vol_increment = vol_increment*100
    vol_range = range(int(initvol),int(vol_max),int(vol_increment))
    
    for sigma in vol_range:
        sigma = sigma/100
        d1_factor = (math.log(st/k)+(r+(sigma**2/2))*(T-t))
        d1 = (1/(sigma * math.sqrt(T-t))) * d1_factor
        d2 = d1 - sigma*math.sqrt(T-t)
        opt_price = N(d1) * st
        opt_price = opt_price - N(d2) * k * math.exp(-r*(T-t))
        options.append(opt_price)
    
    plt.figure(figsize=(9,6))
    plt.grid()
    plt.ylabel('Option Price (€)',fontsize=20)
    plt.xlabel('σ: Volatility (%)',fontsize=20)
    plt.tick_params(axis='both', which='major', labelsize=15)
    plt.plot(vol_range,options,linewidth = '2.5')
    plt.fill_between(vol_range,options, facecolor = 'Blue', alpha = 0.3)
    plt.savefig('Images/I-1_An_Opt_vsSigma',dpi=600)

        
    return opt_price    


def opt_vs_volatility(initvol,vol_max,vol_increment,K, r, s0, t ,T, steps, c = 0):
    
    '''    
    Evaluate option price for a range of values for sigma (the volatility)
    with plot.
    (Question I.1)
    
    Arguments
        
        initvol: Initial value for sigma 
        vol_max: Upper bound for sigma
        vol_increment: Increment taken towards vol_max
        K: strike price
        r: interest rate
        s0: stock price at t=0
        sigma: volatility
        time: maturity
        steps: number of steps taken in the binomial tree
    '''

    options = []   
    an_options = []

    initvol = initvol*100
    vol_max = vol_max*100
    vol_increment = vol_increment*100
    vol_range = range(int(initvol),int(vol_max),int(vol_increment))
    
    for sigma in vol_range:
        sigma = sigma/100
        
        #Black-Scholes
        d1_factor = (math.log(s0/K)+(r+(sigma**2/2))*(T-t))
        d1 = (1/(sigma * math.sqrt(T-t))) * d1_factor
        d2 = d1 - sigma*math.sqrt(T-t)
        opt_price = N(d1) * s0
        opt_price = opt_price - N(d2) * K * math.exp(-r*(T-t))
        an_options.append(opt_price)

        #Binomial tree
        option = b_tree(K, r, s0, sigma, T, steps, c = 0)[1]
        options.append(option)

    plt.figure(figsize=(9,6))
    plt.grid()
    plt.ylabel('Option Price (€)',fontsize=20)
    plt.xlabel('σ: Volatility (%)',fontsize=20)
    plt.tick_params(axis='both', which='major', labelsize=15)
    plt.plot(vol_range,options,color='green',linewidth = '2.5')
    plt.plot(vol_range,an_options,linewidth = '2.5')
    plt.fill_between(vol_range,options, facecolor = 'Green', alpha = 0.3)
    #plt.savefig('Images/I-1_Opt_vsSigma_An',dpi=600)

    plt.figure(figsize=(9,6))
    plt.grid()
    plt.ylabel('Difference (%)',fontsize=20)
    plt.xlabel('σ: Volatility (%)',fontsize=20)
    plt.tick_params(axis='both', which='major', labelsize=15)
    plt.plot(vol_range,[options[i]-an_options[i] for i in range(len(options))],'r',linewidth = 2.5)
    plt.savefig('Images/I-1_OP_vsSigma_Diff')
    
# =============================================================================
#     plt.figure(figsize=(9,6))
#     plt.grid()
#     plt.ylabel('Difference (%)',fontsize=20)
#     plt.xlabel('σ: Volatility (%)',fontsize=20)
#     plt.tick_params(axis='both', which='major', labelsize=15)
#     plt.plot(vol_range,[options[i]-an_options[i] for i in range(len(options))],'r',linewidth = 2.5)
#     plt.plot(vol_range,options,color='green',linewidth = '2.5')
#     plt.plot(vol_range,an_options,linewidth = '2.5')
#     plt.fill_between(vol_range,options, facecolor = 'Green', alpha = 0.3)
# =============================================================================
    
 
    
    fig, ax1 = plt.subplots(figsize=(9,6))

    color = 'tab:blue'
    ax1.set_xlabel('σ: Volatility (%)', fontsize=20)
    ax1.set_ylabel('Option Price (€)', color=color, fontsize=20)  # we already handled the x-label with ax1
    ax1.plot(vol_range,options,color='green',linewidth = '2.5')
    ax1.plot(vol_range,an_options,linewidth = '2.5')
    ax1.fill_between(vol_range,options, facecolor = 'Green', alpha = 0.3)
    ax1.tick_params(axis='y', labelcolor=color, which='major', labelsize=15)
    ax1.grid()
    
    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
    
    color = 'tab:red'
    ax2.set_ylabel('Difference between Methods(€)', color=color,fontsize=20)
    ax2.plot(vol_range,[options[i]-an_options[i] for i in range(len(options))],'r',linewidth = 2.5, alpha = 0.75)
    ax2.grid()
    ax2.tick_params(axis='y', labelcolor=color, which='major', labelsize=15)
    print(options[-1]-an_options[-1])
    fig.tight_layout()  
    plt.savefig('Images/Full_I1',dpi=600)

def options_forXsteps(init,steps_range,increment,K, r, s0, sigma, time):

    '''
    Calculate option value for certain lengths of a binomial tree with plot to
    check convergence.
    (Question I.2)

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
    plt.ylabel('Option Price (€)',fontsize=20)
    plt.xlabel('Number of Steps in the BT',fontsize=20)
    plt.tick_params(axis='both', which='major', labelsize=15)
    plt.plot(steprange,optionprices,color='green')
    plt.savefig('Images/I-2_Opt_vsSteps1750-10000',dpi=600)
    
    #return optionprices


def hedge_vs_volatility(initvol,vol_range,vol_increment,K, r, s0, time, steps):
#hedge_vs_volatility(0.01,2,0.01,99,0.06,100,1,50)
    
    
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
    
    
    hedge_parameters_analytical = []
    hedge_parameters_BT = []
      
    initvol = initvol*100
    vol_range = vol_range*100
    vol_increment = vol_increment*100
    hedge_range = range(int(initvol),int(vol_range),int(vol_increment))
    
    
    for sigma in hedge_range:
    
        dt = time / (steps)
            
        u = math.exp((sigma/100)*math.sqrt(dt))
        d = math.exp((-sigma/100)*math.sqrt(dt))
        
        #Option Prices at Maturity
        C_U = b_tree(K, r, s0, sigma/100, time, steps, c=0)[0][0]
        C_D = b_tree(K, r, s0, sigma/100, time, steps, c=0)[0][1]
        
        #to be fixed#
        deltaA = black_scholes(1/steps,s0,K,time,sigma/100,r)[0]
        deltaB = (C_U - C_D)/(s0*u-s0*d)
        #############
        
        hedge_parameters_analytical.append(deltaA)
        hedge_parameters_BT.append(deltaB)
        
    plt.figure(figsize=(8,6))   
    plt.grid()
    plt.ylabel('Δ: Number of Shares',fontsize=20)
    plt.xlabel('σ: Volatility (%)',fontsize=20)
    plt.tick_params(axis='both', which='major', labelsize=15)
    plt.plot(hedge_range,hedge_parameters_analytical,label='BS')
    plt.plot(hedge_range,hedge_parameters_BT,color='green',label='BT') 
    plt.legend()
    plt.savefig('Images/I-3_Del_vsSigma',dpi=600)


def upstock(init,maxsteps,incrsteps,s0,sigma,time):

    '''
    Graph the value of up stock at maturity given different numbers of steps.

    Arguments
    
        init: initial number of steps considered
        maxsteps: final value for the number of steps
        incrsteps: increment in the number of steps before maturity after each 
                   computation of the up stock price
        s0: stock price at t=0
        sigma: volatility
        time: maturity

    '''
    steps = range(1,5000,50)
    upstocks = []
        
    for stepstaken in steps:

        dt = time / stepstaken
        u = math.exp(sigma*math.sqrt(dt))
        
        upstocks.append(s0 * u**(stepstaken))
        
    plt.figure(figsize=(8,6))
    plt.grid()
    plt.ylabel('Number of Steps in the BT',fontsize=18)
    plt.xlabel('Price of Up-Stock at Maturity',fontsize=18)    
    plt.plot(steps,upstocks)
    
