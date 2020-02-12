# -*- coding: utf-8 -*-
"""
Created on Tue Feb  4 17:38:09 2020

Code by Philippe Nicolau, Lotte Felius & Beau Furnée
"""

# Imports
import math
import matplotlib.pyplot as plt

"""

K = strikeprice
r = risk-free rate
S0 = value stock at time 0
sigma = volatilaty
time = maturity
steps = stepsize
c = {call option = 0, put option = 1}

"""

def b_tree(K, r, s0, sigma, time, steps, c = 0):
    """
    Function to approximate the price of an option
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
            l_option[j] = (p * l_option[j] + (1-p) * l_option[j+1]) * math.exp(-]r * dt)
        l_option.pop()

    return l_option[0]

def options(init,steps_range,increment,K, r, s0, sigma, time):
    challa = []
    for steps in range(init,steps_range,increment):
        challa.append(b_tree(K, r, s0, sigma, time, steps))

    plt.figure(figsize=(8,6))
    plt.grid()
    plt.ylabel('Option Price(€)',fontsize=18)
    plt.xlabel('Step',fontsize=18)
    plt.plot(range(init,steps_range,increment),challa)
    
    return challa


#test:
options(50,500,5,99,0.06,100,0.2,1)

