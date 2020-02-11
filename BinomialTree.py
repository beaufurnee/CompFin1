# -*- coding: utf-8 -*-
"""
Created on Tue Feb  4 17:38:09 2020

Code by Philippe Nicolau, Lotte Felius & Beau Furnée
"""

# Imports
import math

def b_tree(K, r, s0, sigma, time, steps):
    """
    Function to approximate the price of an option
    """
    dt = time / steps
    
    u = math.exp(sigma*math.sqrt(dt))
    d = math.exp(-sigma*math.sqrt(dt))
    
    p = (math.exp(r*dt) - d) / (u - d)
    print(math.exp(r*dt))
    
    # Lists to store the stock prices and option prices
    l_stocks = []
    l_option = []
    
    # Calculate the stock values at expiration 
    for i in range(steps):
        l_stocks.append(s0 * u**((steps-1-i) - i))
        
  
    # Move backwards through the tree calculating option prices
    
    
    return l_stocks