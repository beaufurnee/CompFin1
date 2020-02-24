#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 17 11:22:32 2020

@author: Lotte Felius, Beau Furnée, Philippe Nicolau
"""

import math
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import quad
from BinomialTree import black_scholes
from BinomialTree import b_tree
"""

Create tree with expected stock values

N = steps
S0 = Stock value at t = 0
sigma = volatility
T = time until maturity
K = strike price
r = risk free rate
c = {call = 0, put = 1}

"""


def full_tree(K, r, S0, sigma, time, steps, c = 0):
    """
    Function to approximate price of an American option
    """

    dt = time / steps
    tree = []


    ## Create tree with coordinates (node(i,j))

    for i in range(N + 1):
        tree = [[[i, j + 1] for j in range(i + 1)] for i in range(N + 1)]

    tree[0][0][0] = S0

    u = math.exp(sigma * math.sqrt(dt))
    d = math.exp(-sigma * (math.sqrt(dt)))
    p = ((math.exp(r*dt))-d)/(u-d)

    for timestep in range(N + 1):

        if timestep > 1:
            for node in range(timestep + 1):
                if tree[timestep][node][1] > 1:
                    tree[timestep][node][0] = tree[timestep - 1][node - 1][0] * d
                elif tree[timestep][node][1] == 1:
                    tree[timestep][node][0] = tree[timestep - 1][node][0] * u

        elif timestep == 1:
            for i in range(2):
                if i == 0:
                    tree[timestep][i][0] = tree[timestep - 1][0][0] * u
                elif i == 1:
                    tree[timestep][i][0] = tree[timestep - 1][0][0] * d


    # calculate option values at Maturity

    if c == 1:
        for node in range(N+1):
            tree[N][node][1] = max(0, (K - tree[N][node][0]))
    else:
        for node in range(N+1):
            tree[N][node][1] = max(0, (tree[N][node][0]-K))

    # Calculate all option values and store them

    if c == 1:
        for timestep in range(N, -1, -1):
            for node in range(timestep):
                tree[timestep-1][node][1] = max((K - tree[timestep-1][node][0]), (math.exp((-r*dt))*(
                    (p*(tree[timestep][node][1])) + ((1-p)*(tree[timestep][min(node+1, timestep)][1])))))

    else:
        for timestep in range(N, -1, -1):
            for node in range(timestep):
                tree[timestep - 1][node][1] = max((tree[timestep-1][node][0] - K), (math.exp((-r * dt)) * (
                (p * (tree[timestep][node][1])) + ((1 - p) * (tree[timestep][min(node + 1, timestep)][1])))))

    # Just for representation

    def cp(c):
        if c == 1:
            return "put"
        if c == 0:
            return "call"
    cp = str(cp(c))

    print("Value %s option at time 0 is " % cp + str(tree[0][0][1]))

    # return tree
    return tree[0][0][1]


sigma = 0.2
T = 12/12
N = 50
r = 0.06
s0 = 100
K = 99
put = 1

full_tree(K, r, s0, sigma, T, N, 0)
full_tree(K, r, s0, sigma, T, N, 1)

# tree = full_tree(K, r, s0, sigma, T, N, put)
#
# print("")
#
# tis = 0
# for timestep in tree:
#   print(tis, timestep)
#   tis = tis + 1


def opt_vs_volatility(initvol, vol_max, vol_increment, K, r, s0, sigma, time, steps):
    '''
    Evaluate option price for a range of values for sigma (the volatility) with plot.
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

    option_c = []
    option_p = []

    initvol = initvol * 100
    vol_max = vol_max * 100
    vol_increment = vol_increment * 100
    vol_range = range(int(initvol), int(vol_max), int(vol_increment))


    for sigma in vol_range:
        optionc = full_tree(K, r, s0, sigma / 100, time, steps, c=0)
        optionp = full_tree(K, r, s0, sigma / 100, time, steps, c=1)
        option_c.append(optionc)
        option_p.append(optionp)

    plt.figure(figsize=(8, 6))
    plt.grid()
    plt.ylabel('Option Price (€)', fontsize=20)
    plt.xlabel('σ: Volatility (%)', fontsize=20)
    plt.tick_params(axis='both', which='major', labelsize=15)
    plt.plot(vol_range, option_c, color='yellowgreen')
    plt.plot(vol_range, option_p, color='indianred')
    plt.legend(['Call Option', 'Put Option'], loc='upper left')
    plt.show()
    plt.savefig('Images/AI-1_Opt_vsSigma', dpi=600)

# Plot
opt_vs_volatility(0.01, 6.1, 0.1, K, r, s0, sigma, T, N)

# print(b_tree(K, r, s0, sigma, T, N, c = 1))




