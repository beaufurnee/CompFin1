#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 17 11:22:32 2020

@author: Lotte Felius, Beau Furnée, Philippe Nicolau
"""

import math


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
                tree[timestep-1][node][1] = max((K - tree[timestep-1][node][0]), (math.exp((-r*dt))*((p*(tree[timestep][node][1])) + ((1-p)*(tree[timestep][min(node+1, timestep)][1])))))

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

    return tree


# Example from slides (Binomial trees, slide 22)

# sigma = 0.2
# T = 12/12
# N = 50
# r = 0.06
# s0 = 100
# K = 99
# put = 0

sigma = 0.4
T = 5/12
N = 5
r = 0.1
s0 = 50
K = 50
put = 1

tree = full_tree(K, r, s0, sigma, T, N, put)

print("")

tis = 0
for timestep in tree:
  print(tis, timestep)
  tis = tis + 1