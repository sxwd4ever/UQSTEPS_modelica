# -*- coding: utf-8 -*-
"""
Created on Wed Aug 19 11:43:19 2020

@author: uqxsui
"""

import math

d_c = 2e-3
A_c = math.pi * d_c ** 2 / 8

mdot = 1

G = mdot / A_c

peri_c = d_c * math.pi / 2 + d_c

d_h = 4 * A_c / peri_c

mu = 3.36364e-005

Re = G * d_h / mu