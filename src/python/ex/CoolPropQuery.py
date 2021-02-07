# -*- coding: utf-8 -*-
"""
Created on Sun Aug 16 17:48:28 2020

@author: uqxsui
"""

from __future__ import print_function
import CoolProp as CP
from CoolProp import AbstractState
from CoolProp.CoolProp import PhaseSI, PropsSI, get_global_param_string
import CoolProp.CoolProp as CoolProp
from CoolProp.HumidAirProp import HAPropsSI
import CoolProp.State 
from math import log10, sin
from numpy.core.fromnumeric import mean

from numpy.lib.function_base import average

# P, h -> T ([Pa], [J/kg] -> [degC])
# result = CP.CoolProp.PropsSI("T", "P", 20e6, 'H', 995164, 'CO2') - 273.15

# P, T -> h ([Pa], [K] -> [J/kg])
# result = CP.CoolProp.PropsSI("H", "P", 20e6, 'T', 715 + 273.15, 'CO2')

# P, T -> mu (local viscosity)([Pa], [K] -> [mu])
# result = CP.CoolProp.PropsSI("V", "P", 3616780, 'T', 700 + 273.15, 'CO2')

# result = CP.CoolProp.PropsSI("H", "P", 8e6, 'T', 15 + 273.15, 'CO2')

# print("result=" + str(result))

# result = CP.CoolProp.PropsSI("H", "P", 20e6, 'T', 15 + 273.15, 'CO2')

# print("result=" + str(result))

# result = CP.CoolProp.PropsSI("H", "P", 8e6, 'T', 700 + 273.15, 'CO2')

# print("result=" + str(result))
# 2e7,  1.01549e6
# result = CP.CoolProp.PropsSI("T", "P", 0.8e7, 'H', 1.07878e6, 'CO2') - 273.15
# print("result=" + str(result))

# # 2e7,  1.01549e6
# # 2.01077e7 518847
# result = CP.CoolProp.PropsSI("T", "P", 2.01077e7, 'H', 518847, 'CO2') - 273.15

# print("result=" + str(result))

# 2e7,  1.01549e6
# result = CP.CoolProp.PropsSI("S", "P", 2e7, 'H', 1.26204e6, 'CO2')
# print("result=" + str(result))

# result = CP.CoolProp.PropsSI("H", "P", 9e6, 'S', result, 'CO2')
# print("result=" + str(result))

# result = CP.CoolProp.PropsSI("T", "P", 9e6, 'H', 1.23596e6, 'CO2') - 273.15
# print("result=" + str(result))

# print('high pressure, low T')
# # density
# p, T = 20e6, 273.15 + 162.144

# result = PropsSI("DMASS", "P", p, "T" , T, 'CO2')
# print("result=" + str(result))

# result = PropsSI("DMASS", "P", p, "T" , T, 'Helium')
# print("result=" + str(result))

# # conductivity
# result = PropsSI("CONDUCTIVITY", "P", p, "T" , T, 'CO2')
# print("result=" + str(result))

# result = PropsSI("CONDUCTIVITY", "P", p, "T" , T, 'Helium')
# print("result=" + str(result))

# # 
# result = PropsSI("VISCOSITY", "P", p, "T" , T, 'CO2')
# print("result=" + str(result))

# result = PropsSI("VISCOSITY", "P", p, "T" , T, 'Helium')
# print("result=" + str(result))

# print('\nlow pressure, high T')
# # density
# p, T = 8e6, 273.15 + 600

# result = PropsSI("DMASS", "P", p, "T" , T, 'CO2')
# print("result=" + str(result))

# result = PropsSI("DMASS", "P", p, "T" , T, 'Helium')
# print("result=" + str(result))

# # conductivity
# result = PropsSI("CONDUCTIVITY", "P", p, "T" , T, 'CO2')
# print("result=" + str(result))

# result = PropsSI("CONDUCTIVITY", "P", p, "T" , T, 'Helium')
# print("result=" + str(result))

# # 
# result = PropsSI("VISCOSITY", "P", p, "T" , T, 'CO2')
# print("result=" + str(result))

# result = PropsSI("VISCOSITY", "P", p, "T" , T, 'Helium')
# print("result=" + str(result))


# print('high pressure, High T')
# # density
# p, T = 20e6, 273.15 + 710

# result = PropsSI("DMASS", "P", p, "T" , T, 'CO2')
# print("result=" + str(result))

# result = PropsSI("DMASS", "P", p, "T" , T, 'Helium')
# print("result=" + str(result))

# # conductivity
# result = PropsSI("CONDUCTIVITY", "P", p, "T" , T, 'CO2')
# print("result=" + str(result))

# result = PropsSI("CONDUCTIVITY", "P", p, "T" , T, 'Helium')
# print("result=" + str(result))

# # 
# result = PropsSI("VISCOSITY", "P", p, "T" , T, 'CO2')
# print("result=" + str(result))

# result = PropsSI("VISCOSITY", "P", p, "T" , T, 'Helium')
# print("result=" + str(result))

# print('\nlow pressure, high T')
# # density
# p, T = 8e6, 273.15 + 600

# result = PropsSI("DMASS", "P", p, "T" , T, 'CO2')
# print("result=" + str(result))

# result = PropsSI("DMASS", "P", p, "T" , T, 'Helium')
# print("result=" + str(result))

# # conductivity
# result = PropsSI("CONDUCTIVITY", "P", p, "T" , T, 'CO2')
# print("result=" + str(result))

# result = PropsSI("CONDUCTIVITY", "P", p, "T" , T, 'Helium')
# print("result=" + str(result))

# # 
# result = PropsSI("VISCOSITY", "P", p, "T" , T, 'CO2')
# print("result=" + str(result))

# result = PropsSI("VISCOSITY", "P", p, "T" , T, 'Helium')
# print("result=" + str(result))

# Ra = 50e-6
# D = 2e-3
# Re = 22000
# A = -2.0 * log10(Ra / D / 3.7 + 12 / Re  )
# B = -2.0 * log10(Ra / D / 3.7 + 2.51 * A / Re)
# C = 4.781
# f = 0.25 * ((C - (A - C)** 2 / (B - 2 * A + C))**(-2))
# print(f"result: A={A}, B={B}, f={f}")
# calculate the eta of Compressor of 'D:\sxwd\Projects\UQ\2019.10.01 Thermal Cycle Simulation\doc\UQMECH05_99_CL02_B 1 MW Power Block (RCBC)_2.pdf' 
# h_in = 309.42e3
# h_out = 329.94e3

# kc_cf_z_set = [10, 12, 15]    
# kc_cf_s_set = [1.0, 1.2, 1.5]

# import itertools

# for x in itertools.product(kc_cf_s_set, kc_cf_z_set):
#         print(x)

# print('high pressure, low T')
# # density
# p, T = 20e6, 273.15 + 162.144

# result = PropsSI("DMASS", "P", p, "T" , T, 'CO2')
# print("result=" + str(result))

# # conductivity
# result = PropsSI("CONDUCTIVITY", "P", p, "T" , T, 'CO2')
# print("result=" + str(result))

# # 
# result = PropsSI("VISCOSITY", "P", p, "T" , T, 'CO2')
# print("result=" + str(result))

# print('high pressure, High T')
# # density
# p, T = 20e6, 273.15 + 710

# result = PropsSI("DMASS", "P", p, "T" , T, 'CO2')
# print("result=" + str(result))

# # conductivity
# result = PropsSI("CONDUCTIVITY", "P", p, "T" , T, 'CO2')
# print("result=" + str(result))

# # 
# result = PropsSI("VISCOSITY", "P", p, "T" , T, 'CO2')
# print("result=" + str(result))

# print('\nlow pressure, high T')
# # density
# p, T = 8e6, 273.15 + 600

# result = PropsSI("DMASS", "P", p, "T" , T, 'CO2')
# print("result=" + str(result))

# # conductivity
# result = PropsSI("CONDUCTIVITY", "P", p, "T" , T, 'CO2')
# print("result=" + str(result))

# # 
# result = PropsSI("VISCOSITY", "P", p, "T" , T, 'CO2')
# print("result=" + str(result))

# print('\nlow pressure, low T')
# # density
# p, T = 8e6, 273.15 + 162.144

# result = PropsSI("DMASS", "P", p, "T" , T, 'CO2')
# print("result=" + str(result))

# # conductivity
# result = PropsSI("CONDUCTIVITY", "P", p, "T" , T, 'CO2')
# print("result=" + str(result))

# # 
# result = PropsSI("VISCOSITY", "P", p, "T" , T, 'CO2')
# print("result=" + str(result))


# s = CP.CoolProp.PropsSI("S", "H", h_in, "P", 7.91e6, 'CO2')

# h_isen = CP.CoolProp.PropsSI("H", "S", s, "P", 20.04e6, 'CO2')

# eta =  (h_isen - h_in) / (h_out - h_in)

# print("result=" + str(eta))

# # recompressor
# h_in = 471.73e3
# h_out = 525.39e3

# s = CP.CoolProp.PropsSI("S", "H", h_in, "P", 7.91e6, 'CO2')

# h_isen = CP.CoolProp.PropsSI("H", "S", s, "P", 20.04e6, 'CO2')

# eta =  (h_isen - h_in) / (h_out - h_in)

# print("result=" + str(eta))

# # turbine

# # recompressor
# h_in = 1223.3e3
# h_out = 1078.89e3

# s = CP.CoolProp.PropsSI("S", "H", h_in, "P", 20e6, 'CO2')

# h_isen = CP.CoolProp.PropsSI("H", "S", s, "P", 8e6, 'CO2')

# eta =  (h_in - h_out) / (h_in - h_isen)

# print("result=" + str(eta))

# AS_SAT = AbstractState("HEOS", "CO2")
# AS_SAT.update(CoolProp.PT_INPUTS, 20e6, 273.15+700)
# print("First saturation derivative:", AS_SAT.hmass(), "Pa/K")

def Nu_cor(text, T_b, T_w, p_b, p_w, G, d, verbose=False):
    rho_w = PropsSI("D", "T", T_w, "P", p_w, media)
    rho_b = PropsSI("D", "T", T_b, "P", p_b, media)
    mu_b = PropsSI("VISCOSITY", "T", T_b, "P", p_b, media)
    cp_b = PropsSI("CPMASS", "T", T_b, "P", p_b, media)
    k_b = PropsSI("CONDUCTIVITY", "T", T_b, "P", p_b, media)
    h_w = PropsSI("HMASS", "T", T_w, "P", p_w, media)
    h_b = PropsSI("HMASS", "T", T_b, "P", p_b, media)
    cp_bar = (h_w - h_b) / (T_w - T_b)

    Gr = abs((rho_w - rho_b) * rho_b * g * (d**3) / (mu_b ** 2))
    Re = G * d / mu_b
    Pr = cp_b * mu_b / k_b

    q1 = Gr / (Re ** 2)
    q2 = rho_w / rho_b
    q3 = cp_bar / cp_b

    Nu = 0.124 * (Re ** 0.8) * (Pr ** 0.4) * (q1 ** 0.203) * (q2 ** 0.842) * (q3 ** 0.284)
    # if q2 > 1:
    #     q2 = 1 / q2
    if verbose:
        print(f'{text}: Nu={Nu}, q1={q1}, q2={q2}, q3={q3}')

    return (Nu, q1, q2, q3)

# bcs = [bc1, bc2, ...]
# bcn = (st_hot_in, st_cold_in)
# st = (T [K], p [bar->ba])

def predict_C1(T_h, p_h, T_c, p_c, G, d, medium, verbose=False) -> float:
    '''
        use boundary conditions to predict C1
    '''
    sides ={
        "hot": {
            "text": "for hot side",
            "T_b": T_h,
            "T_w": (T_h + T_c)/2,
            "p_b": p_h, 
            "p_w": p_h,
            "G": G,
            "d": d,
            "medium": medium,
            "verbose": verbose,
        },
        "cold": {
            "text": "for cold side",
            "T_b": T_c,
            "T_w": (T_h + T_c)/2,
            "p_b": p_c, 
            "p_w": p_c,
            "G": G,
            "d": d,
            "medium": medium,
            "verbose": verbose   
        }
    }

    (Nu, q1, q2, q3) = ([], [], [], [])

    for side in sides.values():
        (v1, v2, v3, v4) = self.Nu_cor(**side)
        Nu.append(v1)
        q1.append(v2)
        q2.append(v3)
        q3.append(v4)   

    if verbose:
        print(f"Averaged value: Nu_bar={mean(Nu)}, q1_bar={mean(q1)}, q2_bar={mean(q2)}, q3_bar={mean(q3)}")

    a = 12.12268499
    b = -9.528231534

    return a * mean(q3) + b

bcs = {
        "HT": ((730, 9e6), (500, 22.5e6)),
        "LT": ((630, 9e6), (400, 22.5e6)),
    }

media = "CO2"

g = 9.80665
d = 2e-3 # diameter
d_hyd = 1.27e-3
G = 498.341

for name, bc in bcs.items():
    ((T_h, p_h),(T_c, p_c)) = bc
    C1 = predict_C1(T_h, p_h, T_c, p_c, G, d, "CO2", verbose=True)

    print(f'Predicted C1={C1}')