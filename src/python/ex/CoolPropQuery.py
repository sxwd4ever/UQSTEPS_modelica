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

kc_cf_z_set = [10, 12, 15]    
kc_cf_s_set = [1.0, 1.2, 1.5]

import itertools

for x in itertools.product(kc_cf_s_set, kc_cf_z_set):
        print(x)

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