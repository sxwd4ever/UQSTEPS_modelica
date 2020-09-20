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
from math import sin

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

# result = CP.CoolProp.PropsSI("H", "P", 20e6, 'T', 700 + 273.15, 'CO2')

# print("result=" + str(result))

AS_SAT = AbstractState("HEOS", "CO2")
AS_SAT.update(CoolProp.PT_INPUTS, 20e6, 273.15+700)
print("First saturation derivative:", AS_SAT.hmass(), "Pa/K")