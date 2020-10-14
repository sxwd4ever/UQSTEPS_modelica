# -*- coding: utf-8 -*-
"""
Created on Sun Aug 16 17:48:28 2020

@author: uqxsui
"""

import CoolProp as CP

# P, h -> T ([Pa], [J/kg] -> [degC])
# result = CP.CoolProp.PropsSI("T", "P", 20e6, 'H', 995164, 'CO2') - 273.15

# P, T -> h ([Pa], [K] -> [J/kg])
# result = CP.CoolProp.PropsSI("H", "P", 20e6, 'T', 715 + 273.15, 'CO2')

# P, T -> mu (local viscosity)([Pa], [K] -> [mu])
result = CP.CoolProp.PropsSI("V", "P", 3616780, 'T', 700 + 273.15, 'CO2')


print("result=" + str(result))