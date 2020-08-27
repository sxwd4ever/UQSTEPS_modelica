'''
    system level models for test cases
'''

import numpy as np
import math
from enum import IntEnum
from CoolProp.CoolProp import PropsSI
from physics import Temperature, Pressure, MDot, Velocity, Density, Angle

class StreamType(IntEnum):
    Hot = 0,
    Cold = 1

class DesignParam(object):
    '''
    Structure containing On Design parameters
    '''

    def __init__(self, d_c = 12e-3, p=[9e6, 20e6], T=[500, 300], mdot = [8.3, 8.3], Re = 2000, N_seg = 10, len_seg=12e-3, pitch = 12e-3, phi= math.pi / 4):
        # On design params initialization 
        self.d_c = d_c
        self.p = np.array(p)
        self.T = np.array(T)
        self.mdot = np.array(mdot)
        self.Re = np.array(Re)
        self.N_seg = N_seg
        self.len_seg = len_seg  
        self.pitch = pitch
        self.phi = phi

        # parameters determined in cal_on_design_params
        self.d_h = 0.0
        self.A_c = 0.0
        self.A_f = 0.0
        self.N_ch = 0.0
        self.mu = np.zeros(len(self.mdot))  
        self.G = np.zeros(len(self.mdot)) # mass flux kg/(m^2 * s) should be constant if Re, P, T, d_c are constant

        self.cal_design_params()

    def cal_geo_params(self):
        d_c = self.d_c
        A_c = math.pi * d_c * d_c / 8 
        peri_c = d_c * math.pi / 2 + d_c          
        d_h = 4 * A_c / peri_c

        return (A_c, d_h, peri_c)

    def area_flow(self):
        return self.A_c * self.N_ch

    def cal_design_params(self):
        
        #  0 = hot, 1 = cold, following Enum PipeType
        A_fx = np.zeros(len(StreamType))

        (self.A_c, self.d_h, peri_c) = self.cal_geo_params()

        for i in StreamType:
            self.mu[i] = CP.PropsSI('V', 'P', self.p[i], 'T', self.T[i], "CO2")
            
        A_fx = np.divide(self.mdot, self.mu) * (self.d_h / self.Re)   

        self.A_f = max(A_fx)
        self.N_ch = math.ceil(self.A_f / self.A_c)
        self.G = self.mdot / self.A_f

class OffDesignParam(object):

    def __init__(self, param_des: DesignParam, mdot = [], G = []):

        if mdot == []:
            if G == []:
                raise ValueError('no mdot or G assigned for off design')
            else:
                mdot = np.array(G) * param_des.area_flow()     

        num_stream = len(StreamType)

        self.mdot = np.array(mdot)
        self.Re = np.zeros(num_stream)
        self.G = np.zeros(num_stream)
        self.param_des = param_des
        self.__update()

    def __update(self):
        p_des = self.param_des
        A_f = self.param_des.area_flow() 
        self.G = self.mdot / A_f
        self.Re = np.divide(self.mdot, p_des.mu) * p_des.d_h / A_f

class ThermoState(object):
    '''
    ThermoState of the fluid
    '''
    def __init__(self, p : Pressure = Pressure(90, 'bar'), T = Temperature(730)):
        super().__init__()
        self._p = p
        self._T = T

    @property
    def p(self):
        """The p property."""
        return self._p
    @p.setter
    def p(self, value):
        self._p = value

    @property
    def T(self):
        """The T property."""
        return self._T
    @T.setter
    def T(self, value):
        self._T = value
        

class TestConfig(object):
    '''
        configuration of diffrent test cases
        this class is designed in a flatten way to make the configuration and debug simplier
    '''

    def __init__(self):
        super().__init__()
        # default values of test case - Meshram [2016] Fig 4(b) - zigzag channel @ High Temperature

        # document property         
        self.name = "zigzg channal @ High Temperature"

        # geomerty
        self.d_c = 2e-3
        self.pitch = 12.3e-3
        self.phi = uc_from_deg(180 - 108) / 2

        # on design propery
        self.Re_des = 14500
        self.N_seg = 10
        self.l_cell = 12e-3

        # fluid property
        # hot stream state 
        self.st_hot_in = ThermoState(p=Pressure.from_bar(90), T = Temperature(730))
        self.st_hot_out = ThermoState(p=Pressure.from_bar(90), T = Temperature(576.69))

        # cold stream state
        self.st_cold_in = ThermoState(p = Pressure.from_bar(225), T = Temperature(500))
        self.st_cold_out = ThermoState(p = Pressure.from_bar(225), T = Temperature(639.15))
        
        # mass flow rate
        self.mdot_hot = MDot(10)
        self.mdot_cold = MDot(10)

        self.media_hot = "Hot"
        self.media_cold = "Cold"

        # off design configuration
        self.mdot_hot_odes = MDot(100)
        self.mdot_cold_odes = MDot(100)

        self.u_hot_odes = Velocity(7.564)
        self.u_cold_odes = Velocity(1.876)

        st = self.st_hot_in    
        self.rho_hot_odes = Density(PropsSI('D', 'P', st.p, 'T', st.T, media_hot))

        st = self.st_cold_in
        self.rho_cold_odes = Density(PropsSI('D', 'P', st.p, 'T', st.T, media_cold))

