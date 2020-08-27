'''
    system level models for test cases
'''

import numpy as np
import math
from enum import IntEnum
from CoolProp.CoolProp import PropsSI
from physics import Temperature, Pressure, MDot, Velocity, Density, Angle, Length
from copy import deepcopy

import simplejson as json

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
            self.mu[i] = PropsSI('V', 'P', self.p[i], 'T', self.T[i], "CO2")
            
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
    def p(self) -> Pressure:
        """The p property."""
        return self._p
    @p.setter
    def p(self, value: Pressure):
        self._p = value

    @property
    def T(self) -> Temperature:
        """The T property."""
        return self._T
    @T.setter
    def T(self, value: Temperature):
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
        self.name = "zigzg channal, High Temperature"

        # geomerty
        self.d_c = Length.mm(2)

        self.pitch = Length.mm(12.3)

        self.phi = Angle.deg((180 - 108) / 2)

        # on design propery
        self.Re_des = 14500
        self.N_seg = 10
        self.l_cell = Length.mm(12)

        # fluid property
        # hot stream state 
        self.p_hot_in = p=Pressure.bar(90)
        self.T_hot_in = T = Temperature(730)
        self.p_hot_out = p=Pressure.bar(90)
        self.T_hot_out = T = Temperature(576.69)

        # cold stream state
        self.p_cold_in = p=Pressure.bar(225)
        self.T_cold_in = T = Temperature(500)
        self.p_cold_out = p=Pressure.bar(90)
        self.T_cold_out = T = Temperature(639.15)
        
        self.media_hot = "CO2"
        self.media_cold = "CO2"
        
        # mass flow rate
        self.mdot_hot = MDot(10)
        self.mdot_cold = MDot(10)

        # off design configuration
        self.mdot_hot_odes = MDot(100)
        self.mdot_cold_odes = MDot(100)

        self.u_hot_odes = Velocity(7.564)
        self.u_cold_odes = Velocity(1.876)

        self.rho_hot_odes = Density(PropsSI('D', 'P', self.p_hot_in.mag, 'T', self.T_hot_in.mag, self.media_hot))

        self.rho_cold_odes = Density(PropsSI('D', 'P', self.p_cold_in.mag, 'T', self.T_cold_in.mag, self.media_cold))

    def gen_test_param(self) -> (DesignParam, OffDesignParam):

        param_des = DesignParam(
            d_c=self.d_c.mag, 
            p=[self.p_hot_in.mag, self.p_cold_in.mag], 
            T=[self.T_hot_in.mag, self.T_cold_in.mag], 
            Re=self.Re_des, 
            mdot=[self.mdot_hot.mag, self.mdot_cold.mag],
            N_seg= self.N_seg, 
            len_seg= self.l_cell.mag, 
            pitch = self.pitch.mag, 
            phi = self.phi.mag)

        # param_odes = OffDesignParam(param_des=param_des, G = np.multiply(u_odes, rho_odes))
        param_odes = OffDesignParam(param_des=param_des, mdot = [self.mdot_hot_odes.mag, self.mdot_cold_odes.mag])

        return (param_des, param_odes)

class ParamGroup(object):

    def __init__(self, test_names):
        super().__init__()
        self.test_names = test_names
        self.para_seqs = {}

    def add_para_seq(self, param_name, seq):
        if(len(seq) == len(self.test_names)):
            self.para_seqs[param_name] = seq
        else:
            raise ValueError("number of vals not equal to test number")
            
    def get_test_value(self, idx):
        '''
        get value set for one test
        '''

        test_vals = []

        for name, seq in self.para_seqs.items():
            test_vals.append((name, seq[idx]))

        return test_vals

class ParamSweepSet(object):

    def __init__(self):
        super().__init__()

        self.groups = {}

    def new_group(self, name, test_names):
        g = ParamGroup(test_names)
        self.groups[name] = g
        return g

    def gen_configs(self, config: TestConfig):
        '''
        generate serveral configs based a given config
        '''
        clones = []

        for k,group in self.groups.items():
            test_names = group.test_names
            
            for i in range(0, len(test_names)):
                clone = deepcopy(config)

                test_value = group.get_test_value(i)

                clone.name = clone.name + '@' + test_names[i]

                for name, val in test_value:
                    clone.__setattr__(name, val)

                clones.append(clone)

        return clones


def main():



    print('done')


if __name__ == "__main__":
    main()