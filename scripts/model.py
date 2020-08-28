'''
    system level models for test cases
'''

import numpy as np
import math
from enum import IntEnum
from CoolProp.CoolProp import PropsSI
from physics import Temperature, Pressure, MDot, Velocity, Density, Angle, Length
from copy import deepcopy

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
    """
        Corresponding to a set of Tests to see the effect of param sweep
    """
    def __init__(self, name, test_names, enable = True):
        super().__init__()
        self._name = name
        self.test_names = test_names
        self.para_seqs = {}
        self.enable = enable

    def add_para_seq(self, param_name, seq):
        if(len(seq) == len(self.test_names)):
            self.para_seqs[param_name] = seq
        else:
            raise ValueError("number of vals not equal to test number")
            
    def get_para_value(self, idx):
        '''
        get param values for a test
        '''

        para_vals = []

        for name, seq in self.para_seqs.items():
            para_vals.append((name, seq[idx]))

        return para_vals
class ParamSweepSet(object):
    """
    Data structure of ParamSweepSet

    ParaSweepSet
    |
    |---Group 1
    |   |
    |   |---Test 1
        |   |--- Param 1 : val[t_0], val[t_1], ..., val[t_L]
        |   |--- Param 2 :
        |    ...
        |   |--- Param K
        |---Test 2
        ...
        |---Test J
    ...
    |--- Group I

    A Para Sweep Set contains several param Groups for a batch of simulations
    each param Group contains a set of varied params for diferent Tests. 
    each Test, may alter a set of params while keeping rest params constant, is used
    to see the effect these params and will collect values of interested Params 
    and each Param contains a sequenced, time-based values of this param within one test
    """

    def __init__(self):
        super().__init__()

        self.groups = {}
        self.cur_cfg = None

    def new_group(self, name, test_names):
        g = ParamGroup(name, test_names)
        self.groups[name] = g
        return g

    # def gen_configs(self, cfg_cfg: TestConfig):
    #     '''
    #     generate a batch of configs based a given config
    #     configs will be organized in a dict
    #     '''

    #     clone_dict = {}

    #     for k,group in self.groups():

    #         clone_dict[group] = {}

    #         test_names = group.test_names
            
    #         for i in range(0, len(test_names)):
    #             clone = deepcopy(cfg_cfg)

    #             para_value = group.get_para_value(i)

    #             clone.name = clone.name + '@' + test_names[i]

    #             for name, val in para_value:
    #                 clone.__setattr__(name, val)

    #             clone_dict[group][clone.name] = clone

    #     return clone_dict

class TestResult(object):
    """
    docstring for TestResult.

    TestResult contains Result for one test (TestBatch.TestGroup.Test)

    The result can be saved into an excel file or loaded accordingly. 
    Test config that generates each individual result will be saved as well
    """
    def __init__(self, name, para_design : DesignParam):
        """
        :name of this Test Run (parameter sweep)
        """
        super(TestResult, self).__init__()
        
        self.name = name

        self._result_dict = {}

        self.__init_result_dict(para_design)

        self._idx = 0
        self._keys = []
    
    def __init_result_dict(self, para_design : DesignParam):
                # dict to store solutions under 
        result_dict = {}

        # for all the component's inlet/outlet that needs to be analyzed.
        ports_keys = [
            'source_hot.outlet', 
            'source_cold.outlet', 
            'sink_hot.inlet', 
            'sink_cold.inlet', 
            'pchx.inlet_hot',
            'pchx.outlet_hot',
            'pchx.inlet_cool',
            'pchx.outlet_cool'        
            ]
        
        # common variables for these port
        var_keys = ['h_outflow', 'p', 'm_flow']

        # fill in keys of result dict with above 
        for p_key in ports_keys:
            for val_key in var_keys:
                sol_key = p_key + '.' + val_key
                result_dict[sol_key] = []

        # for each HX cell in PCHE
        N_seg = para_design.N_seg

        var_keys = ['T', 'p', 'h', 'u', 'k', 'rho', 'mu' ,'dp', 'G', 'Re', 'Nu', 'f', 'Q']
        node_keys = ['cell_cold', 'cell_hot']
        for var_key in var_keys:
            for node_key in node_keys:
                for i in range(1, N_seg + 1):        
                    result_dict["pchx.{node}[{idx}].{var}".format(idx = i, var = var_key, node = node_key)] = []

        # some special varible in the model, specify them explicitly like
        # input parameters
        result_dict['pchx.N_seg'] = []  
        result_dict['pchx.length_cell'] = [] 
        result_dict['pchx.phi'] = []  
        result_dict['pchx.d_c'] = []  
        result_dict['pchx.pitch'] = []  
        result_dict['pchx.kim_cor.a'] = []  
        result_dict['pchx.kim_cor.b'] = [] 
        result_dict['pchx.kim_cor.c'] = []  
        result_dict['pchx.kim_cor.d'] = []          

        # calculated initial parameters
        result_dict['pchx.d_h'] = []  
        result_dict['pchx.peri_c'] = [] 
        result_dict['pchx.t_wall'] = []  
        result_dict['pchx.N_ch'] = [] 
        result_dict['pchx.A_c'] = []  
        result_dict['pchx.A_flow'] = [] 
 
        result_dict['pchx.A_stack'] = []                  
        result_dict['pchx.length_ch'] = []     
        result_dict['pchx.Re_hot_start'] = [] 
        result_dict['pchx.Re_cold_start'] = []                           
        result_dict['mdot_hot'] = []  
        result_dict['mdot_cold'] = []  

        self._result_dict = result_dict

    def __cal_row(self, vals = [], eval_str = 'x + y'):

        real_eval_str = eval_str.replace('x','op1').replace('y','op2').replace('z', 'op3')

        r = []
        
        for i in range(0, len(vals[0])):
            op1,op2 = vals[0][i], vals[1][i]

            if(len(vals) == 3):
                op3 = vals[2][i]
            r.append(eval(real_eval_str))

        return r

    def __gen_cal_map(self):
        '''
        map containing all the calculated fields, which are calculated 
        according to the solutions of model
        '''
        # map for calculated values
        cal_map = []
        # recheck efficiency 
        cal_map.append({'key' : 'Q_in', 'vals':['pcm_heater.h_e', 'pcm_heater.h_i', 'pcm_heater.inlet.m_flow'], 'eval_str': 'z * (x - y)'})
        cal_map.append({'key' : 'W_turbine', 'vals':['turbine.h_i', 'turbine.h_ea', 'turbine.inlet.m_flow'], 'eval_str': 'z * (x - y)'})
        cal_map.append({'key' : 'W_comp', 'vals':['pump.h_ea', 'pump.h_i', 'pump.inlet.m_flow'], 'eval_str': 'z * (x - y)'})
        cal_map.append({'key' : 'W_comp_recom', 'vals':['recom_pump.h_ea', 'recom_pump.h_i', 'recom_pump.inlet.m_flow'], 'eval_str': 'z * (x - y)'})
        cal_map.append({'key' : 'W_net', 'vals':['W_turbine', 'W_comp', 'W_comp_recom'], 'eval_str': 'x - y - z'})
        cal_map.append({'key' : 'eta_all_recal', 'vals':['W_net', 'Q_in'], 'eval_str': 'x / y * 100'})
        cal_map.append({'key' : 'eta_all_recal', 'vals':['W_net', 'Q_in'], 'eval_str': 'x / y * 100'})

        # performance
        # LTR
        cal_map.append({'key' : 'LTR_cool_dt', 'vals':['recup_low.crec_out.T','recup_low.crec_in.T'], 'eval_str': 'x - y'})
        cal_map.append({'key' : 'LTR_hot_dt', 'vals':['recup_low.inlet.T','recup_low.outlet.T'], 'eval_str': 'x - y'})
        cal_map.append({'key' : 'LTR_cool_dQ', 'vals':['recup_low.h_cool_e','recup_low.h_cool_i', 'recup_low.crec_in.m_flow'], 'eval_str': '(x - y) * z'})
        cal_map.append({'key' : 'LTR_hot_dQ', 'vals':['recup_low.h_hot_i','recup_low.h_hot_e', 'recup_low.inlet.m_flow'], 'eval_str': '(x - y) * z'})

        # HTR
        cal_map.append({'key' : 'HTR_cool_dt', 'vals':['recup_high.crec_out.T','recup_high.crec_in.T'], 'eval_str': 'x - y'})
        cal_map.append({'key' : 'HTR_hot_dt', 'vals':['recup_high.inlet.T','recup_high.outlet.T'], 'eval_str': 'x - y'})
        cal_map.append({'key' : 'HTR_cool_dQ', 'vals':['recup_high.h_cool_e','recup_high.h_cool_i', 'recup_high.crec_in.m_flow'], 'eval_str': '(x - y) * z'})
        cal_map.append({'key' : 'HTR_hot_dQ', 'vals':['recup_high.h_hot_i','recup_high.h_hot_e', 'recup_high.inlet.m_flow'], 'eval_str': '(x - y) * z'})

        return cal_map

    def update_cal_columns(self, result_dict):
        # fill calculated fields into the result_dict
        cal_map = self.__gen_cal_map()
        result_dict = self._result_dict
        
        # fill in the calculated value by python's eval
        for item in cal_map:
            vals = item['vals']

            data = result_dict

            if len(vals) == 2:
                result_dict[item['key']] = __cal_row([data[vals[0]], data[vals[1]]], eval_str=item['eval_str'])
            elif len(vals) == 3:
                result_dict[item['key']] = __cal_row([data[vals[0]], data[vals[1]], data[vals[2]]], eval_str=item['eval_str'])   

    def __iter__(self):
        self._keys = list(self._result_dict.keys())
        self._idx = 0

        return self

    def __next__(self) -> str:

        if self._idx < len(self._keys):
            key = self._keys[self._idx]
            self._idx += 1
            return key
        else:
            raise StopIteration

    def set_result(self, para, val_seq):

        self._result_dict[para] = val_seq

    def get_values(self, para_name):

        return self._result_dict[para_name]

class TestDataItem(object):

    def __init__(self, cfg : TestConfig):

        self._name = cfg.name

        self.cfg = cfg

        self._para_des, self._para_odes = cfg.gen_test_param()

        self.result = TestResult(cfg.name, self._para_des)

    @property
    def para_des(self):
        """The para_des property."""
        return self._para_des

    @property
    def para_odes(self):
        """The para_odes property."""
        return self._para_odes

    def gen_sim_param(self):
        # set up on design parameters
        (hot, cold) = (StreamType.Hot, StreamType.Cold)
        para_dict = {}

        para = self.para_des
        para_dict["N_ch"] = para.N_ch

        para_dict["T_cold_in"] = para.T[cold]
        para_dict["p_cold_in"] = para.p[cold]

        para_dict["T_hot_in"] = para.T[hot]
        para_dict["p_hot_in"] = para.p[hot]

        para_dict["phi"] = para.phi
        para_dict["pitch"] = para.pitch
        para_dict["d_c"] = para.d_c
        para_dict["length_cell"] = para.len_seg       

        # set up off design parameters
        para = self.para_odes
        para_dict["Re_hot_start"] = para.Re[hot]
        para_dict["Re_cold_start"] = para.Re[cold]
        
        para_dict["mdot_hot"] = para.mdot[hot]
        para_dict["mdot_cold"] = para.mdot[cold]

        params = []

        for k, v in para_dict.items():
            params.append("{0}={1}".format(k, str(v)))

        return params        

class TestDataSet(object):
    """
    Test Data Set, container for all test configs, result set    
    """

    def __init__(self, cfg : TestConfig = None, para_set : ParamSweepSet = None):
        super().__init__()
        self._para_set =  para_set

        self.cfg_ref = cfg

        self.test_data_set = []

        self.idx = 0

    def __iter__(self):

        self.test_data_set = []
        idx = 0

        cfg_ref = self.cfg_ref

        for k, group in self._para_set.groups.items():

            if not group.enable:
                # ignore unable cfg
                continue

            test_names = group.test_names
            
            for i in range(0, len(test_names)):
                clone = deepcopy(cfg_ref)

                para_value = group.get_para_value(i)

                clone.name = clone.name + '@' + test_names[i]

                for name, val in para_value:
                    clone.__setattr__(name, val)
                

                self.test_data_set.append(TestDataItem(cfg=clone))

        return self

    def __next__(self) -> TestDataItem:

        item = None

        if self.idx < len(self.test_data_set):
            item = self.test_data_set[self.idx]
            self.idx += 1
        else:
            raise StopIteration

        return item

    def items(self):  
        return self

def main():

    cfg_ref = TestConfig()

    para_sweep = ParamSweepSet()

    g = para_sweep.new_group("group1", ["Re_des = 5000", "Re_des = 20000"])
    g.add_para_seq("Re_des", [5000, 20000])

    g = para_sweep.new_group("group2", [r"$\dot{m}_{off}$ = 10", r"${\dot{m}_{off}$ = 100"])
    g.add_para_seq("mdot_hot_odes", [MDot(10), MDot(100)])
    g.add_para_seq("mdot_cold_odes", [MDot(10), MDot(100)])

    g = para_sweep.new_group("group3", [r"$(p,T)_{hi}$ = 10 MPa, 450°", r"$(p,T)_{hi}$ = 12 MPa, 300°"])
    g.add_para_seq("p_hot_in", [Pressure.MPa(10), Pressure(12)])
    g.add_para_seq("T_hot_in", [Temperature.degC(730), Temperature.degC(300)])    

    ds = TestDataSet(cfg = cfg_ref, para_set = para_sweep)

    for data_item in ds:
        print(data_item.cfg.name)

    print('done')


if __name__ == "__main__":
    main()