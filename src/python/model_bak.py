'''
    system level models for test cases
'''

import math
from copy import deepcopy
from datetime import datetime as dt
from enum import IntEnum
from typing import Dict, Tuple
from collections import OrderedDict
from collections.abc import MutableMapping

# 3rd or python modules
import numpy as np
from CoolProp.CoolProp import PropsSI
import pandas as pd
from xlwings.main import Book


# my modules
from physics import (Quantity, Angle, Density, Length, MDot, Pressure, Temperature,
                     Velocity)

from utils import ExcelHelper

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

        A_fx =  np.divide(self.mdot, self.mu) * (self.d_h / self.Re)   

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
class TestConstants(object):
    """
    Constants and key definition for test
    """

    d_c = "d_c"

    pitch = "pitch"

    phi = "phi"

    # on design propery
    Re_des = "Re_des"
    N_seg = "N_seg"
    l_cell = "l_cell"

    # fluid property
    # hot stream state 
    p_hot_in = "p_hot_in"
    T_hot_in = "T_hot_in"
    p_hot_out = "p_hot_out"
    T_hot_out = "T_hot_out"

    # cold stream state
    p_cold_in = "p_cold_in"
    T_cold_in = "T_cold_in"
    p_cold_out = "p_cold_out"
    T_cold_out = "T_cold_out"

    media_hot = "media_hot"
    media_cold = "media_cold"

    # mass flow rate
    mdot_hot = "mdot_hot"
    mdot_cold = "mdot_cold"

    # off design configuration
    mdot_hot_odes = "mdot_hot_odes"
    mdot_cold_odes = "mdot_cold_odes"

    u_hot_odes = "u_hot_odes"
    u_cold_odes = "u_cold_odes"

    rho_hot_odes = "rho_hot_odes"

    rho_cold_odes = "rho_cold_odes"

    # Format constants
    # column start index in the data storage file
    DATA_FILE_COL_START = 2 

class ParamSet(MutableMapping):

    def __init__(self, *args, **kwargs):
        self.store = dict()
        self.update(dict(*args, **kwargs))  # use the free update to set keys

    def __getitem__(self, key):
        val =  self.store[self.__keytransform__(key)]

        if issubclass(type(val), Quantity):
            val = val.mag

        return val

    def __setitem__(self, key, value):
        self.store[self.__keytransform__(key)] = value

    def __delitem__(self, key):
        del self.store[self.__keytransform__(key)]

    def __iter__(self):
        return iter(self.store)

    def __len__(self):
        return len(self.store)

    def __keytransform__(self, key):
        return key

    def to_dict(self):
        val_dict = {}

        for k in self.store.keys():
            val_dict[k] = self.__getitem__(k)

        return val_dict

    def from_dict(self, dict_ : dict):
        """
        instansiate a config from a dict_
        """

        for (k, v) in dict_.items():
            self.__setitem__(k, v)

    # def to_data_frame(self):
    #     """
    #     transform this config into a pandas data frame 
    #     """
    #     return pd.DataFrame(self.to_dict(), index=[0])               

class TestConfig(ParamSet):
    '''
        configuration of diffrent test cases
        this class is designed in a flatten way to make the configuration and debug simplier
    '''

    def __init__(self, name = "default config", group_name = "default group"):
        super().__init__()
        # default values of the test case - Meshram [2016] Fig 4(b) - zigzag channel @ High Temperature

        # document property  
        self.name = name

        self._group_name = group_name        

        self.__init_config()

    def __init_config(self):

        pk = TestConstants
        # geomerty
        self[pk.d_c] = Length.mm(2)

        self[pk.pitch] = Length.mm(12.3)

        self[pk.phi] = Angle.deg((180 - 108) / 2)

        # on design propery
        self[pk.Re_des] = 14500
        self[pk.N_seg] = 10
        self[pk.l_cell] = Length.mm(12)

        # fluid property
        # hot stream state 
        self[pk.p_hot_in] = Pressure.bar(90)
        self[pk.T_hot_in] = Temperature(730)
        self[pk.p_hot_out] = Pressure.bar(90)
        self[pk.T_hot_out] = Temperature(576.69)

        # cold stream state
        self[pk.p_cold_in] = Pressure.bar(225)
        self[pk.T_cold_in] = Temperature(500)
        self[pk.p_cold_out] = Pressure.bar(90)
        self[pk.T_cold_out] = Temperature(639.15)

        self[pk.media_hot] = "CO2"
        self[pk.media_cold] = "CO2"

        # mass flow rate
        self[pk.mdot_hot] = MDot(100)
        self[pk.mdot_cold] = MDot(100)

        # off design configuration
        self[pk.mdot_hot_odes] = MDot(100)
        self[pk.mdot_cold_odes] = MDot(100)

        self[pk.u_hot_odes] = Velocity(7.564)
        self[pk.u_cold_odes] = Velocity(1.876)

        p, T = self[pk.p_hot_in], self[pk.T_hot_in]
        self[pk.rho_hot_odes] = Density(PropsSI('D','P', p, 'T', p, self[pk.media_hot]))

        p, T = self[pk.p_cold_in], self[pk.T_cold_in]
        self[pk.rho_cold_odes] = Density(PropsSI('D','P', p, 'T', p, self[pk.media_cold]))

    def gen_design_param(self) -> Tuple[DesignParam, OffDesignParam]:
        pk = TestConstants

        param_des:DesignParam = DesignParam(
            d_c=self[pk.d_c], 
            p=[self[pk.p_hot_in], self[pk.p_cold_in]], 
            T=[self[pk.T_hot_in], self[pk.T_cold_in]], 
            Re=self[pk.Re_des], 
            mdot=[self[pk.mdot_hot], self[pk.mdot_cold]],
            N_seg= self[pk.N_seg], 
            len_seg= self[pk.l_cell], 
            pitch = self[pk.pitch], 
            phi = self[pk.phi])

        # param_odes = OffDesignParam(param_des=param_des, G = np.multiply(u_odes, rho_odes))
        param_odes:OffDesignParam = OffDesignParam(param_des=param_des, mdot = [self.get(pk.mdot_hot_odes), self.get(pk.mdot_cold_odes)])
        return (param_des, param_odes)

    @property
    def group_name(self):
        """The group_name property."""
        return self._group_name

    @group_name.setter
    def group_name(self, value):
        self._group_name = value

    @property
    def full_name(self) -> str :
        """The full_name property."""
        return self._group_name + '@' + self.name


class TestParam(ParamSet):

    def __init__(self, name="", title="", **options) -> None:
        super().__init__()
        
        if title=='':
            title = name

        self.name = name
        self.title = name

        for k, v in options.items():
            if k in TestConstants.__dict__.keys():
                self[k] = v
            else:
                raise KeyError('Invalid para name=' + k)

class ParamGroup(object):
    """
        Corresponding to a set of Tests to see the effect of param sweep
    """
    def __init__(self, name, enable = True, ):
        super().__init__()
        self._name = name
        self.paras = OrderedDict()
        self._idx = 0
        
        self._enable = enable

    def add_test_param(self, name="", title="", **options):
        para = TestParam(name=name, title=title, **options)
        self.paras[name] = para
        return para

    # def __iter__(self):
    #     self._idx = 0
    #     return self

    # def __next__(self) -> TestParam:
    #     if self._idx < len(self.paras):
    #         key = self.paras.values()[self._idx]
    #         self._idx += 1
    #         return key
    #     else:
    #         raise StopIteration


    # def add_para_seq(self, param_name, seq):
    #     if(len(seq) == len(self.test_names)):
    #         self.para_seqs[param_name] = seq
    #     else:
    #         raise ValueError("number of vals not equal to test number")
            
    # def get_para_value(self, idx):
    #     '''
    #     get param values for a test
    #     '''

    #     para_vals = []

    #     for name, seq in self.para_seqs.items():
    #         para_vals.append((name, seq[idx]))

    #     return para_vals
    
    @property
    def enable(self):
        """The enable property."""
        return self._enable
    @enable.setter
    def enable(self, value):
        if self._enable == value:
            return

        self._enable = value

    @property
    def name(self):
        """The name of param groupd """
        return self._name
    @name.setter
    def name(self, value):
        self._name = value

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
        N_seg = int(para_design.N_seg)

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
                result_dict[item['key']] = self.__cal_row([data[vals[0]], data[vals[1]]], eval_str=item['eval_str'])
            elif len(vals) == 3:
                result_dict[item['key']] = self.__cal_row([data[vals[0]], data[vals[1]], data[vals[2]]], eval_str=item['eval_str'])   

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

    def to_dict(self) -> dict:
        dict_ = {}
        for k, v in self._result_dict.items():
            dict_[k]  = v.tolist()[0] # get the row 1 for 1d list
            
        return dict_
    def from_dict(self, dict_:dict):
        self._result_dict.clear()

        for k, v in dict_.items():
            self._result_dict[k] = np.array(v)

class TestDataItem(object):

    def __init__(self, cfg : TestConfig):

        self.cfg = cfg

        self._para_des, self._para_odes = cfg.gen_design_param()

        self._result = TestResult(cfg.name, self._para_des)

    @property
    def para_des(self):
        """The para_des property."""
        return self._para_des

    @property
    def para_odes(self):
        """The para_odes property."""
        return self._para_odes

    @property
    def name(self):
        """The name property."""
        return self.cfg.name
    @name.setter
    def name(self, value):
        self.cfg.name = value

    @property
    def result(self) -> TestResult:
        """The result property."""
        return self._result

    @result.setter
    def result(self, value: TestResult):
        self._result = value

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

class TestDataSet(dict):
    """
    TestDataSet
        |
        |---TestDataGroup 1
        |   |
        |   |---TestDataItem for Test 1
                |
                |--- TestConfig
                    |
                    |--- Input Param 1 value
                    |--- Input Param 2 value
                    ...
                    |--- Input Param K
                |--- TestResult
                    |
                    |--- Output Param 1 values: val[t_0], val[t_1], ..., val[t_N]
                    ...
                    |--- Output Param K' values: val[t_0], val[t_1], ..., val[t_N]

            |---TestDataItem for Test 2
            ...
            |---TestDataItem for Test J
        ...
        |--- TestDataGroup I

    A TestDataSet contains several groups test data
    each TestDataGroup contains a group of TestData obtained by applying a set of varied tests params,
    each TestDataItem i contains one TestConfig and one TestResult for Test i.      
    """

    def __init__(self, cfg: TestConfig, name = dt.now().strftime("Test_%Y_%m_%d_%H_%M_%S")):
        super().__init__()
        
        self.name = name
        if cfg == None:
            raise ValueError("None value for test config")

        self.cfg_ref:TestConfig = cfg

        self.test_items = []

        self._idx = 0
        
        self.test_groups = {}

    def add_para_group(self, group : ParamGroup):
        name_group = group.name

        if name_group in self.test_groups:
            print (' Test group with name={0} exsits already'.format(name_group))
            return

        self.test_groups[name_group] = group

        if not group.enable:
            # ignore unable cfg
            return 

        # generate test item for this group of tests

        cfg_ref:TestConfig = self.cfg_ref

        for para in group.paras.values():
            clone = deepcopy(cfg_ref)

            clone.group_name = name_group

            clone.name = para.name

            for k in para:
                clone[k] = para[k]

            self.test_items.append(TestDataItem(cfg=clone))

    def remove_para_group(self, name):

        if name in self.test_groups:
            g = self.test_groups[name]
            to_move = []
            # remove related tests
            for test in self.test_items:
                if test.cfg.group_name == g.name:
                    to_move.append(test)

            for test in to_move:
                self.test_items.remove(test)

            self.test_groups[name] = None
            return 

        return None

    def addTestDataItem(self, item: TestDataItem):
        
        self.test_items.append(item)

    def __iter__(self):
        self._idx = 0
        return self

    def __next__(self) -> TestDataItem:

        item = None

        if self._idx < len(self.test_items):
            item = self.test_items[self._idx]
            self._idx += 1
        else:
            raise StopIteration

        return item

col_start = TestConstants.DATA_FILE_COL_START

def save_test(ds_test: TestDataSet):
    import xlwings as xw

    # set visible=True for debug purpose
    app = xw.App(visible=True)
    result_len = 10
    wbs = {}

    try:
        for test in ds_test:
            # prepare workbook and sheet
            result = test.result
            cfg = test.cfg
            gname = cfg.group_name

            if gname in wbs:
                wbk = wbs[gname]
                wbk.activate()                
            else:
                wbk = app.books.add()
                wbs[gname] = wbk

            # generate mock results
            for para in result:
                result.set_result(para, np.random.rand(1, result_len))

            sht =  wbk.sheets.add(name = cfg.name[0:30]) # max lenth for a sheet name
            # sht =  wbk.sheets.add()
            (u, l) = (1, col_start)

            ex_helper = ExcelHelper(sht)
            
            (b, r) = ex_helper.write_table(test.cfg.to_dict(), title={"table 1": "Config"}, up=u, left=l, linespacing=True)

            (b, r) = ex_helper.write_table(result.to_dict(), title={"table 2": "Result"}, up=b, left = l, linespacing=True)

            # sht:xw.Sheet = wbk.sheets['Sheet1']
            # sht.delete()
            wbk.save( dt.now().strftime("Test_{0}_%Y_%m_%d_%H_%M_%S.xlsx".format(gname)))
            # wbk.close()

    finally:
        app.quit()

def load_test(root_path:str) -> TestDataSet:

    import glob
    import os

    list_of_files: list[str] = glob.glob(os.sep.join([root_path, '*.xlsx'])) # get all files ends with csv
    
    ds_test:TestDataSet = TestDataSet(TestConfig("Default Config"))

    list_of_files = list(filter( lambda x: not x.startswith('.\\~$'), list_of_files))

    if list_of_files == []:
        raise FileNotFoundError("No result file saved in {0}/out".format(root_path))

    latest_file = max(list_of_files, key = os.path.getctime)

    if latest_file is None:
        raise FileNotFoundError("No results file saved in " + root_path)

    import xlwings as xw

    app = xw.App(visible=False)

    try:
        wbk:xw.Book = Book(latest_file)

        for sht in wbk.sheets:
            if sht.name.startswith('Sheet1'):
                continue

            helper = ExcelHelper(sht)
            u = 1
            offset = 2 # + 1 for linespacing, + 1 is start of next table 

            dict_ = helper.read_dict((u, col_start))

            cfg = TestConfig(name=sht.name, group_name=sht.book.name)
            cfg.from_dict(dict_)

            u = u + len(dict_) + offset

            dict_ = helper.read_dict((u, col_start))
            (p_des, p_odes) = cfg.gen_design_param()

            result = TestResult(sht.name, p_des)
            result.from_dict(dict_)

            test_item = TestDataItem(cfg)
            test_item.result = result

            ds_test.addTestDataItem(test_item)

    finally:
        app.quit()

    return ds_test

def main():

    # wbk.name='my book' 
    cfg_ref = TestConfig()

    ds_test = TestDataSet(cfg = cfg_ref)

    g = ParamGroup("zigzag_HT_High_RE")
    g.add_test_param("Re_des = 20000", Re_des=20000)
    ds_test.add_para_group(g)

    g = ParamGroup("zigzag_HT_diff_Re")
    g.add_test_param("Re_des = 5000", Re_des=5000)
    g.add_test_param("Re_des = 20000", Re_des=20000)
    ds_test.add_para_group(g)

    # g = ds_test.new_para_group("zigzag_HT_diff_mdot", [r"$\dot{m}_{off}$ = 10"])
    g = ParamGroup("zigzag_HT_mdot_10")
    g.add_test_param("MDot_10", mdot_hot_odes=MDot(10), mdot_cold_odes=MDot(10))
    ds_test.add_para_group(g)    


    g = ParamGroup("zigzag_HT_diff_mdot")
    g.add_test_param("MDot_10", title=r"$\dot{m}_{off}$ = 10", mdot_hot_odes=MDot(10), mdot_cold_odes=MDot(10))
    g.add_test_param("MDot_100", title=r"${\dot{m}_{off}$ = 100", mdot_hot_odes=MDot(100), mdot_cold_odes=MDot(100))
    ds_test.add_para_group(g)

    # g.submit()
    # g.withdraw()

    g = ParamGroup("zigzag_HT_diff_pT_in")
    g.add_test_param("pT=10_450", title=r"$(p,T)_{hi}$ = 10 MPa, 450°", p_hot_in=Pressure.MPa(10), T_hot_in=Temperature.degC(730))
    g.add_test_param("pT=12_300", title=r"$(p,T)_{hi}$ = 12 MPa, 300°", p_hot_in=Pressure.MPa(12), T_hot_in=Temperature.degC(300))
    ds_test.add_para_group(g)

    save_test(ds_test)

    load_test(".")
   
    print('done')


if __name__ == "__main__":
    main()
