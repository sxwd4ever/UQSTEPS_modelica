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
import json


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
       
class Variable(object):
    '''
    definition of a variable in simulation
    simulation result is a collection of variables 
    '''
    def __init__(self, key, unit='[1]', val=0.0, text=None):
        self.key = key
        self.unit = unit
        self.val = val
        self.text = key if text == None else text

    def __str__(self):
        return f'{self.val}'

    def __repr__(self):
        return f'{self.text}={self.val} [{self.unit}]'

class TestResult(dict):
    """
    docstring for TestResult.

    TestResult contains Result for one test (TestBatch.TestGroup.Test)

    The result can be saved into an excel file or loaded accordingly. 
    Test config that generates each individual result will be saved as well
    """
    def __init__(self, name, results:Dict[str, Variable]):
        self._idx = 0        
        self.name = name
        self.results = results 
    
    @classmethod
    def from_keyvalues(cls, sol_dict:Dict[str, Variable], values:list):
        """
        :name of this Test Run (parameter sweep)
        """
        keys = list(sol_dict.keys())
        d = OrderedDict() # a copy of values

       # collect data in solutions
        for i in range(len(keys)):
            key = keys[i]
            sol = values[i]
            var = sol_dict[key]
            if not sol is None:
                # copy the sol_dict
                # for now, save the last val, which is the value when system
                # achive equilibrium
                d[key] = Variable(var.key, var.unit, sol[-1], var.text)

        return cls("", d)
   
    def __iter__(self):        
        self._idx = 0
        return self

    def __next__(self) -> str:
        result = None    

        if self._idx < len(self.results):
            result = list(self.results.values())[self._idx]
            self._idx += 1
            return result
        else:
            raise StopIteration   

    def to_dict(self) -> dict:
        d = OrderedDict()
        for var_ in self.results.values():
            d[var_.key] = var_
        return d

class TestConfig(dict):
    def __init__(self, name, params:dict):
        self.name = name
        self.params = params
        self._idx = 0
    
    @classmethod
    def from_json(cls, name, json_str):
        json_dict = json.loads(json_str)
        return TestConfig(name, **json_dict)

    def to_json(self):
        return json.dumps(self.params)

    def __iter__(self):
        self._idx = 0
        return self

    def __next__(self) -> object:
        param = None

        if self._idx < len(self.params):
            param = list(self.params.values())[self._idx]
            self._idx += 1
        else:
            raise StopIteration

        return param        

    def gen_params(self) -> list:
        params = []
        for (k, v) in self.params.items():
            params.append("{0}={1}".format(k, str(v)))
        
        return params
    
    def to_dict(self) -> dict:
        return deepcopy(self.params)    

class TestItem:
    def __init__(self, name, cfg:dict, result:dict):
        self.name = name
        self.cfg = TestConfig(name, cfg)
        self.result = TestResult(name, result)

    @classmethod
    def from_json(cls, name, json_str):
        json_dict = json.loads(json_str)
        return cls(name, **json_dict)

    def to_dict(self):
        d = {}
        d['cfg'] = self.cfg.to_dict()
        d['result'] = self.result.to_dict()
        return d

class TestGroup(dict):
    def __init__(self, name:str, tests:dict):
        self.name = name
        self.tests = OrderedDict()
        self._idx = 0
        for (k,v) in tests.items():
            self.tests[k] = TestItem(k, **v)

    @classmethod
    def from_json(cls, name, json_str):
        json_dict = json.loads(json_str)
        return cls(name, json_dict)
    
    def to_dict(self) -> dict:
        d = {}
        for (k, v) in self.tests.items():
            d[k] = v.to_dict()
        return d

    def __iter__(self):
        self._idx = 0
        return self

    def __next__(self) -> TestItem:

        item = None

        if self._idx < len(self.tests):
            item = list(self.tests.values())[self._idx]
            self._idx += 1
        else:
            raise StopIteration

        return item
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
    def __init__(self, name, groups:dict):
        self.name = name
        self.groups = OrderedDict()
        self._idx = 0

        for (k,v) in groups.items():
            self.groups[k] = TestGroup(k, v)
        # map(lambda (k,v): TestItem.from_json(k, v), tests)

    def __iter__(self):
        self._idx = 0
        return self

    def __next__(self) -> TestGroup:

        item = None
        if self._idx < len(self.groups):
            item = list(self.groups.values())[self._idx]
            self._idx += 1
        else:
            raise StopIteration

        return item

    @classmethod
    def from_json(cls, name, json_str):
        json_dict = json.loads(json_str)
        return cls(name, json_dict)

    def to_dict(self):
        d = {}
        for (k, v) in self.groups.items():
            d[k] = v.to_dict()
        return d

    def to_json(self):
        d = self.to_dict()
        return json.dumps(d, indent=2)


col_start = TestConstants.DATA_FILE_COL_START

def save_test(ds_test: TestDataSet):
    import xlwings as xw

    # set visible=True for debug purpose
    app = xw.App(visible=True)
    result_len = 10
    wbs = {}

    try:
        for group in ds_test:
            # prepare workbook and sheet
            gname = group.name
            if gname in wbs:
                wbk = wbs[gname]
                wbk.activate()                
            else:
                wbk = app.books.add()
                wbs[gname] = wbk

            for test in group:
                result = test.result
                cfg = test.cfg                
                test_name = test.name

                # generate mock results
                # for para in result:
                    # result.set_result(para, np.random.rand(1, result_len))

                sht =  wbk.sheets.add(name = test.name[0:30]) # max lenth for a sheet name
                # sht =  wbk.sheets.add()
                (u, l) = (1, col_start)

                ex_helper = ExcelHelper(sht)
                
                (b, r) = ex_helper.write_table(cfg.params, title={"table 1": "Config"}, up=u, left=l, linespacing=True)

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
