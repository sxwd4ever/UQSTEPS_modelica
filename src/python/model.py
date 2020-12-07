'''
    system level models for test cases
'''

import math
from copy import copy, deepcopy
from datetime import datetime as dt
from enum import IntEnum
from typing import Dict, Tuple, overload
from collections import OrderedDict
from collections.abc import MutableMapping

# 3rd or python modules
import numpy as np
from numpy.lib.utils import deprecate
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

class MyDict(dict):
    
    @classmethod
    def from_json(cls, name, json_str):
        json_dict = json.loads(json_str)
        return cls(name, **json_dict)

    def to_json(self):
        return json.dumps(self)
    
    def to_dict(self) -> dict:
        return self   

class Variable(object):
    '''
    definition of a variable in simulation
    simulation result is a collection of variables 
    '''
    def __init__(self, key, unit='[1]', val=0.0, text=None):
        '''
            key: key in modelica's solution
        '''
        self.key = key
        self.unit = unit
        self.val = val
        self.text = key if text == None else text

    def __str__(self):
        return f'{self.val}'

    def __repr__(self):
        return f'{self.text}={self.val} [{self.unit}]'

class TestResult(MyDict):
    """
    docstring for TestResult.

    TestResult contains Result for one test (TestBatch.TestGroup.Test)

    The result can be saved into an excel file or loaded accordingly. 
    Test config that generates each individual result will be saved as well
    """
    def __init__(self, name, results:Dict[str, Variable]):
        super().__init__(results)        
        self.name = name        
    
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
            v = Variable(var.key, var.unit, -1 , var.text)
            if not sol is None and len(sol) > 0:
                # copy the sol_dict
                # for now, save the last val, which is the value when system
                # achive equilibrium
                v.val = sol[-1]
            d[key] = v

        return cls("", d)
   
    # def __iter__(self):        
    #     self._idx = 0
    #     return self

    # def __next__(self) -> Variable:
    #     result = None    

    #     if self._idx < len(self.results):
    #         result = list(self.results.values())[self._idx]
    #         self._idx += 1
    #         return result
    #     else:
    #         raise StopIteration   

    # def to_dict(self) -> dict:
    #     d = OrderedDict()
    #     for var_ in self.values():
    #         d[var_.key] = var_
    #     return d
class TestDataView(dict):
    '''
        Data view mapping Test result to other data views
    '''

    def __init__(self, name, map:Dict[str, str]) -> None:
        self.name = name
        self.mapping = map

    def maps(self, data:TestResult) -> Dict[str, Variable]:
        view = OrderedDict()
        results = data
        for (src, dst) in self.mapping.items():
            v_new = Variable(dst)
            view[dst] = v_new

            if src in results.keys():
                v = results[src]
                v_new.val = v.val 
                v_new.unit = v.unit
                v_new.text = v.text               

        return view

    def to_dict(self) -> Dict[str, str]:
        return deepcopy(self.mapping)

class TestConfig(MyDict):
    def __init__(self, name, params:dict):
        super().__init__(params)
        self.name = name        
    
    @classmethod
    def from_json(cls, name, json_str):
        json_dict = json.loads(json_str)
        return cls(name, **json_dict)

    def to_json(self):
        return json.dumps(self)

    def gen_params(self) -> list:
        params = []
        for (k, v) in self.items():
            params.append("{0}={1}".format(k, str(v)))
        
        return params
class TestItem:
    def __init__(self, name, cfg:dict, result:dict):
        self.name = name
        self.cfg = TestConfig(name, cfg)
        self.result = TestResult(name, result)
        self.views : Dict[str, TestDataView] = {}

    def add_view(self, view:TestDataView):
        self.views[view.name] = view

    def has_view(self) -> bool:
        return len(self.views) != 0

    @classmethod
    def from_json(cls, name, json_str):
        json_dict = json.loads(json_str)
        return cls(name, **json_dict)

    def to_dict(self):
        d = {}
        d['cfg'] = self.cfg.to_dict()
        d['result'] = self.result.to_dict()
        views = {}
        for (k, v) in self.views.items():
            views[k] = v.to_dict()
        d['views'] = views
        return d

class TestGroup(MyDict):
    def __init__(self, name:str, tests:dict):
        super().__init__(tests)
        self.name = name

    @classmethod
    def from_json(cls, name, json_str):
        json_dict = json.loads(json_str)
        return cls.from_dict(name, **json_dict)

    @classmethod
    def from_dict(cls, name, **json_dict):
        tests = OrderedDict()
        for (k,v) in json_dict.items():
            tests[k] = TestItem(k, **v)

        return cls(name, tests)    

    def to_dict(self) -> dict:
        d = {}
        for (k, v) in self.items():
            d[k] = v.to_dict()
        return d
   
class TestDataSet(MyDict):
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
        super().__init__(groups)
        self.name = name

    def add_view(self, mapping:list, idx_start = 1, view_name_tmpl="default_view"):
        for g in self.values():
            for t in g.values():
                idx = idx_start
                for m in mapping:                
                    t.add_view(TestDataView(f"{view_name_tmpl} {idx}", m))
                    idx += 1

    def to_dict(self):
        d = {}
        for (k, v) in self.items():
            d[k] = v.to_dict()
        return d

    def to_json(self):
        d = self.to_dict()
        return json.dumps(d, indent=2)

    @classmethod
    def from_json(cls, name, json_str):
        json_dict = json.loads(json_str)

        groups = OrderedDict()        
        for (k,v) in json_dict.items():
            groups[k] = TestGroup.from_dict(k, **v)

        return cls(name, groups)

    KEY_NAME = 'KEY'
    KEY_NA = '_' # N/A, for space occupation

    def dicts2table(self, dicts:list, titles:list, titles_src:list, key_name:str = KEY_NAME, value_name:str = None) -> dict:
        ''' 
        use lists of dicts to form a 2D table
        '''

        table = {}
        table['Key'] = []
        for i in range(0, len(titles)):
            table['Key'].append(titles[i])

        for dict_ in dicts:
            name_dict = dict_[key_name]
            table['Key'].append(name_dict)

            for (k, v) in dict_.items():
                if k == key_name:
                    continue
                
                if not k in table.keys():
                    # add a new row
                    table[k] = [] 
                    # and headings for each row
                    for i in range(0, len(titles)):
                        src_name = titles_src[i]
                        val = k # use k of this item as default value
                        if src_name == self.KEY_NA:
                            val = src_name
                        elif src_name == None:
                            val = v
                        elif src_name != 'KEY':
                            val = getattr(v, src_name)
                            
                        table[k].append(val)

                # append the value
                if(value_name != None):
                    table[k].append(getattr(v, value_name))
                else:
                    table[k].append(v)

        return table

    def get_cfgs(self, gname:str) -> dict:
        '''
        get all cfgs as a 2D table
        '''
        if not gname in self.keys():
            raise KeyError(f'no test config with group name={gname}')

        cfgs = []
        for (tname, test) in self[gname].items():
            cfg = test.cfg.to_dict()
            cfg[self.KEY_NAME] = tname
            cfgs.append(cfg)

        return self.dicts2table(cfgs, ['Name', 'Unit', 'Description'], ['KEY', self.KEY_NA, self.KEY_NA]) 

    def set_cfgs(self, gname:str, cfgs:dict):
        ''' 
        use a 2D table to set all the configurations
        '''

        start = 3

        for i in range(start, len(cfgs['Key'])):
            tname = cfgs['Key'][i]
            dict_ = {}
            for (k, v) in cfgs.items():
                if( k == 'Key'):
                    continue
                dict_[k] = v[i]

            cfg:TestConfig = TestConfig(tname,dict_)

            if not gname in self.keys():
                self[gname] = TestGroup(gname, {})

            g:TestGroup = self[gname]

            if not tname in g.keys():
                g[tname] = TestItem(tname, {}, {})

            test:TestItem = g[tname]
            test.cfg = cfg    

    def get_results(self, gname:str) -> dict:
        '''
        get all results as a 2D table
        '''
        if not gname in self.keys():
            raise KeyError(f'no test result with group name={gname}')

        results = []
        for (tname, test) in self[gname].items():
            result = test.result.to_dict()
            result[self.KEY_NAME] = tname
            results.append(result)     

        return self.dicts2table(results, ['Name', 'Unit', 'Description'], ['text', 'unit', self.KEY_NA], value_name='val')   

    def set_results(self, gname:str, results:dict):
        ''' 
        use a 2D table to set all the test values 
        '''

        start = 3

        for i in range(start, len(results['Key'])):
            tname = results['Key'][i]
            dict_ = {}
            for (k, v) in results.items():
                if( k == 'Key'):
                    continue
                dict_[k] = Variable(k, v[1], v[i], v[0])

            result:TestResult = TestResult(tname,dict_)

            if not gname in self.keys():
                self[gname] = TestGroup(gname, {})

            g:TestGroup = self[gname]

            if not tname in g.keys():
                g[tname] = TestItem(tname, {}, {})

            test:TestItem = g[tname]
            test.result = result

    @classmethod
    def gen_cfg(cls, cfg_ref:dict, cfg_offset:dict, keys:list, idx:int, trace:list, result:dict, gname_template:str, cfg_name:str):
        '''
        use iteration method to walk through all the combinations of configuations
        '''
        key = keys[idx]

        node = cfg_offset[key]

        if isinstance(node, dict): # for a dict node, run recurrsion from this point 
            param_keys = node['keys']
            # for each row in the dict (except for key row)
            # use recursion to create a variation
            for (k, v) in node.items(): 
                if k == 'keys':
                    continue
                else:
                    trace_cp = deepcopy(trace) # local copy
                    cfg_name_cp = copy(cfg_name)
                    cfg_name_cp += f'{k},'

                    for i in range(0, len(param_keys)):
                        trace_cp.append((param_keys[i],v[i]))

                    if idx == len(keys) - 1:
                        cls.create_cfg(cfg_ref, keys, trace_cp, result, gname_template, cfg_name_cp)
                    else:
                        cls.gen_cfg(cfg_ref, cfg_offset, keys, idx + 1, trace_cp, result, gname_template, cfg_name_cp)            

        else: # for list, each value in the node's list will create a variation   
            for val in node:
                trace_cp = deepcopy(trace) # local copy
                trace_cp.append((key, val))

                cfg_name_cp = copy(cfg_name)                
                cfg_name_cp += f'{key}={val},'

                if idx == len(keys) - 1:
                    cls.create_cfg(cfg_ref, keys, trace_cp, result, gname_template, cfg_name_cp)
                else:
                    # iterate for deeper level
                    cls.gen_cfg(cfg_ref, cfg_offset, keys, idx + 1, trace_cp, result, gname_template, cfg_name_cp)

    @classmethod    
    def create_cfg(cls, cfg_ref, keys:list, trace:list, result:dict, gname_template:str, cfg_name:str):
        # reach the deepest level, create cfg instance by cloning
        cfg = deepcopy(cfg_ref)  

        (k0, v0) = trace[0] 

        gname = gname_template
        try:
            gname = gname_template % v0 # use first level's key as group name
        except:
            pass

        if not gname in result.keys():
            result[gname]= {}

        root = result[gname]

        for (k,v) in trace: # use variables above level 1 to create item name
            cfg[k] = v

        item = {}
        item["cfg"] = cfg
        item['result'] = { }
        root[cfg_name[0:-1]] = item   

    @classmethod
    def gen_batch_cfg(cls, cfg_ref:dict, cfg_offset:dict, gname_template:str, ds_name = "demoTest"):
        '''
        generate batch cfgs for parameter sweep based on a referenced base config    
        '''

        keys = []
        for k in cfg_offset.keys():
            v = cfg_offset[k]
            if isinstance(v, list) or isinstance(v, dict):
                keys.append(k)

        # keys = ["mdot_main", "T_heater_hot", "T_cooler_cold"]
        trace=[]
        batch_cfg={}    
        cls.gen_cfg(cfg_ref, cfg_offset, keys, 0, trace, batch_cfg, gname_template, "")

        json_cfg_batch = json.dumps(batch_cfg,indent=2)

        ds_test = TestDataSet.from_json(ds_name, json_cfg_batch)
        # print(json_cfg_batch)
        return ds_test   

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
