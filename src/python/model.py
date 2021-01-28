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
class TestConstants(object):
    """
    Constants and key definition for test
    """
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
                # achive equilibrium state
                v.val = sol[-1]
            d[key] = v

        return cls("", d)
   
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
        self.post_data = {}

    def add_view(self, view:TestDataView):
        self.views[view.name] = view

    def has_view(self) -> bool:
        return len(self.views) != 0
    
    def set_post_data(self, key, val):
        self.post_data[key] = val

    def get_post_data(self, key):
        if(key in self.post_data.keys()):
            return self.post_data[key]

        return None

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
        self.ref_data = None

    def add_view(self, mapping_dict:dict):
        for g in self.values():
            for t in g.values():
                for mname, m in mapping_dict.items():                
                    t.add_view(TestDataView(f"{mname}", m))
                    
    def set_ref_data(self, ref_data):
        self.ref_data = ref_data

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

    def dicts2table(self, dicts:Dict[str, dict], titles:list, titles_src:list, key_name:str = KEY_NAME, value_name:str = None) -> dict:
        ''' 
        use lists of dicts to form a 2D table
        '''

        table = {}
        table['Key'] = []
        for i in range(0, len(titles)):
            table['Key'].append(titles[i])

        for name_dict, dict_ in dicts.items():
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

        cfgs = {}
        for (tname, test) in self[gname].items():
            cfg = test.cfg.to_dict()
            cfgs[tname] = cfg

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

        results = {}
        for (tname, test) in self[gname].items():
            results[tname] = test.result.to_dict()

        return self.dicts2table(results, ['Name', 'Unit', 'Description'], ['text', 'unit', self.KEY_NA], value_name='val')   

    def get_views(self, gname:str) -> dict:
        '''
            get all views within a group as dicts of a 2D Table
        '''
        views = OrderedDict()
        if not gname in self.keys():
            raise KeyError(f'no test result with group name={gname}')   
        
        views_dict:dict[str, dict] = OrderedDict() # dict to collect all view result

        for tname, test in self[gname].items():
            for v_name, view in test.views.items():

                view_result = view.maps(test.result)              
                
                if not v_name in views_dict.keys():
                    views_dict[v_name] = {}
                
                views_dict[v_name][tname] = view_result

        for v_name, v_dict in views_dict.items():        
            views[v_name] = self.dicts2table(v_dict, ['Name', 'Unit', 'Description'], ['text', 'unit', self.KEY_NA], value_name='val')   

        return views

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
    def gen_test_item(cls, cfg_offset_dict:dict, keys:list, idx:int, trace:list, result:dict, gname_templ:str, cfg_name:str):
        
        key = keys[idx]

        node = cfg_offset_dict[key]

        for val in node:
            trace_cp = deepcopy(trace) # local copy
            trace_cp.append((key, val))

            cfg_name_cp = copy(cfg_name) 
            if idx > 0: # combine key other than 1st level into cfg_name
                cfg_name_cp += f'{key}={val},'

            if idx == len(keys) - 1: # reach the deepest level, create a cfg_offset instance

                (k0, v0) = trace[0] 

                gname = gname_templ
                try:
                    gname = gname_templ % v0 # use first level's key as group name
                except:
                    pass 

                vals = []

                for (k,v) in trace_cp: 
                    vals.append(v)

                if not gname in result.keys():
                    result[gname] = {}
                    result[gname]['keys'] = keys

                result[gname][cfg_name_cp[0:-1]] = vals               

            else:
                # iterate for deeper level
                cls.gen_test_item(cfg_offset_dict, keys, idx + 1, trace_cp, result, gname_templ, cfg_name_cp)

    @classmethod
    def gen_cfg_offset_table(cls, cfg_offset_dict:OrderedDict, gname_templ:str) -> dict:
        '''
        use a cfg offset dict to create a cfg offset 2D table, which will be used to create all
        configs in a data set 
        '''
        trace=[]
        cfg_offset_table={}    
        cls.gen_test_item(cfg_offset_dict, list(cfg_offset_dict.keys()), 0, trace, cfg_offset_table, gname_templ, "")

        return cfg_offset_table

    @classmethod
    def gen_test_dataset(cls, cfg_ref:dict, cfg_offset_table:dict, ds_name):
        '''
        use cfg_offset_table create all the combinations of configuations
        '''
        result = OrderedDict()

        for gname, gnode in cfg_offset_table.items():

            param_keys = gnode['keys']
            # for each row in the dict (except for key row)
            # use recursion to create a variation
            for (test_name, test_cfg) in gnode.items(): 
                if test_name == 'keys':
                    continue
                else:
                    cfg = deepcopy(cfg_ref)

                    for i in range(0, len(param_keys)):
                        k = param_keys[i]
                        # if k in cfg.keys():
                        cfg[k] = test_cfg[i]

                    item = {}
                    item["cfg"] = cfg
                    item['result'] = { }

                    if not gname in result.keys():
                        result[gname]= {}

                    root = result[gname]  
                    root[test_name] = item  

        json_cfg_batch = json.dumps(result,indent=2)
        ds_test = TestDataSet.from_json(ds_name, json_cfg_batch)  

        return ds_test    



