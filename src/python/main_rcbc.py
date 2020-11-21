from copy import deepcopy
from logging import exception
from math import trunc
from os.path import join, sep
from OMPython import OMCSessionZMQ
from OMPython import ModelicaSystem
from datetime import datetime as dt
from numpy.core import test
from numpy.core.fromnumeric import size
from numpy.lib.function_base import append
import pandas as pd
import xlwings as xw

import csv
import os
from os import path,sep
import json
import numpy as np
# from ex.CoolPropQuery import result

from model import TestDataSet, TestConfig, TestResult, TestConstants
from plotlib import PlotManager, DataSeries, AxisType
from physics import Temperature, Pressure, MDot
from experiments import Experiments
from utils import ExcelHelper

class RCBC_simulation(Experiments):
   pass    



def from_degC(T_c) -> float :
    return T_c + 273.15

def gen_cfg(cfg_ref:dict, cfg_offset:dict, keys:list, idx:int, trace:list, result:dict, gname_template:str):
    '''
    use iteration method to walk through all the combinations of configuations
    '''
    key = keys[idx]

    for val in cfg_offset[key]:
        trace_cp = deepcopy(trace) # local copy
        trace_cp.append((key, val))

        if idx == len(keys) - 1:
            # reach the deepest level, create cfg instance
            cfg = deepcopy(cfg_ref)   

            (k0, v0) = trace_cp[0] # use first level's key as group name
            gname = gname_template % v0
            if not gname in result.keys():
                result[gname]= {}

            root = result[gname]

            cfg_name = ""
            count = 0
            for (k,v) in trace_cp: # use variables above level 1 to create item name
                cfg[k] = v
                if k in keys and count > 0:
                    cfg_name += f'{k}={v},'
                count+=1
            item = {}
            item["cfg"] = cfg
            item['result'] = { }
            root[cfg_name[0:-1]] = item            
            # result.append(cfg)
        else:
            # iterate for deeper level
            gen_cfg(cfg_ref, cfg_offset, keys, idx + 1, trace_cp, result, gname_template)

def gen_batch_cfg(cfg_ref:dict, cfg_offset:dict, gname_template:str) -> TestDataSet:
    '''
    generate batch cfgs for parameter sweep based on a referenced base config    
    '''

    keys = []
    for k in cfg_offset.keys():
        v = cfg_offset[k]
        if isinstance(v, list):
            keys.append(k)

    # keys = ["mdot_main", "T_heater_hot", "T_cooler_cold"]
    trace=[]
    batch_cfg={}    
    gen_cfg(cfg_ref, cfg_offset, keys, 0, trace, batch_cfg, gname_template)

    json_cfg_batch = json.dumps(batch_cfg,indent=2)

    ds_test = TestDataSet.from_json("demoTest2", json_cfg_batch)
    # print(json_cfg_batch)
    return ds_test

def json_IO_test():
    mdot_main_des = 125 

    # json configuration for this experiment  
    # I preferred this compact way of formatting, not the key-value one.  
    test_data_set_demo1 = {    
        "10MW off-Design(50%)": 
        { # group 1
            "T_H=550,T_C=10": 
            { # test 1
                "cfg":{
                    "mdot_main":mdot_main_des * 0.5, # 62.5
                    "mdot_heater_hot":55,
                    "T_heater_hot" : from_degC(550),
                    "T_cooler_cold" : from_degC(10)
                },
                "result": {}
            },
            "T_H=550,T_C=20": { # test 2
                "cfg":{
                    "mdot_main":mdot_main_des * 0.5, # 62.5
                    "mdot_heater_hot":55,
                    "T_heater_hot" : from_degC(550),
                    "T_cooler_cold" : from_degC(20)
                },
                "result": {}

            }, 
            "T_H=600,T_C=10": { # test 3
                "cfg":{
                    "mdot_main":mdot_main_des * 0.5, # 62.5
                    "mdot_heater_hot":55,
                    "T_heater_hot" : from_degC(600),
                    "T_cooler_cold" : from_degC(10)
                },
                "result": {}
            }
        },
        "10MW off-Design(100%)": 
        {# group 2
            "T_H=550,T_C=10":
            { # test 1
                "cfg":
                {
                    "mdot_main":mdot_main_des * 1, # 125
                    "mdot_heater_hot":55,
                    "T_heater_hot" : from_degC(550),
                    "T_cooler_cold" : from_degC(10)
                },
                "result": {}
            },
            "T_H=550,T_C=20": { # test 2                    
                "cfg":
                {
                    "mdot_main":mdot_main_des * 1, # 125
                    "mdot_heater_hot":55,
                    "T_heater_hot" : from_degC(550),
                    "T_cooler_cold" : from_degC(20)
                },
                "result": {}
            }, 
            "T_H=600,T_C=10": { # test 3                    
                "cfg":{
                    "mdot_main":mdot_main_des * 1, # 125
                    "mdot_heater_hot":55,
                    "T_heater_hot" : from_degC(600),
                    "T_cooler_cold" : from_degC(10)
                },
                "result": {}
            }       
        }              
    }
    ds_test = TestDataSet.from_json("demoTest", json.dumps(test_data_set_demo1))



def main(work_root = []):
    # root path of modelica root
    from pathlib import Path

    if work_root == []:
        work_root = os.path.abspath(os.curdir)   

    mdot_main_des = 125

    # referred base cfg
    cfg_ref = {
        "mdot_main":mdot_main_des,
        "mdot_heater_hot": 55,
        "T_heater_hot" : from_degC(550),
        "T_cooler_cold" : from_degC(10)
    }

    # cfg with varied parameters from the base cfg
    cfg_offset = {}    
    cfg_offset["mdot_main"] = list(map(lambda x: x * mdot_main_des/100, [50, 75, 100, 120]))
    cfg_offset["T_heater_hot"] = list(map(lambda x: from_degC(x), [550, 600, 650, 700]))
    cfg_offset["T_cooler_cold"] = list(map(lambda x: from_degC(x), [10, 20, 30, 40]))
    cfg_offset["mdot_heater_hot"] = 55

    ds_test = gen_batch_cfg(cfg_ref, cfg_offset, gname_template='10MW off-Design(%s%%)')
   
    json_str = ds_test.to_json()
    print(json_str)

    ports = ['r01', 'r02', 'r03', 'r04', 'r05', 'r05a',
    'r06', 'r07', 'r08', 'r08a', 'r08b', 'r09', 'r10',
    'rh1', 'rh2', 'rc1', 'rc2']

    props = ['T', 'p', 'h','w']

    sol_keys = []

    for p in ports:
        for prop in props:
            sol_keys.append(f'{p}.{prop}')

    test = RCBC_simulation(
    work_root, 
    "Steps.Cycle.TP_RCBCycleMock", 
    ex_dlls=[
        'D:/Workspace/Steps/lib/MyProps.dll', 
        'D:/Workspace/Steps/lib/libexternalmedialib.dll'],
    modelica_libs=[
        'D:/Workspace/ExternalMedia/Modelica/ExternalMedia 3.2.1/package.mo', 
        'D:/Workspace/ThermoPower/ThermoPower/package.mo',
        'D:/Workspace/SolarTherm/SolarTherm/package.mo',
        "Modelica 3.2.1"]
    )

    test.simulate(sim_ops=[
            'startTime=0', 
            'stopTime=1',
            'stepSize=1',
            'solver=dassl',
            '-nls=homotopy',
            '-lv=LOG_DEBUG,LOG_INIT,LOG_NLS,LOG_NLS_V,LOG_STATS'],
    solution_names=sol_keys, ds_test=ds_test)

    test_result = ds_test.result

    print('All done!')


###
if __name__ == "__main__":
    main()