from copy import deepcopy
import datetime
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

   def save_results(self, ds_test : TestDataSet):
        '''
        save the simulation result into files   
        override save_results() to save all result in one excel file         
        '''
        # create directory for current test batch
        filename = path.join(self.path_out, ds_test.name + '.xlsx')

        xw_app: xw.App = xw.App(visible=False)       

        try:
            cfgs = {}
            results = {}
            for group in ds_test:                
                gname = group.name
                for test in group:
                    cfg = test.cfg.to_dict()
                    for (k, v) in cfg.items():
                        if k in cfgs.keys():
                            cfgs[k].append(v)
                        else:
                            cfgs[k] = [v]

                    result = test.result.to_dict()

                    for (k, v) in result.items():
                        if k in results.keys():
                            results[k].append(v)
                        else:
                            results[k] = [v]
                    
                    # wbk.close()

            wbk = xw_app.books.add()      
            sht = wbk.sheets.add() # 30 - max lenth for a sheet name
            ex: ExcelHelper = ExcelHelper(sht)
            u, l = 1, TestConstants.DATA_FILE_COL_START
            title = {"Table 1": "Test Config"}

            # config
            (b, r) = ex.write_table(cfgs, title=title, up=u, left=l, linespacing=True)

            title = {"Table 2": "Results"}
            (b, r) = ex.write_table(results, title=title, up=b, left=l, linespacing=True)

            wbk.save(filename)                    

        finally:
            xw_app.quit() # close the xlwings app finally  


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

def gen_batch_cfg(cfg_ref:dict, cfg_offset:dict, gname_template:str, ds_name = "demoTest") -> TestDataSet:
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

    ds_test = TestDataSet.from_json(ds_name, json_cfg_batch)
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
    if work_root == []:
        work_root = os.path.abspath(os.curdir)  

    test_mode = True # =True: use mock data to accelerate development
    model_name = "Steps.Cycle.TP_RCBCycleMock" if test_mode else "Steps.Cycle.TP_RCBCycle"

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

    if test_mode:   
        cfg_offset["mdot_main"] = list(map(lambda x: x * mdot_main_des/100, [75, 100]))
        cfg_offset["T_heater_hot"] = list(map(lambda x: from_degC(x), [550, 700]))
    else:
        cfg_offset["mdot_main"] = list(map(lambda x: x * mdot_main_des/100, [50, 75, 100, 120]))        
        cfg_offset["T_heater_hot"] = list(map(lambda x: from_degC(x), [550, 600, 650, 700]))
        cfg_offset["T_cooler_cold"] = list(map(lambda x: from_degC(x), [10, 20, 30, 40]))
        # cfg_offset["mdot_heater_hot"] = 55
    ds_name = '10MW off-Desin sim {:%Y-%m-%d-%H-%M-%S}'.format(datetime.datetime.now())
    ds_test = gen_batch_cfg(cfg_ref, cfg_offset, gname_template='10MW off-Design(%s%%)',ds_name=ds_name)
   
    json_str = ds_test.to_json()
    print(json_str)

    ports = [
        'r01', 'r02', 'r03', 'r04', 'r05', 'r05a',
        'r06', 'r07', 'r08', 'r08a', 'r08b', 'r09', 'r10',
        'rh1', 'rh2', 'rc1', 'rc2']

    props = ['T', 'p', 'h','w']

    sol_keys = []

    for p in ports:
        for prop in props:
            sol_keys.append(f'{p}.{prop}')

    exp:Experiments = RCBC_simulation(
        work_root, 
        model_name=model_name, 
        ex_dlls=[
            'D:/Workspace/Steps/lib/MyProps.dll', 
            'D:/Workspace/Steps/lib/libexternalmedialib.dll'],
        modelica_libs=[
            'D:/Workspace/ExternalMedia/Modelica/ExternalMedia 3.2.1/package.mo', 
            'D:/Workspace/ThermoPower/ThermoPower/package.mo',
            'D:/Workspace/SolarTherm/SolarTherm/package.mo',
            "Modelica 3.2.1"])

    exp.simulate(
        sim_ops=[
            'startTime=0', 
            'stopTime=1',
            'stepSize=1',
            'solver=dassl',
            'nls=homotopy',
            '-lv=LOG_DEBUG,LOG_INIT,LOG_NLS,LOG_NLS_V,LOG_STATS'],
        solution_names=sol_keys, 
        ds_test=ds_test)   

    exp.save_results(ds_test)

    print('All done!')

###
if __name__ == "__main__":
    main()