from collections import OrderedDict
from copy import deepcopy
import datetime
import shutil
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

from model import TestDataSet, TestConfig, TestDataView, TestResult,Variable, TestConstants
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
            cfg_set = OrderedDict()
            # head of table config
            cfg_set['Key'] = ['Name', '', 'Description']

            result_set = OrderedDict()
            result_set['Key'] = ['Name', 'Unit', 'Description']

            view_set = OrderedDict()            
            for group in ds_test:       
                for test in group:
                    cfg = test.cfg.to_dict()
                    for (k, v) in cfg.items():
                        if k in cfg_set.keys():
                            cfg_set[k].append(v)
                        else:
                            # for first column, add key and two '-' to align with results
                            cfg_set[k] = [k, '-', '-', v]

                    result = test.result.to_dict()

                    for (k, v) in result.items():
                        if k in result_set.keys():
                            result_set[k].append(v.val)
                        else:
                            # for first column, add text, unit
                            result_set[k] = [v.text, v.unit, '_', v.val]
                    # wbk.close()       

                    if test.has_view():
                        for view in test.views.values():                          
                            map_result = view.maps(test.result)
                            if view.name in view_set.keys():
                                data_table = view_set[view.name]
                                for (k, v) in map_result.items():
                                    data_table[k].append(v.val)
                            else:                  
                                data_table = {}
                                data_table['Key'] = ['Name', 'Unit', 'Description']
                                for (k, v) in map_result.items():
                                    data_table[k] = [v.text, v.unit, '_', v.val]
                                view_set[view.name] = data_table


            from pathlib import Path

            ex_file=Path(filename)       
            if ex_file.exists():
                shutil.move(filename, filename+".bak")

            wbk = xw_app.books.add()      
            sht = wbk.sheets.add() # 30 - max lenth for a sheet name
            ex: ExcelHelper = ExcelHelper(sht)
            u, l = 1, TestConstants.DATA_FILE_COL_START
            title = {"Table 1": "Test Config"}

            # config
            (b, r) = ex.write_table(cfg_set, title=title, up=u, left=l, linespacing=True)

            title = {"Table 2": "Results"}
            (b, r) = ex.write_table(result_set, title=title, up=b, left=l, linespacing=True)
            
            if view_set: # not empty
                idx = 3
                for (view_name, data_table) in view_set.items():                    
                    title = {f"Table {idx}": view_name}
                    (b, r) = ex.write_table(data_table, title=title, up=b, left=l, linespacing=True)
                    idx += 1

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

    test_mode = False # =True: use mock data to accelerate development
    use_PCHE = True # =True: use PCHE as recuperator

    model_name = "Steps.Cycle.TP_RCBCycleMock" if test_mode else "Steps.Cycle.TP_RCBCycle_PCHE" if use_PCHE else "Steps.Cycle.TP_RCBCycle"

    mdot_main_des = 125

    # referred base cfg
    cfg_ref = {
        "mdot_main":mdot_main_des,
        "mdot_heater_hot": 55,
        "T_heater_hot" : from_degC(800),
        "T_cooler_cold" : from_degC(35),
        "gamma": 0.45
    }

    # cfg with varied parameters from the base cfg
    cfg_offset = {} 
    mapping = []
    if test_mode:   
        cfg_offset["mdot_main"] = list(map(lambda x: x * mdot_main_des/100, [75, 100, 120]))
        cfg_offset["T_heater_hot"] = list(map(lambda x: from_degC(x), [550, 700]))
        mapping = [
        {
            "eta_pb":"eta_cycle",
            "UA_HTR" : "UA_HTR",
            "UA_LTR" : "UA_LTR",
            "UA_cooler" : "UA_cooler",
            "UA_heater" : "UA_heater",
            "W_MC": "W_comp",
            "W_RC": "W_recomp",
            "W_turb": "W_turb",
            "Q_HTR": "Q_HTR",
            "Q_LTR": "Q_LTR",
            "Q_cooler": "Q_cooler",
            "Q_heater": "Q_heater",
            "ex_HTR":"ex_HTR",
            "ex_LTR":"ex_LTR",
            "ex_comp":"ex_comp",
            "ex_recomp":"ex_recomp",
            "ex_turbine":"ex_turbine",
            "ex_cooler":"ex_cooler",
            "ex_heater":"ex_heater",
        },
        {
            "rc1.T":"T_amb",
            "r05.T": "TIT"
        }]
    else:
        # reduced size batch
        # cfg_offset["mdot_main"] = list(map(lambda x: x * mdot_main_des/100, [75, 100]))        
        # cfg_offset["T_heater_hot"] = list(map(lambda x: from_degC(x), [550, 700]))
        # cfg_offset["T_cooler_cold"] = list(map(lambda x: from_degC(x), [30]))
        # full size batch
        # load ratio < 0.75 leads error, use following values instead 
        # cfg_offset["mdot_main"] = list(map(lambda x: x * mdot_main_des/100, [50, 75, 100, 120]))  
        cfg_offset["mdot_main"] = list(map(lambda x: x * mdot_main_des/100, [75, 90, 100, 120]))       
        cfg_offset["T_heater_hot"] = list(map(lambda x: from_degC(x), [550, 600, 650, 700]))
        cfg_offset["T_cooler_cold"] = list(map(lambda x: from_degC(x), [30, 35, 40, 45]))
        # cfg_offset["mdot_heater_hot"] = 55
        # cfg_offset["gamma"] =[0.3, 0.325, 0.35, 0.375, 0.4, 0.45]	
        cfg_offset["gamma"] =[0.3, 0.35, 0.4, 0.45]	
        # src -> dst
        mapping = [
        {
            "eta_pb":"eta_cycle",
            "eta_turb":"eta_turb",
            "eta_MC":"eta_MC",
            "eta_RC":"eta_RC",            
            "UA_HTR" : "UA_HTR",
            "UA_LTR" : "UA_LTR",
            "UA_cooler" : "UA_cooler",
            "UA_heater" : "UA_heater",
            "W_MC": "W_comp",
            "W_RC": "W_recomp",
            "W_turb": "W_turb",
            "W_net": "W_net",
            "Q_HTR": "Q_HTR",
            "Q_LTR": "Q_LTR",
            "Q_cooler": "Q_cooler",
            "Q_heater": "Q_heater",
            "ex_HTR":"ex_HTR",
            "ex_LTR":"ex_LTR",
            "ex_comp":"ex_comp",
            "ex_recomp":"ex_recomp",
            "ex_turbine":"ex_turbine",
            "ex_cooler":"ex_cooler",
            "ex_heater":"ex_heater",
        }]
    
    ds_name = f'10MW off-Design{"_PCHE" if use_PCHE else ""}' + ' sim {:%Y-%m-%d-%H-%M-%S}'.format(datetime.datetime.now())
    ds_test = gen_batch_cfg(cfg_ref, cfg_offset, gname_template='10MW off-Design(mdot=%s)',ds_name=ds_name)
    
    
    for g in ds_test:
        for t in g:
            idx = 1
            for m in mapping:                
                t.add_view(TestDataView(f"performance map {idx}", m))
                idx += 1

    json_str = ds_test.to_json()    
    print(json_str)

    ports = [
        'r01', 'r02', 'r03', 'r04', 'r05',
        'r06', 'r07', 'r08', 'r08_source', 'r08a', 'r08b', 'r09', 'r10',
        'rh1', 'rh2', 'rc1', 'rc2']

    # (propName, unit)
    props = [
        # ('fluid', '1'), # error in parsing this props, since returned name(string) are not digitial values
        ('T', 'K'), 
        ('p', 'Pa'),
        ('rho', 'kg/m3'),
        ('w', 'kg/s'),
        ('h', 'kJ/kg'),
        ('s', 'kJ/kgK')]

    # props = [
    #     ('T','K'), 
    #     ('p', 'Pa'),
    #     ('h', 'kJ/kg'),
    #     ('w', 'kg/s')]

    sol_dict = OrderedDict()

    for p in ports:
        for (prop, unit) in props:
            key = f'{p}.{prop}'
            text = f'{prop}_{p[1:]}'            
            sol_dict[key] = Variable(key, unit, text=text)    

    # specilized solutions
    var_sp = [
        Variable('W_net', 'MW'),
        Variable('Q_heater', 'MW'),
        Variable('eta_pb', '%'),
        Variable('SR', '%'),
        Variable('W_turb', 'MW'),
        Variable('eta_turb', '%'),
        Variable('W_MC', 'MW'),
        Variable('eta_MC', '%'),
        Variable('W_RC', 'MW'),
        Variable('eta_RC', '%'),
        Variable('Q_HTR', 'MW'),
        Variable('dT1_HTR', 'K'),
        Variable('dT2_HTR', 'K'),
        Variable('T_ltmd_HTR', 'K'),
        Variable('UA_HTR', 'MW/(S^2 K)'),
        Variable('Q_LTR', 'MW'),
        Variable('dT1_LTR', 'K'),
        Variable('dT2_LTR', 'K'),
        Variable('T_ltmd_LTR', 'K'),
        Variable('UA_LTR', 'MW/(S^2 K)'),
        Variable('T_heater_hot_out', 'K')        
    ]

    for var in var_sp:
        sol_dict[var.key] = var

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
            'stopTime=100',
            'stepSize=10',
            'solver=dassl',
            '-nls=homotopy',
            '-lv=LOG_DEBUG,LOG_INIT,LOG_NLS,LOG_NLS_V,LOG_STATS'],
        solution_dict=sol_dict, 
        ds_test=ds_test)   

    # exp.save_results(ds_test)

    print('All done!')

###
if __name__ == "__main__":
    main()