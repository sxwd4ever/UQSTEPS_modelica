from collections import OrderedDict
from copy import copy,deepcopy
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
from experiments import Experiment
from utils import ExcelHelper, from_degC

class RCBC_simulation(Experiment):
    pass

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
                    "mdot_main"      : mdot_main_des * 0.5, # 62.5
                    "mdot_heater_hot": 55,
                    "T_heater_hot"   : from_degC(550),
                    "T_cooler_cold"  : from_degC(10)
                },
                "result": {}
            },
            "T_H=550,T_C=20": { # test 2
                "cfg":{
                    "mdot_main"      : mdot_main_des * 0.5, # 62.5
                    "mdot_heater_hot": 55,
                    "T_heater_hot"   : from_degC(550),
                    "T_cooler_cold"  : from_degC(20)
                },
                "result": {}

            }, 
            "T_H=600,T_C=10": { # test 3
                "cfg":{
                    "mdot_main"      : mdot_main_des * 0.5, # 62.5
                    "mdot_heater_hot": 55,
                    "T_heater_hot"   : from_degC(600),
                    "T_cooler_cold"  : from_degC(10)
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
                    "mdot_main"      : mdot_main_des * 1, # 125
                    "mdot_heater_hot": 55,
                    "T_heater_hot"   : from_degC(550),
                    "T_cooler_cold"  : from_degC(10)
                },
                "result": {}
            },
            "T_H=550,T_C=20": { # test 2                    
                "cfg": 
                {
                    "mdot_main"      : mdot_main_des * 1, # 125
                    "mdot_heater_hot": 55,
                    "T_heater_hot"   : from_degC(550),
                    "T_cooler_cold"  : from_degC(20)
                },
                "result": {}
            }, 
            "T_H=600,T_C=10": { # test 3                    
                "cfg":{
                    "mdot_main"      : mdot_main_des * 1, # 125
                    "mdot_heater_hot": 55,
                    "T_heater_hot"   : from_degC(600),
                    "T_cooler_cold"  : from_degC(10)
                },
                "result": {}
            }       
        }              
    }
    ds_test = TestDataSet.from_json("demoTest", json.dumps(test_data_set_demo1))


def gen_performance_map(mdot_heat, mdot_main, T_HTF_in, T_amb, work_root=[], test_mode = False, use_PCHE = True):
    '''
        parameter sweep function to generate the RCBC performance map

        mdot_heat: array containing candidates of HTF (hot transfer fluid) mass flow rates, kg/s
        mdot_main: array containing candidates of main working fluid (s-CO2) mass flow rates, kg/s
        T_HTF_in : array containing candidates of HTF inlet temperatures, K
        T_amb    : array containing candidates of cooler's cooling fluid, K.

        work_root: work root for current simulation, Default = [] and current directory will be used.  
        test_mode: Flag for test mode, for debug purpose only. Default = false
        use_PCHE : if use PCHE as recuperator (HTR, LTR). Default = true
    '''

   # root path of modelica root
    if work_root == []:
        work_root = os.path.abspath(os.curdir)  

     # =True: use mock data to accelerate development
     # =True: use PCHE as recuperator

    model_name = "Steps.Cycle.TP_RCBCycleMock" if test_mode else "Steps.Cycle.TP_RCBCycle_PCHE" if use_PCHE else "Steps.Cycle.TP_RCBCycle"

    mdot_main_des = 125

    # referred base cfg
    cfg_ref = {
        "mdot_main"      : mdot_main_des,
        "mdot_heater_hot": 90,             # 55,
        "T_heater_hot"   : from_degC(800),
        "T_cooler_cold"  : from_degC(35),
        "gamma"          : 0.45
    }

    # cfg with varied parameters from the base cfg
    cfg_offset_dict = OrderedDict()
    mapping = []
    if test_mode:   
        cfg_offset_dict["mdot_main"] = list(map(lambda x: x * mdot_main_des/100, [75, 100, 120]))
        cfg_offset_dict["T_heater_hot"] = list(map(lambda x: from_degC(x), [550, 700]))
        
        # src -> dst
        mapping = {
            "performance map":{
                "eta_pb"    : "eta_cycle",
                "eta_turb"  : "eta_turb",
                "eta_MC"    : "eta_MC",
                "eta_RC"    : "eta_RC",
                "UA_HTR"    : "UA_HTR",
                "UA_LTR"    : "UA_LTR",
                "UA_cooler" : "UA_cooler",
                "UA_heater" : "UA_heater",
                "W_MC"      : "W_comp",
                "W_RC"      : "W_recomp",
                "W_turb"    : "W_turb",
                "W_net"     : "W_net",
                "Q_HTR"     : "Q_HTR",
                "Q_LTR"     : "Q_LTR",
                "Q_cooler"  : "Q_cooler",
                "Q_heater"  : "Q_heater",
                "ex_HTR"    : "ex_HTR",
                "ex_LTR"    : "ex_LTR",
                "ex_comp"   : "ex_comp",
                "ex_recomp" : "ex_recomp",
                "ex_turbine": "ex_turbine",
                "ex_cooler" : "ex_cooler",
                "ex_heater" : "ex_heater",
            },
            "demo view":{
                "r_cooler_cin.T": "T_amb",
                "r_turb_in.T"   : "TIT"
            }
        }        
    else:
        # # reduced size batch

        # cfg_offset_dict["mdot_main"] = list(map(lambda x: x * mdot_main_des/100, [75, 100]))        
        # cfg_offset_dict["T_heater_hot"] = list(map(lambda x: from_degC(x), [650, 700]))
        # cfg_offset_dict["T_cooler_cold"] = list(map(lambda x: from_degC(x), [35]))

        # full size batch
        # load ratio < 0.75 leads error, use following values instead 
        # cfg_offset["mdot_main"] = list(map(lambda x: x * mdot_main_des/100, [50, 75, 100, 120]))  
        
        cfg_offset_dict["T_heater_hot"]    = T_HTF_in
        cfg_offset_dict["mdot_heater_hot"] = mdot_heat
        cfg_offset_dict["mdot_main"]       = mdot_main
        cfg_offset_dict["T_cooler_cold"]   = T_amb
        # cfg_offset["gamma"] =[0.3, 0.325, 0.35, 0.375, 0.4, 0.45]	
        # cfg_offset_dict["gamma"] =[0.3, 0.35, 0.4, 0.45]	       

        # src -> dst
        mapping = {
            "performance map":{
                "eta_pb"    : "eta_cycle",
                "eta_turb"  : "eta_turb",
                "eta_MC"    : "eta_MC",
                "eta_RC"    : "eta_RC",
                "UA_HTR"    : "UA_HTR",
                "UA_LTR"    : "UA_LTR",
                "UA_cooler" : "UA_cooler",
                "UA_heater" : "UA_heater",
                "W_MC"      : "W_comp",
                "W_RC"      : "W_recomp",
                "W_turb"    : "W_turb",
                "W_net"     : "W_net",
                "Q_HTR"     : "Q_HTR",
                "Q_LTR"     : "Q_LTR",
                "Q_cooler"  : "Q_cooler",
                "Q_heater"  : "Q_heater",
                "ex_HTR"    : "ex_HTR",
                "ex_LTR"    : "ex_LTR",
                "ex_comp"   : "ex_comp",
                "ex_recomp" : "ex_recomp",
                "ex_turbine": "ex_turbine",
                "ex_cooler" : "ex_cooler",
                "ex_heater" : "ex_heater",
            }
        }
    
    ds_name = f'10MW off-Design{"_PCHE" if use_PCHE else ""}' + ' sim {:%Y-%m-%d-%H-%M-%S}'.format(datetime.datetime.now())
    # ds_test = TestDataSet.gen_batch_cfg(cfg_ref, cfg_offset, gname_template='10MW off-Design(mdot=%s)', ds_name=ds_name)
    cfg_table = TestDataSet.gen_cfg_offset_table(cfg_offset_dict, gname_templ='10MW off-Design(mdot=%s)')  
    
    ds_test = TestDataSet.gen_test_dataset(cfg_ref, cfg_table, ds_name)
    
    ds_test.add_view(mapping)

    json_str = ds_test.to_json()    
    print(json_str)
    
    ports = [
        'r_cooler_hout', 'r_comp_out', 'r_HTR_cin', 'r_HTR_cout', 'r_heater_cout',
        'r_HTR_hin', 'r_HTR_hout', 'r_LTR_hout', 'r_cooler_hin', 'r_recomp_in', 
        'r_recomp_out', 'r_LTR_cout', 'r_heater_hin', 'r_heater_hout', 
        'r_cooler_cin', 'r_cooler_cout']

    # (propName, unit)
    props = [
        # ('fluid', '1'), # error in parsing this props, since returned name(string) are not digital values
        ('T', 'K'), 
        ('p', 'Pa'),
        ('rho', 'kg/m3'),
        ('w', 'kg/s'),
        ('h', 'kJ/kg'),
        ('s', 'kJ/kgK')]

    sol_dict = OrderedDict()

    for p in ports:
        for (prop, unit) in props:
            key = f'{p}.{prop}'
            text = f'{prop}_{p[1:]}'            
            sol_dict[key] = Variable(key, unit, text=text)    

    # specialized solutions
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

    exp:Experiment = RCBC_simulation(
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

    exp.simulate_batch(
        sim_ops=[
            'startTime=0', 
            'stopTime=10',
            'stepSize=2',
            'solver=dassl',
            '-nls=homotopy',
            '-newtonXTol=1e-6',
            '-newtonFTol=1e-6',
            '-lv=LOG_DEBUG,LOG_INIT,LOG_NLS,LOG_STATS'],
        solution_dict=sol_dict, 
        ds_test=ds_test)   

    exp.save_results(ds_test)

    print('All done!')


if __name__ == "__main__":

    mdot_main_des = 125 # Working fluid mass flow rate design point. 

    # call the parameter sweep function to generate performance map
    gen_performance_map(
        mdot_heat = [30, 35, 40, 45],
        mdot_main = list(map(lambda x: x * mdot_main_des/100, [75, 90, 100, 120])),
        T_HTF_in  = list(map(lambda x: from_degC(x), [550, 600, 650, 700])),
        T_amb     = list(map(lambda x: from_degC(x), [30, 35, 40, 45]))
    )