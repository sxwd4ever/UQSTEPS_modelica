from collections import OrderedDict
from copy import deepcopy
import datetime
import shutil
from logging import exception
from math import trunc
from os.path import join, sep
from typing import Tuple
from OMPython import OMCSessionZMQ
from OMPython import ModelicaSystem
from datetime import datetime as dt
from numpy.core import test
from numpy.core.fromnumeric import size
from numpy.lib.function_base import append
from enum import Enum

import os
from os import path, sep

import numpy as np
from model import TestDataSet,TestConfig, TestItem,Variable, TestConstants
from plotlib import PlotManager, DataSeries, AxisType
from physics import Temperature, Pressure, MDot
from experiments import Experiment, PCHEExperiment
from utils import from_degC, from_bar

class MarchionniTest(PCHEExperiment):
    '''
    script for test against Marchionni [2019]
    '''
    def gen_sol_dict(self, cfg:dict) -> dict:
        sol_dict = super().gen_sol_dict(cfg)

                # specilized solutions
        var_sp = [               
            Variable('T_hot_out_act', 'K'),
            Variable('T_cold_out_act', 'K'),
            Variable('dp_hot_act', 'pa'),
            Variable('dp_cold_act', 'pa')
        ]

        for v in var_sp:
            sol_dict[v.key] = v

        return sol_dict
        
    def post_process(self, test: TestItem):
        return super().post_process(test)

    def gen_plot_cfgs(self, values, meta_cfg):
        '''
            use template to generate all plot configs 
        '''
        return {
            'T_comparison':{
                'series':[{
                    'name' : 'T_hot',
                    'cs': 'r-o',
                    'x':
                    {
                        'value': values['x_hot'] * 1e3,
                        'axis_range': meta_cfg['axis_x'],
                        'label': meta_cfg['label_x']
                    },
                    'y':
                    {
                        'value': values['T_hot'] - 273.15,
                        'axis_range': meta_cfg['axis_T'],
                        'label': meta_cfg['label_T']
                    }
                },
                {
                    'name' : 'T_cold',
                    'cs': 'b-v',
                    'x':
                    {
                        'value': values['x_cold'] * 1e3,
                        'axis_range': meta_cfg['axis_x'],
                        'label': meta_cfg['label_x']
                    },
                    'y':
                    {
                        'value': values['T_cold'] - 273.15,
                        'axis_range': meta_cfg['axis_T'],
                        'label': meta_cfg['label_T']
                    }
                }],
                "imgfile" : meta_cfg["imgfile_T"]
            },
            'dp_comparison':{
                'series':[{
                    'name' : 'dp_hot',
                    'cs': 'r-o',
                    'x':
                    {
                        'value': values['x_hot'] * 1e3,
                        'axis_range': meta_cfg['axis_x'],
                        'label': meta_cfg['label_x']
                    },
                    'y':
                    {
                        'value': values['dp_hot'],
                        'axis_range': meta_cfg['axis_dp'],
                        'label': meta_cfg['label_dp']
                    }
                },
                {
                    'name' : 'dp_cold',
                    'cs': 'b-v',
                    'x':
                    {
                        'value': values['x_cold'] * 1e3,
                        'axis_range': meta_cfg['axis_x'],
                        'label': meta_cfg['label_x']
                    },
                    'y':
                    {
                        'value': values['dp_cold'],
                        'axis_range': meta_cfg['axis_dp'],
                        'label': meta_cfg['label_dp']
                    }
                }],
                "imgfile" : meta_cfg["imgfile_dp"]
            },  
            'gamma_comparison':{
                'series':[{
                    'name' : 'gamma_hot',
                    'cs': 'r-o',
                    'x':
                    {
                        'value': values['x_hot'] * 1e3,
                        'axis_range': meta_cfg['axis_x'],
                        'label': meta_cfg['label_x']
                    },
                    'y':
                    {
                        'value': values['hc_hot'] / 1e3,
                        'axis_range': meta_cfg['axis_gamma'],
                        'label': meta_cfg['label_gamma']
                    }
                },
                {
                    'name' : 'gamma_cold',
                    'cs': 'b-v',
                    'x':
                    {
                        'value': values['x_cold'] * 1e3,
                        'axis_range': meta_cfg['axis_x'],
                        'label': meta_cfg['label_x']
                    },
                    'y':
                    {
                        'value': values['hc_cold'] / 1e3,
                        'axis_range': meta_cfg['axis_gamma'],
                        'label': meta_cfg['label_gamma']
                    }
                }],
                "imgfile" : meta_cfg["imgfile_gamma"]
            }                  
        } 

class ExpType(Enum):
    LOAD_PRE_EXP = 1,
    VS_CFD = 2,  # "aginst marchionni's 1D_vs_3D"
    FULL_SCALE = 3, # full scale off design in Sec. 4 of [Marchionni 2019]

def main(work_root = []):
    # root path of modelica root
    if work_root == []:
        work_root = os.path.abspath(os.curdir)  

    # parameters initialization for this simulation
    exp_type:ExpType = ExpType.VS_CFD
    model_name = "Steps.Test.TestTP_PCHE_Marchionni"  
    exp_name = 'Test-Marchionni {:%Y-%m-%d-%H-%M-%S}'.format(datetime.datetime.now())            
    # end of parameter initialization

    exp:Experiment = MarchionniTest(
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

    plot_metacfg_base = {
        "axis_x": [0, 272],
        "axis_T" : [80, 440], # degC  
        "axis_dp" : [0, 10], # kPa
        "axis_gamma" : [1.7, 3.5], # kW/(m^2K)
        "label_x" : 'Z (m)',
        "label_T" : 'T (K)',
        "label_dp" : '$\\Delta~p~(kPa)$ (K)',  
        "label_gamma" : '$Local heat transfer co (kW/(m^2 K)$',
        "imgfile_T": "Marchionni_Fig_04_T_1D_vs_3D.png",
        "imgfile_dp": "Marchionni_Fig_05_dp_1D_vs_3D.png",  
        "imgfile_gamma": "Marchionni_Fig_06_gamma_1D_vs_3D.png"                  
    }

    plot_metacfg_offest = {} 
    
    if exp_type == ExpType.VS_CFD:

        # referred base cfg
        cfg_ref = {
            # geometry parameters
            "N_ch": 1e4, # "channel number"
            "N_seg": 20, # "segments number"
            "D_ch": 2e-3, # "channel diameter, semi circular tube"
            "L_fp": 272e-3, # "channel flow path length"
            "L_pitch": 12.3e-3, # "pitch length"
            "a_phi": 36.0, # "pitch angle, degree"
            "H_ch": 3.26e-3, # "Height of the solid domain, containing one cold tube and one hot tube"
            "W_ch": 1.27e-3 * 2, # "Width of the solid domain"
            # boundary conditon
            "G_hot_in": 509.3, # "hot inlet velocity m/s";
            "G_cold_in": 509.3, # "cold inlet velocity m/s";
            "p_hot_in": from_bar(75), # "hot inlet pressure";
            "p_cold_in":from_bar(150), # "cold inlet pressure";
            "T_hot_in": from_degC(400), # "hot inlet temperature, K";
            "T_hot_out": from_degC(140), # "cold outlet temperature, K";
            "T_cold_in": from_degC(100), # "cold inlet temperature, K";
            "T_cold_out": from_degC(300), # "cold outlet temperature, K";
            "kc_dp": 1 # "pressure drop correction coefficient"
        } 

        # cfg with varied parameters from the base cfg
        cfg_offset = OrderedDict() 
        cfg_offset['Group 1'] = {
                # values in Table 3 in Meshram [2016]
                # Only with kc_dp = 1 at High T and kc_dp = 2 at Low T can I get good agreement with 
                # Meshram's DP result. 
                "keys" : ["a_phi", "kc_dp"],   
                "a_phi=5" : [5.0, 1],
                "a_phi=10" : [10.0, 1],
                "a_phi=20" : [20.0, 1],                     
                "a_phi=30" : [30.0, 1],
                "a_phi=35" : [35.0, 1],
                "a_phi=40" : [40.0, 1],
                "a_phi=45" : [45.0, 1]
            }          

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

        ds_exp = TestDataSet.gen_test_dataset(cfg_ref, cfg_offset, exp_name)

        ds_exp.add_view(mapping, view_name_tmpl="performance map")

        json_str = ds_exp.to_json()    
        print(json_str)    

        exp.simulate(
            sim_ops=[
                'startTime=0', 
                'stopTime=10',
                'stepSize=2',
                'solver=dassl',
                '-nls=homotopy',
                '-lv=LOG_DEBUG,LOG_INIT,LOG_NLS,LOG_NLS_V,LOG_STATS'],
            solution_dict=exp.gen_sol_dict(cfg_ref), 
            ds_test=ds_exp,
            append_save=False)   

        exp.save_results(ds_exp)
        exp.plot_results(ds_exp, plot_metacfg_base, plot_metacfg_offest)        

    elif exp_type == ExpType.FULL_SCALE:
        # referred base cfg
        cfg_ref = {
            # geometry parameters
            "N_ch": 54 * 42, # "channel number"
            "N_seg": 20, # "segments number"
            "D_ch": 2e-3, # "channel diameter, semi circular tube"
            "L_fp": 1012e-3, # "channel flow path length"
            "L_pitch": 12.3e-3, # "pitch length"
            "a_phi": 45.0, # "pitch angle, degree"
            "H_ch": 3.26e-3, # "Height of the solid domain, containing one cold tube and one hot tube"
            "W_ch": 1.27e-3 * 2, # "Width of the solid domain"
            # boundary conditon
            "mdot_hot_in": 2.06,  # "hot inlet mass flow rate kg/s";
            "mdot_cold_in": 2.06,  # "cold inlet mass flow rate kg/s";
            "p_hot_in": from_bar(75), # "hot inlet pressure";
            "p_cold_in":from_bar(125), # "cold inlet pressure";
            "T_hot_in": from_degC(344.3), # "hot inlet temperature, K";
            "T_hot_out": from_degC(81), # "cold outlet temperature, K";
            "T_cold_in": from_degC(72.9), # "cold inlet temperature, K";
            "T_cold_out": from_degC(283), # "cold outlet temperature, K";
            "kc_dp": 1 # "pressure drop correction coefficient"
        } 

        # cfg with varied parameters from the base cfg
        cfg_offset = {} 
        cfg_offset['phi=5'] = {
                # values in Table 3 in Meshram [2016]
                # Only with kc_dp = 1 at High T and kc_dp = 2 at Low T can I get good agreement with 
                # Meshram's DP result. 
                "keys" : ["mdot_hot_in", "mdot_cold_in", "T_cold_in", "a_phi"],   
                "Design point" : [2.06, 2.06, from_degC(72.9), 5],
                "off-design #1" : [1.57, 1.57, from_degC(72.9), 5],
                "off-design #2" : [2.09, 2.09, from_degC(87.5), 5],                  
                "off-design #3" : [2.09, 2.09,  from_degC(62), 5],
                "off-design #4" : [2.62, 2.62, from_degC(72.9), 5]
            }  
        cfg_offset['phi=10'] = {
                # values in Table 3 in Meshram [2016]
                # Only with kc_dp = 1 at High T and kc_dp = 2 at Low T can I get good agreement with 
                # Meshram's DP result. 
                "keys" : ["mdot_hot_in", "mdot_cold_in", "T_cold_in"],   
                "Design point" : [2.06, 2.06, from_degC(72.9), 10],
                "off-design #1" : [1.57, 1.57, from_degC(72.9), 10],
                "off-design #2" : [2.09, 2.09, from_degC(87.5), 10],                  
                "off-design #3" : [2.09, 2.09,  from_degC(62), 10],
                "off-design #4" : [2.62, 2.62, from_degC(72.9), 10]
            }  

        cfg_offset['phi=15'] = {
                # values in Table 3 in Meshram [2016]
                # Only with kc_dp = 1 at High T and kc_dp = 2 at Low T can I get good agreement with 
                # Meshram's DP result. 
                "keys" : ["mdot_hot_in", "mdot_cold_in", "T_cold_in", "a_phi"],   
                "Design point" : [2.06, 2.06, from_degC(72.9), 15],
                "off-design #1" : [1.57, 1.57, from_degC(72.9), 15],
                "off-design #2" : [2.09, 2.09, from_degC(87.5), 15],                  
                "off-design #3" : [2.09, 2.09,  from_degC(62), 15],
                "off-design #4" : [2.62, 2.62, from_degC(72.9), 15]
            }                       

        ds_exp = TestDataSet.gen_test_dataset(cfg_ref, cfg_offset, ds_name=exp_name)        

        json_str = ds_exp.to_json()    
        print(json_str)    

        exp.simulate(
            sim_ops=[
                'startTime=0', 
                'stopTime=10',
                'stepSize=2',
                'solver=dassl',
                '-nls=homotopy',
                '-lv=LOG_DEBUG,LOG_INIT,LOG_NLS,LOG_NLS_V,LOG_STATS'],
            solution_dict=exp.gen_sol_dict(cfg_ref), 
            ds_test=ds_exp,
            append_save=False)   

        exp.save_results(ds_exp)                   

    elif exp_type == ExpType.LOAD_PRE_EXP:
        # should search the exp_name 
        exp_name = "Test-Marchionni 2020-12-09-15-04-22"
        ds_exp = exp.load_results(exp_name, dir_name=exp_name)

        exp.plot_results(ds_exp, plot_metacfg_base, plot_metacfg_offest)   

    print('All done!')

###
if __name__ == "__main__":
    main()