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

import os
from os import path, sep

import numpy as np
from model import TestDataSet,TestConfig, TestItem,Variable, TestConstants
from plotlib import PlotManager, DataSeries, AxisType
from physics import Temperature, Pressure, MDot
from experiments import Experiment, PCHEExperiment
from utils import ExcelHelper, from_degC, from_bar,mkdir_filepath

class MeshramTest(PCHEExperiment):
    '''
    script for test against Meshram [2016]
    '''

    def gen_plot_cfgs(self, values, meta_cfg):
        '''
            use template to generate all plot configs 
        '''
        return {
            'T_comparison':{
                'series':[{
                    'name' : 'T_hot',
                    'cs': 'r-s',
                    'x':
                    {
                        'value': values['x_hot'],
                        'axis_range': meta_cfg['axis_x'],
                        'label': meta_cfg['label_x']
                    },
                    'y':
                    {
                        'value': values['T_hot'],
                        'axis_range': meta_cfg['axis_T'],
                        'label': meta_cfg['label_T']
                    }
                },
                {
                    'name' : 'T_cold',
                    'cs': 'b--s',
                    'x':
                    {
                        'value': values['x_cold'],
                        'axis_range': meta_cfg['axis_x'],
                        'label': meta_cfg['label_x']
                    },
                    'y':
                    {
                        'value': values['T_cold'],
                        'axis_range': meta_cfg['axis_T'],
                        'label': meta_cfg['label_T']
                    }
                },
                {
                    'name' : 'dp_hot',
                    'cs': 'r-^',
                    'x':
                    {
                        'value': values['x_hot'],
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
                    'cs': 'b--^',
                    'x':
                    {
                        'value': values['x_cold'],
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
                "imgfile" : meta_cfg["imgfile"]
            }                            
        }

def main(work_root = []):
    # root path of modelica root
    if work_root == []:
        work_root = os.path.abspath(os.curdir)  

    run_sim = False # True: run the simulation, False: load pre-saved exp results
    model_name = "Steps.Test.TestTP_PCHE_Meshram"  
    # referred base cfg
    cfg_ref = {
        # geometry parameters
        "N_ch": 1e4, # "channel number"
        "N_seg": 10, # "segments number"
        "D_ch": 2e-3, # "channel diameter, semi circular tube"
        "L_fp": 200e-3, # "channel flow path length"
        "L_pitch": 12e-3, # "pitch length"
        "a_phi": 36.0, # "pitch angle, degree"
        "H_ch": 3.2e-3, # "Height of the solid domain, containing one cold tube and one hot tube"
        "W_ch": 2.5e-3, # "Width of the solid domain"
        # boundary conditon
        "u_hot_in": 7.564, # "hot inlet velocity m/s";
        "u_cold_in": 1.876, # "cold inlet velocity m/s";
        "p_hot_in": from_bar(90), # "hot inlet pressure";
        "p_cold_in":from_bar(225), # "cold inlet pressure";
        "T_hot_in": 730, # "hot inlet temperature, K";
        "T_hot_out": 576.69, # "cold outlet temperature, K";
        "T_cold_in": 500, # "cold inlet temperature, K";
        "T_cold_out": 639.15, # "cold outlet temperature, K";
        "kc_dp": 1 # "pressure d"
    }

    # cfg with varied parameters from the base cfg
    cfg_offset = {} 
    cfg_offset['Group 1'] = {
            # values in Table 3 in Meshram [2016]
            # Only with kc_dp = 1 at High T and kc_dp = 2 at Low T can I get good agreement with 
            # Meshram's DP result. 
            "keys" : ["T_cold_in", "T_cold_out", "T_hot_in", "T_hot_out", "u_cold_in", "u_hot_in", "a_phi", "kc_dp"],
            "zigzag High T" : [500, 639.15, 730, 576.69, 1.876, 7.564, 36.0, 1],
            # for 'zigzag Low T', I can not retrieve values match Meshram's Fig. 5 (a) with
            # u_hot_in = 1.345, which makes mdot_hot_in = 1.614 << mdot_cold_in = 5.488. 
            # Based on the facts that for other cases, mdot_hot_in ~ mdot_cold_in, 
            # I choose mdot_hot_in = 4.501 according to u_cold_in/u_hot_in = 5.584 
            # in 'Straight Low T' the senario and get good match
            "zigzag Low T": [400, 522.23, 630, 466.69, 0.806, 4.501, 36.0, 2], 
            "Straight High T": [500, 615.48, 730, 601.83, 1.518, 6.118, 5, 1],
            "Straigth Low T": [400, 498.45, 630, 494.37, 0.842, 4.702, 5, 2]
        }
    mapping = [] 

    exp:Experiment = MeshramTest(
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

    if run_sim:
        exp_name = 'Test-Meshram {:%Y-%m-%d-%H-%M-%S}'.format(datetime.datetime.now())            
        ds_exp = TestDataSet.gen_batch_cfg(cfg_ref, cfg_offset, gname_template="Meshram_Test", ds_name=exp_name)
    
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
    else:
        # should search the exp_name 
        exp_name = "Test-Meshram 2020-12-09-11-49-05"
        ds_exp = exp.load_results(exp_name, dir_name=exp_name)

    plot_metacfg_base = {
            "axis_x" : [0, 200e-3],
            "axis_T" : [400.0, 750.0],   
            "axis_dp" : [0, 80],
            "label_x" : 'Z (m)',
            "label_T" : 'T (K)',
            "label_dp" : '$\\Delta~p~(kPa)$ (K)'
        }

    plot_metacfg_offest = {
        "zigzag High T" :
        {
           "axis_T" : [400.0, 750.0], 
           "imgfile" : "Meshram_Fig_05_HT.png"
        },        
        "zigzag Low T" :
        {
            "axis_T" : [300.0, 650.0], 
            "imgfile" : "Meshram_Fig_05_LT.png" 
        },
        "Straight High T" :
        {
            "axis_T" : [400.0, 750.0], 
            "imgfile" : "Meshram_Fig_04_HT.png"
        },
        "Straigth Low T" :
        {
            "axis_T" : [300.0, 650.0], 
            "imgfile" : "Meshram_Fig_04_LT.png"
        }
    }    

    exp.plot_results(ds_exp, plot_metacfg_base, plot_metacfg_offest)

    print('All done!')

###
if __name__ == "__main__":
    main()