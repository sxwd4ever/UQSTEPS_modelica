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

import CoolProp as CP
from CoolProp.CoolProp import PhaseSI, PropsSI, get_global_param_string
import CoolProp.CoolProp as CoolProp

import os
from os import path, sep

import numpy as np
from model import TestDataSet,TestConfig, TestItem, TestResult,Variable
from plotlib import PlotManager, DataSeries, AxisType
from physics import Temperature, Pressure, MDot
from experiments import Experiment, PCHEExperiment
from utils import from_degC, from_bar, from_kPa

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
            Variable('dp_cold_act', 'pa'),          
            Variable('rho_bar_hot'),
            Variable('rho_bar_cold')
            # Variable('dp_hot_act_m', 'pa'),
            # Variable('dp_cold_act_m', 'pa')
        ]

        for v in var_sp:
            sol_dict[v.key] = v

        return sol_dict
        
    def post_process(self, test: TestItem, ds_exp:TestDataSet):
        super().post_process(test, ds_exp)

        if ds_exp.ref_data == None:
            return

        # revese hc
        hc_hot = test.get_post_data('hc_hot').tolist()
        hc_cold = test.get_post_data('hc_cold').tolist()

        test.set_post_data('hc_hot', np.array(hc_hot[::-1]))
        test.set_post_data('hc_cold', np.array(hc_cold[::-1]))

        if not test.name in ds_exp.ref_data.keys():
            return

        ref_data = ds_exp.ref_data[test.name]

        keys = ref_data['keys']
        row_names = ['1D', 'OEM']

        err_dict = {}

        for rname in row_names:
            row = ref_data[rname]

            for i in range(0, len(keys)):
                k = keys[i] + "_act"
                val = test.result[k].val
                err_dict[f'{k}_{rname}_err'] = (val - row[i]) / row[i] * 100

        for (k, v) in err_dict.items():
            test.result[k] = Variable(k, '%', val=v)

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
                        'value': values['x'] * 1e3,
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
                        'value': values['x'] * 1e3,
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
                        'value': values['x'] * 1e3,
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
                        'value': values['x'] * 1e3,
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
    # exp flags 
    use_rho_bar = -1.0 # > 1.0 use rho_bar for dp calculation 
    same_kc_cf = False # if kc_cf for hs and cs are identical or kc_cf_cold = 2 * kc_cf_hot
    sweep_on_phi = False # parameter sweep on phi or kc_cfs, alter Modelica codes accordingly

    exp_name = 'Test-Marchionni_{}{}{} {:%Y-%m-%d-%H-%M-%S}'.format(exp_type.name, '_avg_rho_dp' if use_rho_bar > 0 else '', '_same_kc_cf' if same_kc_cf else '_diff_kc_cf', datetime.datetime.now())            
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

    # src->dst 
    mapping = {
        "Comparison of T and dp": {
            "T_hot_out_act": "T_hot_out_act",
            "T_cold_out_act": "T_cold_out_act",
            "T_hot_out_act_OEM_err": "T_hot_out_err",
            "T_cold_out_act_OEM_err": "T_cold_out_err",            
            "dp_hot_act": "dp_hot_act",
            "dp_cold_act": "dp_cold_act",
            "dp_hot_act_OEM_err": "dp_hot_err",
            "dp_cold_act_OEM_err": "dp_cold_err"
        }
    }

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

    ref_data= {
            "Design point" :{
                "keys" : ['dp_hot','T_hot_out', 'dp_cold', 'T_cold_out'],
                '1D': [from_kPa(131), from_degC(81.4), from_kPa(119), from_degC(283.0)],
                'OEM': [from_kPa(130), from_degC(80.5), from_kPa(120), from_degC(284.9)]
            },
            "off-design #1" :{
                "keys" : ['dp_hot','T_hot_out', 'dp_cold', 'T_cold_out'],
                '1D': [from_kPa(80), from_degC(79.5), from_kPa(73), from_degC(285.4)],
                'OEM': [from_kPa(79), from_degC(78.6), from_kPa(74), from_degC(287.2)]
            },            
            "off-design #2" :{
                "keys" : ['dp_hot','T_hot_out', 'dp_cold', 'T_cold_out'],
                '1D': [from_kPa(146), from_degC(99.4), from_kPa(138), from_degC(293.2)],
                'OEM': [from_kPa(145), from_degC(99.7), from_kPa(139), from_degC(294.5)]
            }, 
            "off-design #3" :{
                "keys" : ['dp_hot','T_hot_out', 'dp_cold', 'T_cold_out'],
                '1D': [from_kPa(123), from_degC(67.7), from_kPa(104), from_degC(267.7)],
                'OEM': [from_kPa(122), from_degC(66.6), from_kPa(106), from_degC(269.3)]
            },  
            "off-design #4" :{
                "keys" : ['dp_hot','T_hot_out', 'dp_cold', 'T_cold_out'],
                '1D': [from_kPa(205), from_degC(83.6), from_kPa(183), from_degC(280.1)],
                'OEM': [from_kPa(202), from_degC(82.7), from_kPa(184), from_degC(282.3)]
            }                       
        }

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
            # "kc_dp": 1 # "pressure drop correction coefficient"
        } 
        
        cfg_offset_base = {
                "keys" : ["a_phi"],   
                "a_phi=0" : [-1],
                "a_phi=5" : [5.0]      
        }

        rho_bar = [1.0, 1.0]

        if use_rho_bar:
            T_flow = [(cfg_ref['T_hot_in'], cfg_ref['T_hot_out']), (cfg_ref['T_cold_in'], cfg_ref['T_cold_out'])]
            p = [cfg_ref['p_hot_in'], cfg_ref['p_cold_in']]
            for i in range(0, len(rho_bar)):
                (T_in, T_out) = T_flow[i]
                rho_bar[i] = PropsSI("DMASS", "P", p[i], "T" , (T_in + T_out) / 2, 'CO2')

        kc_cfs = [1, 1.2, 1.5, 1.8]
        cfg_offset = {} 

        # cfg with varied parameters from the base cfg
        for kc_cf in kc_cfs:
            cfg_offset_cp = deepcopy(cfg_offset_base)
            cfg_offset_cp['keys'].extend(['Cf_C1','Cf_C2','use_rho_bar', 'rho_bar_hot', 'rho_bar_cold'])
            for k, v in cfg_offset_cp.items():
                if k == 'keys':
                    continue
                v.extend([kc_cf, kc_cf * (1 if same_kc_cf else 2), use_rho_bar])
                v.extend(rho_bar)

            cfg_offset[f'kc_cf={kc_cf}'] = cfg_offset_cp

        ds_exp = TestDataSet.gen_test_dataset(cfg_ref, cfg_offset, exp_name)
        # add referred data for comparison
        ds_exp.set_ref_data(ref_data)        
        ds_exp.add_view(mapping)

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
        N_ch = 54 * 42
        A_c = 1.57e-6 # mm^2 -> m^2
        A_stack = N_ch * A_c # area of cross section for a stack containg all channels
        # referred base cfg
        cfg_ref = {
            # geometry parameters
            "N_ch": N_ch, # "channel number"
            "N_seg": 20, # "segments number"
            "D_ch": 2e-3, # "channel diameter, semi circular tube"
            "L_fp": 1012e-3, # "channel flow path length"
            "L_pitch": 12.3e-3, # "pitch length"
            "a_phi": 45.0, # "pitch angle, degree"
            "H_ch": 3.26e-3, # "Height of the solid domain, containing one cold tube and one hot tube"
            "W_ch": 1.27e-3 * 2, # "Width of the solid domain"
            # boundary conditon
            "G_in": 2.06 / A_stack,  # "inlet mass flux rate kg/(m^2 s)";            
            "p_hot_in": from_bar(75), # "hot inlet pressure";
            "p_cold_in":from_bar(125), # "cold inlet pressure";
            "T_hot_in": from_degC(344.3), # "hot inlet temperature, K";
            "T_hot_out": from_degC(81), # "cold outlet temperature, K";
            "T_cold_in": from_degC(72.9), # "cold inlet temperature, K";
            "T_cold_out": from_degC(283), # "cold outlet temperature, K";
            # "kc_dp": 1 # "pressure drop correction coefficient"
        } 

        cfg_offset_base = {
                "keys" : ["G_in", "T_cold_in"],   
                "Design point" : [2.06 / A_stack, from_degC(72.9)],
                "off-design #1" : [1.57 / A_stack, from_degC(72.9)],
                "off-design #2" : [2.09 / A_stack, from_degC(87.5)],                  
                "off-design #3" : [2.09 / A_stack,  from_degC(62)],
                "off-design #4" : [2.62 / A_stack, from_degC(72.9)]        
        }

        if sweep_on_phi:
            # cfg with varied parameters from the base cfg
            phis = [5, 10, 15, 20, 25, 30, 35, 40, 45]
            # phis = [5, 10, 15]

            cfg_offset = {} 

            for phi in phis:
                cfg_offset_cp = deepcopy(cfg_offset_base)
                cfg_offset_cp['keys'].append('a_phi')
                for k, v in cfg_offset_cp.items():
                    if k == 'keys':
                        continue
                    v.append(phi)

                cfg_offset[f'phi={phi}'] = cfg_offset_cp
        else:
            # kc_cfs = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
            kc_cfs_1 = [7, 7.3, 7.5, 7.8, 8] 
            kc_cfs = []
            kc_cfs.extend(kc_cfs_1)
            kc_cfs.extend(kc_cfs_1 * 2)

            cfg_offset = {} 
            l = []            
            for kc_cf in kc_cfs:
                cfg_offset_cp = deepcopy(cfg_offset_base)
                cfg_offset_cp['keys'].extend(['kc_cf_hot','kc_cf_cold'])
                for k, v in cfg_offset_cp.items():
                    if k == 'keys':
                        continue
                    v.extend([kc_cf, kc_cf * (1 if same_kc_cf else 2)])

                cfg_offset[f'kc_cf={kc_cf}_{"same" if same_kc_cf else "diff"}'] = cfg_offset_cp

        ds_exp = TestDataSet.gen_test_dataset(cfg_ref, cfg_offset, ds_name=exp_name)        

        # add referred data for comparison
        ds_exp.set_ref_data(ref_data)
        ds_exp.add_view(mapping)

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
        exp_name = "Test-Marchionni_VS_CFD_diff_kc_cf 2021-01-26-16-00-50"
        ds_exp:TestDataSet = TestDataSet(name=exp_name, groups={})

        ds_exp.set_ref_data(ref_data)    
        exp.load_results(exp_name, dir_name=exp_name, ds_exp=ds_exp,has_view=True)
        exp.plot_results(ds_exp, plot_metacfg_base, plot_metacfg_offest)   

    print('All done!')

###
if __name__ == "__main__":
    main()