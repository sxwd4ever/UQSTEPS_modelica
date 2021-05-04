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
from experiments import Experiment, PCHEExperiment, ErrorFunc, ParamFitting
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
        
    def post_process(self, test: TestItem, data_ref=None):
        super().post_process(test, data_ref)

        if data_ref == None:
            return

        # revese hc
        hc_hot = test.get_post_data('hc_hot').tolist()
        hc_cold = test.get_post_data('hc_cold').tolist()

        test.set_post_data('hc_hot', np.array(hc_hot[::-1]))
        test.set_post_data('hc_cold', np.array(hc_cold[::-1]))
        data_ref_cur = data_ref[test.name]
        keys = data_ref_cur['keys']
        row_names = ['1D', 'OEM']

        err_dict = {}

        for rname in row_names:
            row = data_ref_cur[rname]

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

class MarchionniExpFitting(ParamFitting):

    def prepare_data(self, values):
        values_new = deepcopy(values)
        keys = ['T_hot', 'dp_hot', 'T_cold', 'dp_cold']


        # extract 11 points out of 21 results to align with data referenced
        for k in keys:
            values_new[k] = values[k][0::2]

        return values_new
class ExpType(Enum):
    LOAD_PRE_EXP = 1,
    VS_CFD = 2,  # "aginst marchionni's 1D_vs_3D"
    FULL_SCALE = 3, # full scale off design in Sec. 4 of [Marchionni 2019]
    Fitting = 4 # paramters fitting

def cal_rho_bar(cfg_ref) -> Tuple[float, float]:
    '''
        rho bar calculation according to the boundary conditions set in cfg
    '''
    rho_bar = [1.0, 1.0]    

    T_flow = [(cfg_ref['T_hot_in'], cfg_ref['T_hot_out']), (cfg_ref['T_cold_in'], cfg_ref['T_cold_out'])]
    p = [cfg_ref['p_hot_in'], cfg_ref['p_cold_in']]
    for i in range(0, len(rho_bar)):
        (T_in, T_out) = T_flow[i]
        rho_bar[i] = PropsSI("DMASS", "P", p[i], "T" , (T_in + T_out) / 2, 'CO2')

    return tuple(rho_bar)

def main(work_root = []):
    # root path of modelica root
    if work_root == []:
        work_root = os.path.abspath(os.curdir)  

    # parameters initialization for this simulation
    exp_type:ExpType = ExpType.FULL_SCALE
    model_name = "Steps.Test.TPComponents.TestTP_PCHE_Marchionni"  
    # exp flags 
    use_rho_bar = -1.0 # > 1.0 use rho_bar for dp calculation 

    exp_name = 'Test-Marchionni_{}{} {:%Y-%m-%d-%H-%M-%S}'.format(exp_type.name, '_avg_rho_dp' if use_rho_bar > 0 else '', datetime.datetime.now())            
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
    sim_ops=[
        'startTime=0', 
        'stopTime=10',
        'stepSize=2',
        'solver=dassl',
        '-nls=homotopy',
        '-lv=LOG_DEBUG,LOG_INIT,LOG_NLS,LOG_NLS_V,LOG_STATS']

    if exp_type == ExpType.VS_CFD:
        single_run = True
        # referred base cfg
        cfg_base = {
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
        }         

        rho_bar = cal_rho_bar(cfg_base)
        cfg_offset = {}
        
        if not single_run: # in the parameter sweep way
            
            cfg_offset_base = {
                    "keys" : ["a_phi"],   
                    "a_phi=0" : [-1],
                    "a_phi=5" : [5.0]      
            }

            kc_cfs = [1, 1.2, 1.5, 1.8]

            # cfg with varied parameters from the base cfg
            for kc_cf in kc_cfs:
                cfg_offset_cp = deepcopy(cfg_offset_base)
                cfg_offset_cp['keys'].extend(['Cf_C1','Cf_C2','use_rho_bar', 'rho_bar_hot', 'rho_bar_cold'])
                for k, v in cfg_offset_cp.items():
                    if k == 'keys':
                        continue
                    v.extend([kc_cf, kc_cf * 2, use_rho_bar])
                    v.extend(rho_bar)

                cfg_offset[f'kc_cf={kc_cf}'] = cfg_offset_cp
        
        else: # sigle run to see result
            cfg_offset_base = {
                    "keys" : ['a_phi', 'Cf_C1','Cf_C2', 'use_rho_bar'],   
                    "High T" : [-1, 1.504692615, -1, use_rho_bar],
                    "Low T" : [-1, 1.660627977, -1, use_rho_bar]      
            }

            cfg_offset_cp = deepcopy(cfg_offset_base)
            cfg_offset_cp['keys'].extend(['rho_bar_hot', 'rho_bar_cold'])

            for k, v in cfg_offset_cp.items():
                if k == 'keys':
                    continue
                v.extend(rho_bar)

            cfg_offset[f'single_run'] = cfg_offset_cp

        ds_exp = TestDataSet.gen_test_dataset(cfg_base, cfg_offset, exp_name)
        # add referred data for comparison
        ds_exp.set_ref_data(ref_data)        
        ds_exp.add_view(mapping)

        json_str = ds_exp.to_json()    
        print(json_str)    

        exp.simulate_batch(
            sim_ops = sim_ops,
            solution_dict=exp.gen_sol_dict(cfg_base), 
            ds_test=ds_exp,
            append_save=False)   

        exp.save_results(ds_exp)
        exp.plot_results(ds_exp, plot_metacfg_base, plot_metacfg_offest)        

    elif exp_type == ExpType.FULL_SCALE:
        N_ch = 54 * 42
        A_c = 1.57e-6 # mm^2 -> m^2
        A_stack = N_ch * A_c # area of cross section for a stack containg all channels
        # referred base cfg
        cfg_base = {
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
            "G_hot_in": 2.06 / A_stack,  # "inlet mass flux rate kg/(m^2 s)";            
            "p_hot_in": from_bar(75), # "hot inlet pressure";
            "p_cold_in":from_bar(125), # "cold inlet pressure";
            "T_hot_in": from_degC(344.3), # "hot inlet temperature, K";
            "T_hot_out": from_degC(81), # "cold outlet temperature, K";
            "T_cold_in": from_degC(72.9), # "cold inlet temperature, K";
            "T_cold_out": from_degC(283), # "cold outlet temperature, K";
            # "kc_dp": 1 # "pressure drop correction coefficient"
        } 

        cfg_offset_base = {
                "keys" : ["G_hot_in", "T_cold_in"],   
                "Design point" : [2.06 / A_stack, from_degC(72.9)],
                "off-design #1" : [1.57 / A_stack, from_degC(72.9)],
                "off-design #2" : [2.09 / A_stack, from_degC(87.5)],                  
                "off-design #3" : [2.09 / A_stack,  from_degC(62)],
                "off-design #4" : [2.62 / A_stack, from_degC(72.9)]        
        }
        sweep_on_phi = False # parameter sweep on phi or kc_cfs, alter Modelica codes accordingly

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
            # kc_cfs_1 = [7, 7.3, 7.5, 7.8, 8] 
            kc_cfs_1 = [7, 7.3]  # quick version
            kc_cfs = []
            kc_cfs.extend(kc_cfs_1)
            # kc_cfs.extend(kc_cfs_1 * 2)

            cfg_offset = {} 
            l = []            
            for kc_cf in kc_cfs:
                cfg_offset_cp = deepcopy(cfg_offset_base)
                cfg_offset_cp['keys'].extend(['Cf_C1_hot','Cf_C1_cold'])
                for k, v in cfg_offset_cp.items():
                    if k == 'keys':
                        continue
                    v.extend([kc_cf, kc_cf * 2])

                cfg_offset[f'kc_cf={kc_cf}'] = cfg_offset_cp

        ds_exp = TestDataSet.gen_test_dataset(cfg_base, cfg_offset, ds_name=exp_name)        

        # add referred data for comparison
        ds_exp.set_ref_data(ref_data)
        ds_exp.add_view(mapping)

        json_str = ds_exp.to_json()    
        print(json_str)    

        exp.simulate_batch(
            sim_ops=[
                'startTime=0', 
                'stopTime=10',
                'stepSize=2',
                'solver=dassl',
                '-nls=homotopy',
                '-lv=LOG_DEBUG,LOG_INIT,LOG_NLS,LOG_NLS_V,LOG_STATS'],
            solution_dict=exp.gen_sol_dict(cfg_base), 
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

    elif exp_type == ExpType.Fitting:

        data_ref_full = {
            # 	0	0.016	0.032	0.048	0.064	0.08	0.096	0.112	0.128	0.144	0.16
            # Temperature [K]
            "T_hs_3D" : [673.9419, 650.9382, 630.6764, 610.4142, 588.7818, 569.0686, 549.6293, 531.0123, 512.6700, 495.4240, 479.2734],
            "T_hs_1D" : [670.3775, 650.3897, 630.6760, 610.4139, 592.0714, 573.1803, 554.5635, 536.4951, 518.9749, 502.0031, 487.7732],
            "T_cs_3D" : [528.0776, 507.2670, 487.8280, 469.2110, 452.2394, 435.8161, 420.7637, 407.3566, 394.4977, 383.2836, 374.2622],
            "T_cs_1D" : [519.8520, 501.5092, 483.4409, 465.9206, 449.7713, 433.8962, 419.3922, 406.2588, 393.9483, 383.0084, 374.8103],
            # (accumulated) pressdure drop, [kPa]
            "dp_hs_3D" : [0.0000, 0.6471, 1.3399, 2.0785, 2.8552, 3.6776, 4.5381, 5.4520, 6.4040, 7.4093, 8.3384],
            "dp_hs_1D" : [0.0000, 1.0968, 2.1402, 3.1380, 4.0976, 5.0038, 5.8720, 6.6944, 7.4787, 8.2248, 8.9253],
            "dp_cs_3D" : [3.7729, 3.4520, 3.1463, 2.8178, 2.4665, 2.0998, 1.7027, 1.2904, 0.8399, 0.3437, 0.0000],
            "dp_cs_1D" : [4.0549, 3.5358, 3.0473, 2.5968, 2.1463, 1.7340, 1.3445, 0.9779, 0.6189, 0.3056, 0.0152],
            "h_hs_3D" : [1.9840, 1.9608, 1.9719, 1.9926, 2.0202, 2.0382, 2.0589, 2.0988, 2.0976, 2.1307, 2.8428],
            "h_hs_1D" : [1.8852, 1.9032, 1.9212, 1.9419, 1.9640, 1.9847, 2.0081, 2.0302, 2.0537, 2.0771, 2.0896],
            "h_cs_3D" : [3.1872, 2.1268, 2.0529, 1.9748, 1.9352, 1.9024, 1.8929, 1.9013, 1.8973, 1.9249, 1.9593],
            "h_cs_1D" : [2.7152, 2.6248, 2.4727, 2.3480, 2.2521, 2.1768, 2.1248, 2.0879, 2.0702, 2.0579, 2.0581]}

        cfg_offset_full = {
                "keys" : ['a_phi', 'use_rho_bar'],   
                "Fitting_CFD" : [-1, use_rho_bar]
        }
        # referred base cfg
        cfg_base = {
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
        }   

        suffix_train = [("Fitting_CFD", "3D")] # use Marchionni's CFD result for params fitting

        case_dict = {
            name: {
                "cfg": {
                    "keys": [x for x in cfg_base.keys()] + cfg_offset_full['keys'],
                    "train": [x for x in cfg_base.values()] + cfg_offset_full[name]
                },
                "data_ref": {
                    "train": {
                        "T_hs": data_ref_full["T_hs_" + suffix],
                        "dp_hs": data_ref_full["dp_hs_" + suffix],
                        "T_cs": data_ref_full["T_cs_" + suffix],
                        "dp_cs": data_ref_full["dp_cs_" + suffix]
                    }
                }}
            for (name, suffix) in suffix_train}

        dim_y = 1

        for l in range(0, dim_y):
            # add mapping for calculated error
            # src -> dest
            k = f'err_fun_{l}'
            mapping["Comparison of T and dp"][k] = k

        for name, case in case_dict.items():
            fitting = MarchionniExpFitting(
                exp, case["cfg"], mapping, sim_ops, case["data_ref"], errfunc=ErrorFunc.Dp)

            fitting.run_fitting(pt=(1,1,1),max_steps=50)

    print('All done!')

###
if __name__ == "__main__":
    main()