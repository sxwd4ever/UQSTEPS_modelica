from collections import OrderedDict
from copy import deepcopy
import datetime
import itertools
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
from ex.random_local_search import random_local_search_2d
from model import TestDataSet,TestConfig, TestItem,Variable, TestConstants
from plotlib import PlotManager, DataSeries, AxisType
from physics import Temperature, Pressure, MDot
from experiments import Experiment, PCHEExperiment
from utils import ExcelHelper, from_degC, from_bar,mkdir_filepath

class MeshramTest(PCHEExperiment):
    '''
    script for test against Meshram [2016]
    '''

    def gen_sol_dict(self, cfg: dict) -> dict:
        sol_dict = super().gen_sol_dict(cfg)

        # specilized solutions
        # for KimPCHEHeatTransferFV
        # var_sp = [               
        #     Variable('HE.fluidFlow.heatTransfer.kim_cor_a.y'),
        #     Variable('HE.fluidFlow.heatTransfer.kim_cor_b.y'),
        #     Variable('HE.fluidFlow.heatTransfer.kim_cor_c.y'),
        #     Variable('HE.fluidFlow.heatTransfer.kim_cor_d.y'),
        #     Variable('HE.gasFlow.heatTransfer.kim_cor_a.y'),
        #     Variable('HE.gasFlow.heatTransfer.kim_cor_b.y'),
        #     Variable('HE.gasFlow.heatTransfer.kim_cor_c.y'),
        #     Variable('HE.gasFlow.heatTransfer.kim_cor_d.y')
        # ]   
        # for MarchionniPCHEHeatTransferFV
        var_sp = [
            Variable('HE.fluidFlow.heatTransfer.C1'),
            Variable('HE.fluidFlow.heatTransfer.C2'),
            Variable('HE.gasFlow.heatTransfer.C1'),
            Variable('HE.gasFlow.heatTransfer.C2')    
        ]

        for v in var_sp:
            sol_dict[v.key] = v 

        return sol_dict

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
                        'value': values['x'],
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
                        'value': values['x'],
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
                        'value': values['x'],
                        'axis_range': meta_cfg['axis_x'],
                        'label': meta_cfg['label_x']
                    },
                    'y':
                    {
                        'value': values['dp_hot'],
                        'axis_range': meta_cfg['axis_dp'],
                        'label': meta_cfg['label_dp']
                    },
                    "ax_type" : AxisType.Secondary
                },
                {
                    'name' : 'dp_cold',
                    'cs': 'b--^',
                    'x':
                    {
                        'value': values['x'],
                        'axis_range': meta_cfg['axis_x'],
                        'label': meta_cfg['label_x']
                    },
                    'y':
                    {
                        'value': values['dp_cold'],
                        'axis_range': meta_cfg['axis_dp'],
                        'label': meta_cfg['label_dp']
                    },
                    "ax_type" : AxisType.Secondary
                }],
                "imgfile" : meta_cfg["imgfile"]
            }                            
        }

class ExpType(Enum):
    LOAD_PRE_EXP = 1,
    VS_CFD = 2,  # "aginst mesharm's 1D_vs_3D"
    FITTING = 3, # Hyper parameters fitting with CFD data

def main(work_root = []):
    # root path of modelica root
    if work_root == []:
        work_root = os.path.abspath(os.curdir)  

    exp_type = ExpType.FITTING

    same_kc_cf = True # if kc_cf for hs and cs are identical or kc_cf_cold = 2 * kc_cf_hot
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
        "kc_dp": 2
    }

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

    # modelica simulation options
    sim_ops=[
        'startTime=0', 
        'stopTime=10',
        'stepSize=2',
        'solver=dassl',
        '-nls=homotopy',
        '-lv=LOG_DEBUG,LOG_INIT,LOG_NLS,LOG_NLS_V,LOG_STATS']

    if exp_type == ExpType.LOAD_PRE_EXP or exp_type == ExpType.VS_CFD:
        # cfg with varied parameters from the base cfg
        # for 'Zigzag Low T', I can not retrieve values match Meshram's Fig. 5 (a) with
        # u_hot_in = 1.345, which makes mdot_hot_in = 1.614 << mdot_cold_in = 5.488. 
        # Based on the facts that for other cases, mdot_hot_in ~ mdot_cold_in, 
        # I choose mdot_hot_in = 4.501 according to u_cold_in/u_hot_in = 5.584 
        # in 'Straight Low T' the senario and get good match    

        kc_cf_z_set = [10, 12, 15]    
        kc_cf_s_set = [1.0, 1.2, 1.5]
    
        cfg_offset_base = {
                # values in Table 3 in Meshram [2016]
                # Only with kc_dp = 1 at High T and kc_dp = 2 at Low T can I get good agreement with 
                # Meshram's DP result. 
                "keys" : ["T_cold_in", "T_cold_out", "T_hot_in", "T_hot_out", "u_cold_in", "u_hot_in", "a_phi","L_fp"],
                "Zigzag High T" : [500, 639.15, 730, 576.69, 1.876, 7.564, 36.0, 0.16],
                "Zigzag Low T": [400, 522.23, 630, 466.69, 0.806, 4.501, 36.0, 0.16], 
                "Straight High T": [500, 615.48, 730, 601.83, 1.518, 6.118, 0, 0.2],
                "Straigth Low T": [400, 498.45, 630, 494.37, 0.842, 4.702, 0, 0.2]
        }
        tests = ["Zigzag High T", "Zigzag Low T", "Straight High T", "Straigth Low T"]
        cfg_offset = {} 
        for (kc_cf_z, kc_cf_s) in itertools.product(kc_cf_z_set, kc_cf_s_set):
            cfg_offset_cp = deepcopy(cfg_offset_base)

            cfg_offset_cp['keys'].extend(['kc_cf_hot','kc_cf_cold'])
            
            for test in tests:
                kc_cf = kc_cf_z
                if "Straight" in test:
                    kc_cf = kc_cf_s
                # add kc_cf_hot, kc_cf_cold
                cfg_offset_cp[test].extend([kc_cf, kc_cf * (1 if same_kc_cf else 2)])

            cfg_offset[f'Z(kc_cf={kc_cf_z}) S(kc_cf={kc_cf_s})'] = cfg_offset_cp
                                
        mapping = {} 

        if exp_type == ExpType.VS_CFD:
            exp_name = 'Test-Meshram{} {:%Y-%m-%d-%H-%M-%S}'.format('_same_kc_cf' if same_kc_cf else '_diff_kc_cf', datetime.datetime.now())            
            ds_exp = TestDataSet.gen_test_dataset(cfg_ref, cfg_offset, ds_name=exp_name)
        
            ds_exp.add_view(mapping)

            json_str = ds_exp.to_json()    
            print(json_str)

            exp.simulate(
                sim_ops=sim_ops,
                solution_dict=exp.gen_sol_dict(cfg_ref), 
                ds_test=ds_exp,
                append_save=False)   

            exp.save_results(ds_exp)
        elif exp_type == ExpType.LOAD_PRE_EXP:
            # should search the exp_name 
            exp_name = "Test-Meshram_same_kc_cf 2021-01-27-14-40-51"
            ds_exp = exp.load_results(exp_name, dir_name=exp_name)
        else:
            raise ValueError("Invalid Experiment type")

        plot_metacfg_base = {
                "axis_x" : [0, 200e-3],
                "axis_T" : [400.0, 750.0],   
                "axis_dp" : [0, 80],
                "label_x" : 'X (m)',
                "label_T" : 'T (K)',
                "label_dp" : '$\\Delta~p~(kPa)$'
            }

        plot_metacfg_offest = {
            "Zigzag High T" :
            {
                "axis_x" : [0, 160e-3],
                "axis_T" : [400.0, 750.0], 
                "imgfile" : "Meshram_Fig_05_HT.png"
            },        
            "Zigzag Low T" :
            {
                "axis_x" : [0, 160e-3],            
                "axis_T" : [300.0, 650.0], 
                "imgfile" : "Meshram_Fig_05_LT.png" 
            },
            "Straight High T" :
            {
                "axis_T" : [400.0, 750.0],
                "axis_dp" : [0, 5], 
                "imgfile" : "Meshram_Fig_04_HT.png"
            },
            "Straigth Low T" :
            {
                "axis_T" : [300.0, 650.0], 
                "axis_dp" : [0, 5], 
                "imgfile" : "Meshram_Fig_04_LT.png"
            }
        }    

        exp.plot_results(ds_exp, plot_metacfg_base, plot_metacfg_offest)

        print('All done!')

    elif exp_type == ExpType.FITTING:        

        (d_kc_h, d_kc_c) = (0,0)

        test = "Straight High T"  

        cfg_offset_base = {
                # values in Table 3 in Meshram [2016]
                # Only with kc_dp = 1 at High T and kc_dp = 2 at Low T can I get good agreement with 
                # Meshram's DP result. 
                "keys" : ["T_cold_in", "T_cold_out", "T_hot_in", "T_hot_out", "u_cold_in", "u_hot_in", "a_phi","L_fp"],
                test: [500, 615.48, 730, 601.83, 1.518, 6.118, 0, 0.2],
        }        

        fitting  = ParamFitting(exp, cfg_ref, cfg_offset_base, {}, sim_ops, test)

        (h_pt, h_eval) = param_random_local_search(
            func_para=fitting.evaluate, 
            pt=(0.1, 0.1), 
            max_steps=20,
            num_samples=2, 
            steplength=0.5)

class ParamFitting(object):

    def __init__(self, exp, cfg_ref, cfg_offset_base, mapping, sim_ops, test) -> None:
        super().__init__() 
        self.exp = exp
        self.cfg_ref = cfg_ref
        self.cfg_offset_base = cfg_offset_base
        self.mapping = mapping
        self.sim_ops = sim_ops
        self.data_ref = self.gen_reference_data()
        self.test = test

    def evaluate(self, x) -> Tuple[float,float]:

        (d_kc_h, d_kc_c) = x

        ds_exp = self.gen_experiment_ds(1 + d_kc_h, 1+ d_kc_c, self.cfg_ref, self.cfg_offset_base,self.mapping)
        sol_dict=self.exp.gen_sol_dict(self.cfg_ref)

        self.exp.simulate(
            sim_ops=self.sim_ops,
            solution_dict=sol_dict, 
            ds_test=ds_exp,
            append_save=False)  

        values = None

        # calculate errors
        # only 1 test suppose to have
        for g in ds_exp.values():
            for test in g.values():
                values = test.post_data
                # only 1 test suppose to have
                break
            break

        T_hs = values['T_hot']
        dp_hs = values['dp_hot']
        T_cs = values['T_cold']
        dp_cs = values['dp_cold']

        data_ref = self.gen_reference_data()['train']

        err_h = 0.0
        err_c = 0.0

        for i in range(0, len(T_hs)):
            
            err_h += \
                self.cal_err(T_hs[i],data_ref['T_hs'][i]) + \
                self.cal_err(dp_hs[i], data_ref['dp_hs'][i])

            err_c += \
                self.cal_err(T_cs[i], data_ref['T_cs'][i]) + \
                self.cal_err(dp_cs[i], data_ref['dp_cs'][i])

        return (err_h, err_c)
    
    def cal_err(self, x, y):
        if abs(y) < 1e-10:
            return 0
        
        return ((x - y) / y) ** 2

    def gen_experiment_ds(self, kc_cf_h, kc_cf_c, cfg_ref, cfg_offset_base, mapping) -> TestDataSet:
              

        cfg_offset = {} 
        cfg_offset_cp = deepcopy(cfg_offset_base)

        cfg_offset_cp['keys'].extend(['kc_cf_hot','kc_cf_cold'])
        
        # add kc_cf_hot, kc_cf_cold
        cfg_offset_cp[self.test].extend([kc_cf_h, kc_cf_c])

        cfg_offset[f'kc_cf_h={kc_cf_h} kc_cf_c={kc_cf_c}'] = cfg_offset_cp

        exp_name = 'Test-Meshram_diff_kc_cf {:%Y-%m-%d-%H-%M-%S}'.format(datetime.datetime.now())            
        ds_exp = TestDataSet.gen_test_dataset(cfg_ref, cfg_offset, ds_name=exp_name)

        ds_exp.add_view(mapping)  

        return ds_exp

    def gen_reference_data(self) ->dict :

        data_ref_full = {
                # 	0	0.016	0.032	0.048	0.064	0.08	0.096	0.112	0.128	0.144	0.16
                "T_hs_ZHT" : [730.0439, 714.2321, 697.9086, 681.5861, 665.7743, 649.9631, 634.663, 619.8756, 604.5755, 590.2993, 577.0465],
                "T_hs_ZLT" : [630.0433, 611.1623, 592.792, 574.4221, 556.5635, 539.2176, 523.9174, 508.1066, 493.3176, 480.0652, 466.8118],
                "T_hs_SHT" : [730.2377, 716.1961, 702.6746, 689.6731, 676.6716, 664.1902, 651.1887, 639.2273, 626.7459, 614.7845, 602.3031],
                "T_hs_SLT" : [630.4124, 614.433, 599.4845, 585.567, 571.134, 557.732, 544.3299, 531.4433, 519.0722, 506.701, 494.8454],
                "T_cs_ZHT" : [639.4728, 624.1737, 609.3853, 594.0863, 578.7862, 565.0217, 551.2567, 537.4927, 524.2394, 512.01, 500.8034],
                "T_cs_ZLT" : [522.5872, 506.7758, 490.9645, 476.1765, 461.388, 448.6474, 437.4408, 426.7464, 417.0754, 408.9395, 400.803],
                "T_cs_SHT" : [616.3447, 602.8232, 590.3418, 578.3804, 566.419, 554.4577, 543.0163, 532.6152, 522.214, 511.2927, 500.8915],
                "T_cs_SLT" : [498.9691, 486.5979, 474.2268, 463.4021, 452.0619, 442.268, 432.9897, 424.2268, 415.9794, 408.2474, 401.0309],
                "dp_hs_ZHT" : [0.1168, 3.0365, 7.4745, 12.146, 15.5328, 20.438, 25.4599, 29.1971, 34.1022, 40.4088, 45.8978],
                "dp_hs_ZLT" : [-0.117, 2.1138, 5.629, 9.1445, 11.8422, 15.7079, 19.6903, 22.7387, 26.7211, 31.7548, 36.5547],
                "dp_hs_SHT" : [-0.0074, 0.2901, 0.6098, 0.9295, 1.2715, 1.6135, 1.9703, 2.3346, 2.6988, 3.0853, 3.5905],
                "dp_hs_SLT" : [0, 0.2132, 0.4412, 0.6912, 0.9412, 1.1912, 1.4706, 1.7426, 2.0368, 2.3309, 2.7353],
                "dp_cs_ZHT" : [16.1168, 14.5985, 12.9635, 11.4453, 10.1606, 8.5255, 6.6569, 5.3723, 3.5036, 1.5182, 0.1168],
                "dp_cs_ZLT" : [11.2114, 10.0554, 9.0159, 7.9765, 7.2875, 6.1313, 4.975, 4.0526, 2.5459, 1.1562, 0.1168],
                "dp_cs_SHT" : [1.4169, 1.2693, 1.1439, 1.0185, 0.8858, 0.7455, 0.6053, 0.4577, 0.3101, 0.155, 0],
                "dp_cs_SLT" : [1.0147, 0.9191, 0.8309, 0.7426, 0.6544, 0.5515, 0.4559, 0.3529, 0.2426, 0.125, 0]
            }

        # High T range (train) v.s. low T range (test)
        # straight channel PCHE
        suffix_train = "SHT"
        suffix_test = "SLT"

        data_ref = {
            "train":{
                "T_hs" : data_ref_full["T_hs_" + suffix_train],
                "dp_hs" : data_ref_full["dp_hs_" + suffix_train],
                "T_cs" : data_ref_full["T_cs_" + suffix_train],
                "dp_cs" : data_ref_full["dp_cs_" + suffix_train]            
            },
            "test": {
                "T_hs" : data_ref_full["T_hs_" + suffix_test],
                "dp_hs" : data_ref_full["dp_hs_" + suffix_test],
                "T_cs" : data_ref_full["T_cs_" + suffix_test],
                "dp_cs" : data_ref_full["dp_cs_" + suffix_test] 
            }
        }

        return data_ref

def param_random_local_search(func_para,pt,max_steps,num_samples,steplength):
    '''
        Concurrently train C1 for hot/cold side due to simulation is a time consuming process
    '''
    
    # starting point evaluation
    current_eval = func_para(pt)
    current_pt = pt
    
    # loop over max_its descend until no improvement or max_its reached
    pt_history = [current_pt]
    eval_history = [current_eval]
    for i in range(max_steps):
        # loop over num_samples, randomly sample direction and evaluate, move to best evaluation
        swap = 0
        keeper_pt = current_pt
        
        # check if diminishing steplength rule used
        if steplength == 'diminish':
            steplength_temp = 1/(1 + i)
        else:
            steplength_temp = steplength
        
        for j in range(num_samples):            
            # produce direction
            import random as rnd
            import math

            # theta = np.random.rand(1)
            theta = rnd.uniform(0,1)

            # x = steplength_temp*np.cos(2*np.pi*theta)
            # y = steplength_temp*np.sin(2*np.pi*theta)
            x = steplength_temp*math.cos(2*math.pi*theta)
            y = steplength_temp*math.sin(2*math.pi*theta)
            new_pt = np.asarray([x,y])
            temp_pt = deepcopy(keeper_pt) 
            new_pt += temp_pt
            
            # evaluate new point
            new_eval = func_para(new_pt)
            x_new = [current_pt[0], current_pt[1]]
            y_old = [current_eval[0], current_eval[1]]
            y_min = y_old
            y_new = [new_eval[0], new_eval[1]]
            
            # evaluate results respectively
            update = False
            for k in range(0, len(y_new)):
                if y_new[k] < y_old[k]:
                    x_new[k] = new_pt[k]
                    y_min[k] = y_new[k]
                    update = True

            if update:
                current_pt = (x_new[0], x_new[1])
                current_eval = (y_min[0], y_min[1])
                swap = 1
        
        # if nothing has changed
        if swap == 1:
            pt_history.append(current_pt)
            eval_history.append(current_eval)
    
    # translate to array, reshape appropriately
    pt_history = np.asarray(pt_history)
    # pt_history.shape = (np.shape(pt_history)[0],np.shape(pt_history)[1])

    eval_history = np.asarray(eval_history)
    # eval_history.shape = (np.shape(eval_history)[0])

    return pt_history,eval_history   

if __name__ == "__main__":
    main()