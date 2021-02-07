from collections import OrderedDict
from copy import deepcopy
import datetime
import itertools
from logging import exception
from math import trunc
from os.path import join, sep
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
from experiments import Experiment, PCHEExperiment, ErrorFunc, ParamFitting
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
        N_seg = cfg['N_seg']
        keys = [
            'HE.fluidFlow.heatTransfer.q1[%s]',
            'HE.fluidFlow.heatTransfer.q2[%s]',
            'HE.fluidFlow.heatTransfer.q3[%s]',
            # 'HE.fluidFlow.heatTransfer.C2[%s]',
            'HE.gasFlow.heatTransfer.q1[%s]',
            'HE.gasFlow.heatTransfer.q2[%s]',
            'HE.gasFlow.heatTransfer.q3[%s]',
            # 'HE.gasFlow.heatTransfer.C2[%s]'
        ]

        for k in keys:
            for i in range(1, N_seg+1):
                key = k % i
                sol_dict[key] = Variable(key)

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
                }
                ],
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

    exp_type = ExpType.VS_CFD

    model_name = "Steps.Test.TestTP_PCHE_Meshram"  
    # referred base cfg
    cfg_base = {
        # geometry parameters
        "N_ch": 1e5, # "channel number"
        "N_seg": 10, # "segments number"
        "D_ch": 2e-3, # "channel diameter, semi circular tube"
        "L_fp": 200e-3, # "channel flow path length"
        "L_pitch": 15e-3, # "pitch length"
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
        # "kc_dp": 2
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

    # data view mapping - for summury generation
    mapping = {
        "Calculated values": {
            "Re_bar_hot": "Re_bar_hot",
            "Re_bar_cold": "Re_bar_cold",
            "Pr_bar_hot": "Pr_bar_hot",
            "Pr_bar_cold": "Pr_bar_cold",            
            "mu_bar_hot": "mu_bar_hot",
            "mu_bar_cold": "mu_bar_cold",
            "cp_bar_hot": "cp_bar_hot",
            "cp_bar_cold": "cp_bar_cold",
            "k_bar_hot": "k_bar_hot",
            "k_bar_cold": "k_bar_cold",
            "Re_tilde_hot": "Re_tilde_hot",
            "Re_tilde_cold": "Re_tilde_cold",          
            "Pr_tilde_hot": "Pr_tilde_hot",
            "Pr_tilde_cold": "Pr_tilde_cold",
        }
    }

    if exp_type == ExpType.LOAD_PRE_EXP or exp_type == ExpType.VS_CFD:
        # cfg with varied parameters from the base cfg
        # for 'Zigzag Low T', I can not retrieve values match Meshram's Fig. 5 (a) with
        # u_hot_in = 1.345, which makes mdot_hot_in = 1.614 << mdot_cold_in = 5.488. 
        # Based on the facts that for other cases, mdot_hot_in ~ mdot_cold_in, 
        # I choose mdot_hot_in = 4.501 according to u_cold_in/u_hot_in = 5.584 
        # in 'Straight Low T' the senario and get good match    

        # comparison among different parameters
        if exp_type == ExpType.VS_CFD:

            cfg_offset = {} 
            mapping = {} 
            exp_name = ""

            # prepare config offset and experiment name accordingly
            single_run = True
            if not single_run:
                kc_cf_z_set = [10, 12]    
                kc_cf_s_set = [1.0, 1.2]        
            
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
                
                for (kc_cf_z, kc_cf_s) in itertools.product(kc_cf_z_set, kc_cf_s_set):
                    cfg_offset_cp = deepcopy(cfg_offset_base)

                    cfg_offset_cp['keys'].extend(['Cf_C1','Cf_C2'])
                    
                    for test in tests:
                        kc_cf = kc_cf_z
                        if "Straight" in test:
                            kc_cf = kc_cf_s
                        # add kc_cf_hot, kc_cf_cold
                        cfg_offset_cp[test].extend([kc_cf, kc_cf * 2])

                    cfg_offset[f'Z(kc_cf={kc_cf_z}) S(kc_cf={kc_cf_s})'] = cfg_offset_cp

                exp_name = 'Test-Meshram {:%Y-%m-%d-%H-%M-%S}'.format(datetime.datetime.now())            
                
            else: # single run 
                # Gnielinski cor for straight channel and Ngo's cor for zigzag channels
                # different C1, C2 applied for HT and LT respectively. 
                # For Straight Channel, C2 is a dummy variable since C2 = 1 / C1
                # cfg_offset[f"vs Meshram's CFD"] = {
                #         # values in Table 3 in Meshram [2016]
                #         "keys" : ["T_cold_in", "T_cold_out", "T_hot_in", "T_hot_out", "u_cold_in", "u_hot_in", "a_phi","L_fp", "Cf_C1", "Cf_C2"],
                #         "Zigzag High T" : [500, 639.15, 730, 576.69, 1.876, 7.564, 36.0, 0.16, 1.464907987, 1.382382861],
                #         "Zigzag Low T": [400, 522.23, 630, 466.69, 0.806, 4.501, 36.0, 0.16, 1.28633131, 2.647075995], 
                #         "Straight High T": [500, 615.48, 730, 601.83, 1.518, 6.118, -1, 0.2, 1.504692615,1],
                #         "Straigth Low T": [400, 498.45, 630, 494.37, 0.842, 4.702, -1, 0.2, 1.660627977, 1]
                # }

                # Gnielinski cor for straight channel and Liao's Nusselt number cor + Ngo's friction factor cor for zigzag channels
                # different C1, C2 applied for HT and LT respectively. 
                # For Straight Channel, C2 is a dummy variable since C2 = 1 / C1
                # cfg_offset[f"vs Meshram's CFD"] = {
                #         # values in Table 3 in Meshram [2016]
                #         "keys" : ["T_cold_in", "T_cold_out", "T_hot_in", "T_hot_out", "u_cold_in", "u_hot_in", "a_phi","L_fp", "Cf_C1", "Cf_C2"],
                #         "Zigzag High T" : [500, 639.15, 730, 576.69, 1.876, 7.564, 36.0, 0.16, 2.146423845, 1.382382861],
                #         "Zigzag Low T": [400, 522.23, 630, 466.69, 0.806, 4.501, 36.0, 0.16, 1.345870392, 2.647075995], 
                #         "Straight High T": [500, 615.48, 730, 601.83, 1.518, 6.118, -1, 0.2, 1.504692615, 1],
                #         "Straigth Low T": [400, 498.45, 630, 494.37, 0.842, 4.702, -1, 0.2, 1.660627977, 1]
                # }                

                # My correlation for zigzag PCHE
                G = 498.341 # should be calculated as G = mdot / A

                (C1_HT, C2_HT) = exp.predict_cor_coff(730, cfg_base['p_hot_in'], 500, cfg_base['p_cold_in'], G, cfg_base['D_ch'], "CO2")
                (C1_LT, C2_LT) = exp.predict_cor_coff(630, cfg_base['p_hot_in'], 400, cfg_base['p_cold_in'], G, cfg_base['D_ch'], "CO2")      

                cfg_offset[f"vs Meshram's CFD"] = {
                        # values in Table 3 in Meshram [2016]
                        "keys" : ["T_cold_in", "T_cold_out", "T_hot_in", "T_hot_out", "u_cold_in", "u_hot_in", "a_phi","L_fp", "Cf_C1", "Cf_C2"],
                        "Zigzag High T" : [500, 639.15, 730, 576.69, 1.876, 7.564, 36.0, 0.16, C1_HT, C2_HT],
                        "Zigzag Low T": [400, 522.23, 630, 466.69, 0.806, 4.501, 36.0, 0.16, C1_LT, C2_LT], 
                        "Straight High T": [500, 615.48, 730, 601.83, 1.518, 6.118, -1, 0.2, 1.504692615, 1],
                        "Straigth Low T": [400, 498.45, 630, 494.37, 0.842, 4.702, -1, 0.2, 1.660627977, 1]
                } 

                exp_name = 'Test-Meshram_single {:%Y-%m-%d-%H-%M-%S}'.format(datetime.datetime.now())            
            
            ds_exp = TestDataSet.gen_test_dataset(cfg_base, cfg_offset, ds_name=exp_name)

            ds_exp.add_view(mapping)

            json_str = ds_exp.to_json()    
            print(json_str)

            exp.simulate_batch(
                sim_ops=sim_ops,
                solution_dict=exp.gen_sol_dict(cfg_base), 
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

        cfg_offset_full = {
                # values in Table 3 in Meshram [2016]
                # Only with kc_dp = 1 at High T and kc_dp = 2 at Low T can I get good agreement with 
                # Meshram's DP result. 
                "keys" : ["T_cold_in", "T_cold_out", "T_hot_in", "T_hot_out", "u_cold_in", "u_hot_in", "a_phi","L_fp"],
                "Zigzag High T" : [500, 639.15, 730, 576.69, 1.876, 7.564, 36.0, 0.16],
                "Zigzag Low T": [400, 522.23, 630, 466.69, 0.806, 4.501, 36.0, 0.16], 
                "Straight High T": [500, 615.48, 730, 601.83, 1.518, 6.118, -1, 0.2],
                "Straigth Low T": [400, 498.45, 630, 494.37, 0.842, 4.702, -1, 0.2]
        }        

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

        suffix_straight = [
            ("Straight High T", "SHT"),
            ("Straigth Low T", "SLT")]

        suffix_zigzag = [
            ("Zigzag High T", "ZHT"),
            ("Zigzag Low T", "ZLT")]

        suffix_train_full = deepcopy(suffix_straight)
        suffix_train_full.extend(deepcopy(suffix_zigzag))

        suffix_train = suffix_zigzag 
        # suffix_train = suffix_train_full
        # suffix_train = suffix_straight 

        # case_dict = {
        #     name: {
        #         "cfg": {
        #             "keys": [x for x in cfg_base.keys()] + ["T_cold_in", "T_cold_out", "T_hot_in", "T_hot_out", "u_cold_in", "u_hot_in", "a_phi", "L_fp"],
        #             "train": [x for x in cfg_base.values()] + cfg_offset_full[name]
        #         },
        #         "data_ref": {
        #             "train": {
        #                 "T_hs": data_ref_full["T_hs_" + suffix],
        #                 "dp_hs": data_ref_full["dp_hs_" + suffix],
        #                 "T_cs": data_ref_full["T_cs_" + suffix],
        #                 "dp_cs": data_ref_full["dp_cs_" + suffix]
        #             }
        #         }}
        #     for (name, suffix) in suffix_train}

        # create cfg table containing all tests in one, each row contains values for one cfg instance
        # cfg_table = {
        #     "keys": [x for x in cfg_base.keys()] + ["T_cold_in", "T_cold_out", "T_hot_in", "T_hot_out", "u_cold_in", "u_hot_in", "a_phi", "L_fp"]
        # }
        # # append all the cfg rows 
        # cfg_table.update({
        #         name : [x for x in cfg_base.values()] + cfg_offset_full[name]
        #     for (name, suffix) in suffix_train})

        # data_ref = {
        #     name: {
        #         "T_hs": data_ref_full["T_hs_" + suffix],
        #         "dp_hs": data_ref_full["dp_hs_" + suffix],
        #         "T_cs": data_ref_full["T_cs_" + suffix],
        #         "dp_cs": data_ref_full["dp_cs_" + suffix]
        #     }
        #     for (name, suffix) in suffix_train}

        # case_dict = {
        #     "cfg": cfg_table,
        #     "data_ref": data_ref
        # }

        # separate test cases into two groups
        case_dict = {
            name:{
                "cfg": {
                    "keys": [x for x in cfg_base.keys()] + ["T_cold_in", "T_cold_out", "T_hot_in", "T_hot_out", "u_cold_in", "u_hot_in", "a_phi", "L_fp"],
                    name : [x for x in cfg_base.values()] + cfg_offset_full[name]
                },
                "data_ref": {
                    name: {
                        "T_hs": data_ref_full["T_hs_" + suffix],
                        "dp_hs": data_ref_full["dp_hs_" + suffix],
                        "T_cs": data_ref_full["T_cs_" + suffix],
                        "dp_cs": data_ref_full["dp_cs_" + suffix]
                    }
                }
            }
            for (name, suffix) in suffix_train
        }
            
        dim_y = 1

        for l in range(0, dim_y):
            # add mapping for calculated error
            # src -> dest
            k = f'err_fun_{l}'
            mapping["Calculated values"][k] = k

        for case in case_dict.values():

            fitting = ParamFitting(
                exp, case, mapping, sim_ops, errfunc=ErrorFunc.T)

            fitting.run_fitting(pt=(1,1,1), max_steps=50)

if __name__ == "__main__":
    main()
