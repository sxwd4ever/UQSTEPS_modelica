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
from experiments import Experiments
from utils import ExcelHelper, from_degC, from_bar

class MeshramTest(Experiments):
    '''
    script for test against Meshram [2016]
    '''

    def gen_plot_manager(self, test: TestItem, plot_cfg:dict) -> PlotManager:
        cfg = test.cfg
        result = test.result

        N_seg = int(cfg['N_seg'])
        L_fp = cfg['L_fp']
        len_seg = L_fp / N_seg   

        axis_x = plot_cfg["axis_x"]
        axis_T = plot_cfg["axis_T"]
        axis_dp = plot_cfg["axis_dp"]
        label_x = plot_cfg["label_x"]
        label_T =  plot_cfg["label_T"]
        label_dp = plot_cfg["label_dp"]  

        T_hot = []
        T_cold = []
        # global pressure drop R. T. the inlet temperature 
        dp_hot = [0] 
        dp_cold = [0]
        # started, fix pressure of each stream
        dp_hot_cur = 0
        dp_cold_cur = 0

        for i in range(0, N_seg):
            idx = i + 1            
            T_hot.append(result['HE.gasFlow.heatTransfer.T[%s]' % idx].val)
            T_cold.append(result['HE.fluidFlow.heatTransfer.T[%s]' % idx].val)
            # pa -> kPa
            dp_hot_cur += result['HE.gasFlow.heatTransfer.dp[%s]' % idx].val / 1e3
            dp_cold_cur += result['HE.fluidFlow.heatTransfer.dp[%s]' % idx].val / 1e3
            dp_hot.append(dp_hot_cur)
            dp_cold.append(dp_cold_cur)
        
        T_cold.reverse()
        dp_cold.reverse()

        x_values = np.arange(0, len_seg * (N_seg) , len_seg) 

        plot_cfgs = {
            'T_hot': {
                'x' : x_values,
                'y' : np.array(T_hot),
                'range_x': axis_x,
                'range_y': axis_T,
                'cs': 'r-s',
                'label' : [label_x, label_T]
            },
            'T_cold': {
                'x' : x_values + len_seg,
                'y' : np.array(T_cold),
                'range_x': axis_x,
                'range_y': axis_T,
                'cs': 'b--s',
                'label' : [label_x, label_T]
            },
            'dp_hot': {
                'x' : x_values,
                'y' : np.array(dp_hot),
                'range_x': axis_x,
                'range_y': axis_dp,
                'cs': 'r-^',
                'label' : [label_x, label_dp],
                'ax_type': AxisType.Secondary
            },
            'dp_cold': {
                'x' : x_values + len_seg,
                'y' : np.array(dp_cold),
                'range_x': axis_x,
                'range_y': axis_dp,
                'cs': 'b--^',
                'label' : [label_x, label_dp],
                'ax_type': AxisType.Secondary
            }  
        }

        plot = PlotManager(title=cfg.name)

        for (k, v) in plot_cfgs.items():
            plot.add(DataSeries(name=k, **v))

        return plot  
    
    def load_result(self, exp_name) -> TestDataSet:
        '''
        Load simulation result from the latest, saved file
        '''
        import glob
        import xlwings as xw

        # find test_files 
        file_name = sep.join([self.path_out, exp_name+'.xlsx'])

        from pathlib import Path

        data_file = Path(file_name)

        if not data_file.exists():
            raise FileNotFoundError(f"data file {file_name} not found")

        app = xw.App(visible=False)
        ds_exp:TestDataSet = TestDataSet('empty test', {})

        try:
            wbk:xw.Book = xw.Book(file_name,read_only=True)
            col_start = TestConstants.DATA_FILE_COL_START

            for sht in wbk.sheets:
                if sht.name.startswith('Sheet'):
                    continue
                
                gname = sht.name

                helper = ExcelHelper(sht)
                u = 1
                offset = 2 # + 1 for linespacing, + 1 is start of next table 

                cfg_dict = helper.read_dict((u, col_start))
                u = u + len(cfg_dict) + offset
                ds_exp.set_cfgs(gname, cfg_dict)

                result_dict = helper.read_dict((u, col_start))
                u = u + len(result_dict) + offset
                ds_exp.set_results(gname, result_dict)   

        finally:
            app.quit()

        return ds_exp    

def gen_exp(exp_name:str ) -> TestDataSet:
# referred base cfg
    cfg_ref = {
        # geometry parameters
        "N_ch": 1e4, # "channel number"
        "N_seg": 10, # "segments number"
        "D_ch": 2e-3, # "channel diameter, semi circular tube"
        "L_fp": 200e-3, # "channel flow path length"
        "L_pitch": 12.3e-3, # "pitch length"
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
    
    
    ds_exp = TestDataSet.gen_batch_cfg(cfg_ref, cfg_offset, gname_template="Meshram_Test", ds_name=exp_name)
    
    # ds_test.add_view(mapping, view_name_tmpl="performance map")

    json_str = ds_exp.to_json()    
    print(json_str)

    return ds_exp

def gen_sol_dict() -> dict:
    ports = [
        ('HE.fluidFlow.heatTransfer', 'cold'),
        ('HE.fluidFlow.heatTransfer', 'cold'),
        ('HE.gasFlow.heatTransfer', 'hot'),
        ('HE.gasFlow.heatTransfer', 'hot')]

    # (propName, unit)
    props = [
        # ('fluid', '1'), # error in parsing this props, since returned name(string) are not digitial values
        ('T[%s]', 'K'), 
        ('dp[%s]', 'Pa'),
        ('Re[%s]', '[1]'),
        ('rho[%s]', 'kg/m3')]

    sol_dict = OrderedDict()

    for (p, p_short) in ports:
        for (prop, unit) in props:
            for n_seg in range(1, 11):
                prop_name = f'{prop}' % n_seg
                key = f'{p}.{prop_name}'
                text = f'{prop_name}_{p_short}'            
                sol_dict[key] = Variable(key, unit, text=text)   

    # specilized solutions
    var_sp = [               
        Variable('HE.waterIn.m_flow', 'kg/s'),
        Variable('HE.gasIn.m_flow', 'kg/s'),
        Variable('HE.gasIn.m_flow', 'kg/s'),
        Variable('HE.fluidFlow.heatTransfer.phi'),
        Variable('HE.fluidFlow.heatTransfer.pitch'),   
        Variable('HE.fluidFlow.heatTransfer.kim_cor_a.y'),
        Variable('HE.fluidFlow.heatTransfer.kim_cor_b.y'),
        Variable('HE.fluidFlow.heatTransfer.kim_cor_c.y'),
        Variable('HE.fluidFlow.heatTransfer.kim_cor_d.y'),
        Variable('HE.gasFlow.heatTransfer.phi'),
        Variable('HE.gasFlow.heatTransfer.pitch'),        
        Variable('HE.gasFlow.heatTransfer.kim_cor_a.y'),
        Variable('HE.gasFlow.heatTransfer.kim_cor_b.y'),
        Variable('HE.gasFlow.heatTransfer.kim_cor_c.y'),
        Variable('HE.gasFlow.heatTransfer.kim_cor_d.y')
    ]

    for v in var_sp:
        sol_dict[v.key] = v

    return sol_dict    

def main(work_root = []):
    # root path of modelica root
    if work_root == []:
        work_root = os.path.abspath(os.curdir)  

    run_sim = True # True: run the simulation, False: load pre-saved exp results
    model_name = "Steps.Test.TestTP_PCHE_Meshram"  

    exp:Experiments = MeshramTest(
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
        ds_exp = gen_exp(exp_name)
        exp.simulate(
            sim_ops=[
                'startTime=0', 
                'stopTime=10',
                'stepSize=2',
                'solver=dassl',
                '-nls=homotopy',
                '-lv=LOG_DEBUG,LOG_INIT,LOG_NLS,LOG_NLS_V,LOG_STATS'],
            solution_dict=gen_sol_dict(), 
            ds_test=ds_exp,
            append_save=False)   

        exp.save_results(ds_exp)
    else:
        # should search the exp_name 
        exp_name = "Test-Meshram 2020-12-06-12-28-34"
        ds_exp = exp.load_result(exp_name)

    plot_cfg_base = {
            "axis_x" : [0, 200e-3],
            "axis_T" : [400.0, 750.0],   
            "axis_dp" : [0, 80],
            "label_x" : 'Z (m)',
            "label_T" : 'T (K)',
            "label_dp" : '$\\Delta~p~(kPa)$ (K)'
        }

    plot_offest = {
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

    for group in ds_exp.values():
        for test in group.values():
            plot_cfg = deepcopy(plot_cfg_base)
            for (k, v) in plot_offest[test.name].items():
                plot_cfg[k] = v # replace with varied cfg

            srcfile = plot_cfg["imgfile"]
            dest_dir = sep.join([exp.path_out, exp_name])

            from pathlib import Path

            dir = Path(dest_dir)

            if not dir.exists():
                os.mkdir(dir)
            
            destfile = sep.join([dest_dir, "{0}_w_{1}.png".format(srcfile, test.name)])

            srcfile = sep.join([exp.path_pics, srcfile])
            exp.gen_plot_manager(test, plot_cfg).draw(img_file=srcfile, dest_file=destfile)

    print('All done!')

###
if __name__ == "__main__":
    main()