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

class MarchionniTest(Experiments):
    '''
    script for test against Marchionni [2019]
    '''

    def gen_plot_manager(self, test: TestItem, plot_cfg:dict) -> PlotManager:
        cfg = test.cfg
        result = test.result

        N_seg = int(cfg['N_seg'])
        L_fp = cfg['L_fp']
        len_seg = L_fp / N_seg   

        axis_x = plot_cfg["axis_x"]
        axis_T = plot_cfg["axis_T"]
        label_x = plot_cfg["label_x"]
        label_T =  plot_cfg["label_T"]

        T_hot = []
        T_cold = []

        for i in range(0, N_seg):
            idx = i + 1            
            T_hot.append(result['HE.gasFlow.heatTransfer.T[%s]' % idx].val)
            T_cold.append(result['HE.fluidFlow.heatTransfer.T[%s]' % idx].val)
        
        T_cold.reverse()

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
        "kc_dp": 1 # "pressure d"
    }

    # cfg with varied parameters from the base cfg
    cfg_offset = {} 
    cfg_offset['Group 1'] = {
            # values in Table 3 in Marchionni [2016]
            # Only with kc_dp = 1 at High T and kc_dp = 2 at Low T can I get good agreement with 
            # Marchionni's DP result. 
            "keys" : ["a_phi"],            
            "a_phi=30" : [30.0],
            "a_phi=35" : [35.0],
            "a_phi=40" : [40.0],
            "a_phi=45" : [45.0]
        }
    
    ds_exp = TestDataSet.gen_batch_cfg(cfg_ref, cfg_offset, gname_template="Marchionni_Test", ds_name=exp_name)
    
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
            for n_seg in range(1, 21):
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
    model_name = "Steps.Test.TestTP_PCHE_Marchionni"  

    exp:Experiments = MarchionniTest(
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
        exp_name = 'Test-Marchionni {:%Y-%m-%d-%H-%M-%S}'.format(datetime.datetime.now())            
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
        exp_name = "Test-Marchionni 2020-12-06-12-28-34"
        ds_exp = exp.load_result(exp_name)

    plot_cfg_base = {
            "axis_x" : [0, 272e-3],
            "axis_T" : [from_degC(80), from_degC(440)],   
            "label_x" : 'Z (m)',
            "label_T" : 'T (K)',
            "imgfile" : "Marchionni_Fig_04_T_1D_vs_3D.png"
        }

    plot_offest = {}    

    for group in ds_exp.values():
        for test in group.values():
            plot_cfg = deepcopy(plot_cfg_base)
            if test.name in plot_offest.keys():
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