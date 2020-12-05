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

from model import TestDataSet,TestItem,Variable, TestConstants
from plotlib import PlotManager, DataSeries, AxisType
from physics import Temperature, Pressure, MDot
from experiments import Experiments
from utils import ExcelHelper, from_degC, from_bar

class MeshramTest(Experiments):
    '''
    script for test against Meshram [2016]
    '''

    def gen_plot_manager(self, test: TestItem) -> PlotManager:
        cfg = test.cfg
        N_seg = cfg['N_seg']
        L_fp = cfg['L_fp']
        len_seg = L_fp / N_seg

        result = test.result
          
        # axis_x=[[0, 0.12], [0, 0.16]][zigzag]
        axis_x=[0, L_fp]
        axis_T=[400.0, cfg['T_hot_in']]
        axis_dp=[0, 80]

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

        label_x = 'Z (m)'
        label_T = 'T (K)'
        label_dp = '$\\Delta~p~(kPa)$ (K)'

        plot = PlotManager(title=cfg.name)
        plot.add(DataSeries(name = 'T_hot', x = x_values, y = np.array(T_hot), range_x=axis_x, range_y=axis_T, cs = 'r-s', label=[label_x, label_T]))
        plot.add(DataSeries(name = 'T_cold', x = x_values + len_seg, y = np.array(T_cold), range_x=axis_x, range_y=axis_T, cs = 'b--s', label=[label_x, label_T]))
        plot.add(DataSeries(name = 'dp_hot', x = x_values, y = np.array(dp_hot), range_x=axis_x, range_y=axis_dp, cs = 'r-^', label=[label_x, label_dp]), ax_type=AxisType.Secondary)
        plot.add(DataSeries(name = 'dp_cold', x = x_values + len_seg, y = np.array(dp_cold), range_x=axis_x, range_y=axis_dp, cs = 'b--^', label=[label_x, label_dp]), ax_type=AxisType.Secondary)

        return plot  


def main(work_root = []):
    # root path of modelica root
    if work_root == []:
        work_root = os.path.abspath(os.curdir)  

    model_name = "Steps.Test.TestTP_PCHE_Meshram"     

    # referred base cfg
    cfg_ref = {
        # geometry parameters
        "N_ch": 1e4, # "channel number"
        "N_seg": 10, # "segments number"
        "D_ch": 2e-3, # "channel diameter, semi circular tube"
        "L_fp": 200e-3, # "channel flow path length"
        "pitch": 12.3e-3, # "pitch length"
        "phi": 36, # "pitch angle °"
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
    }

    # cfg with varied parameters from the base cfg
    cfg_offset = {} 
    mapping = []
    cfg_offset['Group 1'] = {
            # values in Table 3 in Meshram [2016]
            "keys" : ["T_cold_in", "T_hot_in", "u_cold_in", "u_hot_in"],
            "zigzag High T" : [500, 730, 1.876, 7.564],
            "zigzag Low T": [400, 630, 0.806, 1.345],
            "Straight High T": [500, 730, 1.518, 6.118],
            "Straigth Low T": [400, 630, 0.842, 4.702]
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
    
    ds_name = 'Test-Meshram {:%Y-%m-%d-%H-%M-%S}'.format(datetime.datetime.now())
    ds_test = TestDataSet.gen_batch_cfg(cfg_ref, cfg_offset, gname_template="Meshram_Test", ds_name=ds_name)
    
    # ds_test.add_view(mapping, view_name_tmpl="performance map")

    json_str = ds_test.to_json()    
    print(json_str)

    ports = [
        'HE.fluidFlow.heatTransfer',
        'HE.fluidFlow.heatTransfer',
        'HE.gasFlow.heatTransfer',
        'HE.gasFlow.heatTransfer']

    # (propName, unit)
    props = [
        # ('fluid', '1'), # error in parsing this props, since returned name(string) are not digitial values
        ('T[%s]', 'K'), 
        ('dp[%s]', 'Pa'),
        ('Re[%s]', 'kg/m3')]

    # props = [
    #     ('T','K'), 
    #     ('p', 'Pa'),
    #     ('h', 'kJ/kg'),
    #     ('w', 'kg/s')]

    sol_dict = OrderedDict()

    for p in ports:
        for (prop, unit) in props:
            for n_seg in range(1, 11):
                key = f'{p}.{prop}' % n_seg
                text = f'{prop}_{p[1:]}'            
                sol_dict[key] = Variable(key, unit, text=text)    

    # specilized solutions
    # var_sp = [
    #     Variable('W_net', 'MW'),
    #     Variable('Q_heater', 'MW'),
    #     Variable('eta_pb', '%'),
    #     Variable('SR', '%'),
    #     Variable('W_turb', 'MW'),
    #     Variable('eta_turb', '%'),
    #     Variable('W_MC', 'MW'),
    #     Variable('eta_MC', '%'),
    #     Variable('W_RC', 'MW'),
    #     Variable('eta_RC', '%'),
    #     Variable('Q_HTR', 'MW'),
    #     Variable('dT1_HTR', 'K'),
    #     Variable('dT2_HTR', 'K'),
    #     Variable('T_ltmd_HTR', 'K'),
    #     Variable('UA_HTR', 'MW/(S^2 K)'),
    #     Variable('Q_LTR', 'MW'),
    #     Variable('dT1_LTR', 'K'),
    #     Variable('dT2_LTR', 'K'),
    #     Variable('T_ltmd_LTR', 'K'),
    #     Variable('UA_LTR', 'MW/(S^2 K)'),
    #     Variable('T_heater_hot_out', 'K')        
    # ]

    # for var in var_sp:
    #     sol_dict[var.key] = var

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

    exp.simulate(
        sim_ops=[
            'startTime=0', 
            'stopTime=10',
            'stepSize=2',
            'solver=dassl',
            '-nls=homotopy',
            '-lv=LOG_DEBUG,LOG_INIT,LOG_NLS,LOG_NLS_V,LOG_STATS'],
        solution_dict=sol_dict, 
        ds_test=ds_test,
        append_save=False)   

    exp.save_results(ds_test)

    zigzag = 1
    imgfile = ["Meshram_Fig_04.png", "Meshram_Fig_05.png"][zigzag]    
    imgfile = sep.join([exp.path_pics, imgfile])

    for group in ds_test.values():
        for test in group.values():
            destfile = sep.join([exp.path_out, "Meshram_Fig{0}b_compare_{1}.png".format(4 + zigzag, test.name)])
            exp.gen_plot_manager(test).draw(img_file=imgfile, dest_file=destfile)

    print('All done!')

###
if __name__ == "__main__":
    main()