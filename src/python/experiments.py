from collections import OrderedDict
from logging import exception
import datetime
from math import trunc
from copy import deepcopy
from os.path import join, sep
import shutil
import sys
from OMPython import OMCSessionZMQ
from OMPython import ModelicaSystem
from datetime import datetime as dt
from numpy.core import test
import pandas as pd
import xlwings as xw
from typing import Tuple, List
from numpy.core.fromnumeric import mean

from os import path,sep
import os
import numpy as np
from enum import Enum

# from ex.CoolPropQuery import result

from model import TestDataSet, TestGroup, Variable, TestResult, TestItem, TestConstants
from CoolProp.CoolProp import PhaseSI, PropsSI, get_global_param_string
from plotlib import PlotManager, DataSeries, AxisType
from physics import Temperature, Pressure, MDot
from utils import ExcelHelper, mkdir_filepath

class Experiment(object):
        
    def __init__(self, work_root, model_name, ex_dlls=[], modelica_libs=[]):
        '''
            ex_dlls : external pre-compiled dlls such as cool prop
            modelica_libs: 3rd party modelica libraries like ThermoPower, ExternalMedia
        '''
        super().__init__()
        # common variable defination
        self.path_root = work_root
        from pathlib import Path
        parent_dir = Path(self.path_root).parent
        self.path_model = sep.join([parent_dir.__str__(), "Modelica", "Steps"])
        self.path_workspace = sep.join([self.path_root, "build"])
        self.path_pics = sep.join([self.path_workspace, "pics"])
        self.path_out = sep.join([self.path_workspace, "out"])
        self.model_name = model_name
        self.modelica_libs = modelica_libs
        self.ex_dlls = ex_dlls

        self.model_loaded = False

        self.prepare_workspace()       

    def gen_sol_dict(self, cfg:dict) -> dict:
        '''
            generate the solution dict for a experiment
        '''
        return {}

    def prepare_workspace(self):
        '''
        prepare the workspace for simulation, the directory structure is
        - Steps (currernt work directory) -> path_root
        |--.vscode              - directory for vscode's configurations
        |--docs                 - documents
        |--lib                  - compiled libraries
        |--build                - -> path_workspace, workspace for modelica simulation, all
        |   |                        the relative paths in this program are relative to point
        |   |--out               - Final output of the Simulation
        |   |  |--Test_Batch 1   - Simulation output data for test 1
        |   |  ...               - Simulation output data for test ...
        |   |--pics              - Temp storage dir for referenced pic 
        |--src
        |   |--c                 - c++ sources
            |--Modelica          - Modelica models
                |--Steps                - Directory for modelica models codes
            |--scripts           - shell scripts for compiling and tasks
            |--Python            - Python scripts for parameters sweep          
        |.env                   - configuration file for vscode python extension
        |.gitignore             - git ignore file
        ''' 

        import os
        import shutil      

        # setup working directory
        if not os.path.exists(self.path_workspace):    
            os.mkdir(self.path_workspace)

        os.chdir(self.path_workspace)

        if os.path.exists(self.path_pics):
            shutil.rmtree(self.path_pics)     

        # copy referenced pics for output figure 
        shutil.copytree(sep.join([self.path_root, 'pics']), self.path_pics)

        # directory for output        
        if not os.path.exists(self.path_out):
            os.mkdir(self.path_out)        

        # copy dlls such as 'cool prop.dll or libexternalmedia.dll'
        lib_dlls = self.ex_dlls
        for lib_dll in lib_dlls: 
            if os.path.exists(lib_dll):  
                dll_name = str.split(lib_dll,'/')[-1]
                try:
                    shutil.copyfile(lib_dll, sep.join([".", dll_name])) # completa target name needed 
                except:
                    print("Unexpected error in copying dlls:", sys.exc_info()[0])
                    pass
                        
            else:              
                print(f'Dll {lib_dll} not Exists')

    def load_model(self):

        if not self.model_loaded:

            print(f'loading model{self.model_name}...\n')        
            self.model = ModelicaSystem(sep.join([self.path_model,"package.mo"]), self.model_name,self.modelica_libs) 
            print(f'model{self.model_name} loaded, prepare to simulate\n')

            self.model_loaded = True

    def simulate_batch(self, sim_ops, solution_dict, ds_test:TestDataSet, append_save=True):         

        for group in ds_test.values():
            for (test_name, test) in group.items():

                print(f' ready to run test={test_name} in group={group.name}\n')

                self.simulate(sim_ops=sim_ops, solution_dict=solution_dict, test=test, data_ref=ds_test.ref_data)               

                if(append_save): # append the result each time we finish the simulation
                    # use a copy for saving. avoid concurrent iteration on same object
                    self.save_results(deepcopy(ds_test))
        
        print('all simulation(s) done, results saved in Test Data Set')        

    def simulate(self, sim_ops, solution_dict, test:TestItem, data_ref=None):

        sol_keys = list(solution_dict.keys())
        solutions = np.ones((len(sol_keys),1)) * -1 # default value        

        try:

            if not self.model_loaded:
                self.load_model() # load model once and once only to save time for batch execution

            self.model.setParameters(test.cfg.gen_params())
            self.model.setSimulationOptions(sim_ops)
            self.model.simulate()
            solutions = self.model.getSolutions(sol_keys) 
        except:
            pass

        test.result = TestResult.from_keyvalues(solution_dict, solutions) 

        self.post_process(test, data_ref)

    def post_process(self, test:TestItem, data_ref):
        pass

    def save_results(self, ds_test : TestDataSet):
        '''
        save the simulation result into files   
        override save_results() to save all result in one excel file         
        '''
        # create directory for current test batch   
        from pathlib import Path
        dir_name = path.join(self.path_out, ds_test.name)
        dir = Path(dir_name)

        if not dir.exists():
            os.mkdir(dir)

        filename = path.join(dir_name, ds_test.name + '.xlsx')

        xw_app: xw.App = xw.App(visible=False)       

        try:
            wbk = xw_app.books.add() 
                    
            for (gname, group) in ds_test.items():     
                cfg_set =  ds_test.get_cfgs(gname)  
                result_set = ds_test.get_results(gname)
                view_set = ds_test.get_views(gname)
                # a sheet for each group

                sht = wbk.sheets.add(name = group.name[0:30]) # 30 - max lenth for a sheet name
                ex: ExcelHelper = ExcelHelper(sht)
                u, l = 1, TestConstants.DATA_FILE_COL_START
                title = {"Table 1": "Test Config"}
                # config
                (b, r) = ex.write_table(cfg_set, title=title, up=u, left=l, linespacing=True)

                idx = 2
                if view_set: # not empty
                    for (view_name, data_table) in view_set.items():                    
                        title = {f"Table {idx}": view_name}
                        (b, r) = ex.write_table(data_table, title=title, up=b, left=l, linespacing=True)
                        idx += 1

                title = {f"Table {idx}": "Detailed Results"}
                (b, r) = ex.write_table(result_set, title=title, up=b, left=l, linespacing=True)

            from pathlib import Path

            ex_file=Path(filename)       
            if ex_file.exists():
                shutil.move(filename, filename+".bak")

            wbk.save(filename)                    

        finally:
            xw_app.quit() # close the xlwings app finally 

    def load_results(self, exp_name, dir_name=None, ds_exp:TestDataSet = None, has_view=False) -> TestDataSet:
        '''
        Load simulation result from the latest, saved file
        '''
        import glob
        import xlwings as xw

        # find test_files 
        if dir_name != None:
            file_name = sep.join([self.path_out, dir_name, exp_name+'.xlsx'])    
        else:
            file_name = sep.join([self.path_out, exp_name+'.xlsx'])

        from pathlib import Path

        data_file = Path(file_name)

        if not data_file.exists():
            raise FileNotFoundError(f"data file {file_name} not found")

        app = xw.App(visible=False)
        
        if(ds_exp == None):
            ds_exp = TestDataSet(exp_name, {})

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

                # read config dict
                cfg_dict = helper.read_dict((u, col_start))
                u = u + len(cfg_dict) + offset
                ds_exp.set_cfgs(gname, cfg_dict)

                if has_view:
                    # read & neglect summary result, set detailed result instead
                    cfg_dict = helper.read_dict((u, col_start))
                    u = u + len(cfg_dict) + offset          

                # read detailed result
                result_dict = helper.read_dict((u, col_start))
                u = u + len(result_dict) + offset
                ds_exp.set_results(gname, result_dict)   

            for g in ds_exp.values():
                for t in g.values():
                    self.post_process(t,ds_exp)

        finally:
            app.quit()

        return ds_exp  

    def plot_results(self, ds_exp:TestDataSet, plot_metacfg, plot_metacfg_offset):

        for g in ds_exp.values():
            for test in g.values():
                cfg = test.cfg
                values = test.post_data

                plot_metacfg_cp = deepcopy(plot_metacfg)

                if test.name in plot_metacfg_offset.keys():
                    for (k, v) in plot_metacfg_offset[test.name].items():
                        plot_metacfg_cp[k] = v # replace with varied cfg

                plot_cfgs = self.gen_plot_cfgs(values, plot_metacfg_cp)       

                for (k1, plot_cfg) in plot_cfgs.items():

                    plot = PlotManager(title=cfg.name)

                    for series in plot_cfg['series']:
                        plot.add(DataSeries(**series))

                    src_file = sep.join([self.path_pics, plot_cfg["imgfile"]])
                    dest_file = sep.join([self.path_out, ds_exp.name, g.name, "{0}_w_t={1}.png".format(plot_cfg["imgfile"], test.name)])

                    plot.draw(src_file, dest_file)   

    def gen_plot_cfgs(self, values, meta_cfg) -> dict:
        '''
            use template to generate all plot configs 
        '''
        return {}

    def run(self,  simulate = True):
        pass    

    def find_latest_test(self, path_) -> str:
        return ""        
class PCHEExperiment(Experiment):

    def post_process(self, test:TestItem, data_ref):
        '''
            post process the test result after simulation
            data mapping and re-collection can be applied in this function
        '''
 
        cfg = test.cfg
        result = test.result

        N_seg = int(cfg['N_seg'])
        N_nodes = N_seg + 1
        L_fp = cfg['L_fp']
        len_seg = L_fp / N_seg

        T_hot = []
        T_cold = []
        T_vol_hot = []
        T_vol_cold = []
        gamma_hot = []
        gamma_cold = []
        hc_hot = []
        hc_cold = []
        dp_hot = [0] 
        dp_cold = [0]
        # started, fix pressure of each stream
        dp_hot_cur = 0
        dp_cold_cur = 0   

        Re_hot = []
        Re_cold = []
        Pr_hot = []
        Pr_cold = []

        for i in range(0, N_seg):
            idx = i + 1            
            T_hot.append(result['HE.gasFlow.heatTransfer.T[%s]' % idx].val)
            T_cold.append(result['HE.fluidFlow.heatTransfer.T[%s]' % idx].val)

            T_vol_hot.append(result['HE.gasFlow.heatTransfer.T_vol[%s]' % idx].val)
            T_vol_cold.append(result['HE.fluidFlow.heatTransfer.T_vol[%s]' % idx].val) 

            # pa -> kPa
            dp_hot_cur += result['HE.gasFlow.heatTransfer.dp[%s]' % idx].val / 1e3
            dp_cold_cur += result['HE.fluidFlow.heatTransfer.dp[%s]' % idx].val / 1e3
            # accumulated pressure drop, starts from 0
            dp_hot.append(dp_hot_cur)
            dp_cold.append(dp_cold_cur)

            hc_hot.append(result['HE.gasFlow.heatTransfer.hc[%s]' % idx].val)
            hc_cold.append(result['HE.fluidFlow.heatTransfer.hc[%s]' % idx].val)

            # if i <= N_seg - 2: # len(gamma) = len(T) - 1
            gamma_hot.append(result['HE.gasFlow.heatTransfer.gamma[%s]' % idx].val)
            gamma_cold.append(result['HE.fluidFlow.heatTransfer.gamma[%s]' % idx].val)

            Re_hot.append(result['HE.gasFlow.heatTransfer.Re[%s]' % idx].val)
            Re_cold.append(result['HE.fluidFlow.heatTransfer.Re[%s]' % idx].val)
            Pr_hot.append(result['HE.gasFlow.heatTransfer.Pr[%s]' % idx].val)
            Pr_cold.append(result['HE.fluidFlow.heatTransfer.Pr[%s]' % idx].val)           
        
        # last node
        T_hot.append(result['HE.gasFlow.heatTransfer.T[%s]' % N_nodes].val)
        T_cold.append(result['HE.fluidFlow.heatTransfer.T[%s]' % N_nodes].val)    
        # dp_hot.append(dp_hot_cur)
        # dp_cold.append(dp_cold_cur)    

        T_cold.reverse()
        T_vol_cold.reverse()
        dp_cold.reverse()  
        hc_cold.reverse()  
        gamma_cold.reverse()        
        
        x_values = np.arange(0, len_seg * (N_nodes) , len_seg) 
        # x_values_gamma = np.arange(0, len_seg * (N_seg - 1) , len_seg) 

        test.set_post_data('x',x_values)
        test.set_post_data('x_hot', x_values[:N_seg]) # N_seg = N_nodes - 1
        test.set_post_data('T_hot',np.array(T_hot))
        test.set_post_data('dp_hot',np.array(dp_hot))
        test.set_post_data('x_cold', x_values[1:]) # 1:N_nodes
        test.set_post_data('T_cold', np.array(T_cold))
        test.set_post_data('dp_cold',np.array(dp_cold))
        # test.set_post_data('x_values_gamma',np.array(x_values_gamma))  
        test.set_post_data('hc_hot',np.array(hc_hot))
        test.set_post_data('hc_cold',np.array(hc_cold))      
        test.set_post_data('gamma_hot',np.array(gamma_hot))
        test.set_post_data('gamma_cold',np.array(gamma_cold))

        results_cal = []

        # actual averaged Re and Pr numbers
        results_cal.append(('Re_bar_hot', '[1]', sum(Re_hot) / len(Re_hot)))
        results_cal.append(('Re_bar_cold', '[1]', sum(Re_cold) / len(Re_cold)))

        results_cal.append(('Pr_bar_hot', '[1]', sum(Pr_hot) / len(Pr_hot)))
        results_cal.append(('Pr_bar_cold', '[1]', sum(Pr_cold) / len(Pr_cold)))

        # calculate averaged values given boundury conditions
        T_hot_bar = (cfg['T_hot_in'] + cfg['T_hot_out']) / 2
        T_cold_bar = (cfg['T_cold_in'] + cfg['T_cold_out']) / 2
        p_hot_bar = cfg['p_hot_in']
        p_cold_bar = cfg['p_cold_in']

        from CoolProp.CoolProp import PhaseSI, PropsSI
        import math
        
        medium = "CO2"
        mu_bar_hot = PropsSI("VISCOSITY", "P", p_hot_bar, "T", T_hot_bar, medium)
        mu_bar_cold = PropsSI("VISCOSITY", "P", p_cold_bar, "T", T_cold_bar, medium)
        cp_bar_hot = PropsSI("CPMASS", "P", p_hot_bar, "T", T_hot_bar, medium)
        cp_bar_cold = PropsSI("CPMASS", "P", p_cold_bar, "T", T_cold_bar, medium)
        k_bar_hot = PropsSI("CONDUCTIVITY", "P", p_hot_bar, "T", T_hot_bar, medium)
        k_bar_cold = PropsSI("CONDUCTIVITY", "P", p_cold_bar, "T", T_cold_bar, medium)

        results_cal.append(('mu_bar_hot', '[Pa S]', mu_bar_hot))
        results_cal.append(('mu_bar_cold', '[Pa S]', mu_bar_cold))

        results_cal.append(('cp_bar_hot', '[J/m^3-K]', cp_bar_hot))
        results_cal.append(('cp_bar_cold', '[J/m^3-K]', cp_bar_cold))

        results_cal.append(('k_bar_hot', '[J/m^3-K]', k_bar_hot))
        results_cal.append(('k_bar_cold', '[J/m^3-K]', k_bar_cold))

        # hydraulic diameter
        d_h = cfg['D_ch'] * math.pi / ( math.pi + 2)

        # calculated channel-wise properties

        results_cal.append(('Re_tilde_hot', '[1]', result['HE.gasFlow.heatTransfer.G[1]'].val * d_h / mu_bar_hot))
        results_cal.append(('Re_tilde_cold', '[1]', result['HE.fluidFlow.heatTransfer.G[1]'].val * d_h / mu_bar_cold))

        results_cal.append(('Pr_tilde_hot', '[1]', cp_bar_hot * mu_bar_hot / k_bar_hot))
        results_cal.append(('Pr_tilde_cold', '[1]', cp_bar_cold * mu_bar_cold / k_bar_cold))

        for (k, u, v) in results_cal:
            test.result[k] = Variable(k,u,v)


    def gen_sol_dict(self, cfg:dict) -> dict:
        ports = [
            ('HE.fluidFlow.heatTransfer', 'cold'),
            ('HE.fluidFlow.heatTransfer', 'cold'),
            ('HE.gasFlow.heatTransfer', 'hot'),
            ('HE.gasFlow.heatTransfer', 'hot')]

        # (propName, unit)
        N_seg =  cfg['N_seg'] # volumes number
        N_nodes = N_seg + 1 # nodes number
        props = [
            # ('fluid', '1'), # error in parsing this props, since returned name(string) are not digitial values
            ('T[%s]', 'K', N_nodes), 
            ('T_vol[%s]', 'K', N_seg), 
            ('dp[%s]', 'Pa', N_seg),
            ('Re[%s]', '[1]', N_seg),
            ('Pr[%s]', '[1]', N_seg),
            ('rho[%s]', 'kg/m3', N_seg),
            ('hc[%s]', 'W/m2-K', N_seg),
            ('gamma[%s]', 'W/m2-K', N_seg),
            ('G[%s]', 'kg/m^2-s', N_seg)]

        sol_dict = OrderedDict()

        for (p, p_short) in ports:
            for (prop, unit, n) in props:
                for idx in range(1, n + 1):
                    prop_name = f'{prop}' % idx
                    key = f'{p}.{prop_name}'
                    text = f'{prop_name}_{p_short}'            
                    sol_dict[key] = Variable(key, unit, text=text)   

        # specilized solutions
        var_sp = [               
            Variable('HE.waterIn.m_flow', 'kg/s'),
            Variable('HE.gasIn.m_flow', 'kg/s'),
            Variable('HE.fluidFlow.heatTransfer.phi'),
            Variable('HE.fluidFlow.heatTransfer.pitch'),   
            Variable('HE.gasFlow.heatTransfer.phi'),
            Variable('HE.gasFlow.heatTransfer.pitch') 
        ]

        for v in var_sp:
            sol_dict[v.key] = v

        return sol_dict    

    def Nu_cor(self, text, T_b, T_w, p_b, p_w, G, d, medium, verbose=False):
        g = 9.80665 # gravity acceleration

        rho_w = PropsSI("D", "T", T_w, "P", p_w, medium)
        rho_b = PropsSI("D", "T", T_b, "P", p_b, medium)
        mu_b = PropsSI("VISCOSITY", "T", T_b, "P", p_b, medium)
        cp_b = PropsSI("CPMASS", "T", T_b, "P", p_b, medium)
        k_b = PropsSI("CONDUCTIVITY", "T", T_b, "P", p_b, medium)
        h_w = PropsSI("HMASS", "T", T_w, "P", p_w, medium)
        h_b = PropsSI("HMASS", "T", T_b, "P", p_b, medium)
        cp_bar = (h_w - h_b) / (T_w - T_b)

        Gr = abs((rho_w - rho_b) * rho_b * g * (d**3) / (mu_b ** 2))
        Re = G * d / mu_b
        Pr = cp_b * mu_b / k_b

        q1 = Gr / (Re ** 2)
        q2 = rho_w / rho_b
        q3 = cp_bar / cp_b

        Nu = 0.124 * (Re ** 0.8) * (Pr ** 0.4) * (q1 ** 0.203) * (q2 ** 0.842) * (q3 ** 0.284)
        # if q2 > 1:
        #     q2 = 1 / q2
        if verbose:
            print(f'{text}: Nu={Nu}, q1={q1}, q2={q2}, q3={q3}')

        return (Nu, q1, q2, q3)

    # bcs = [bc1, bc2, ...]
    # bcn = (st_hot_in, st_cold_in)
    # st = (T [K], p [bar->ba])

    def predict_cor_coff(self, T_h, p_h, T_c, p_c, G, d, medium, verbose=False) -> Tuple[float, float]:
        '''
            use boundary conditions to predict C1
        '''
        sides ={
            "hot": {
                "text": "for hot side",
                "T_b": T_h,
                "T_w": (T_h + T_c)/2,
                "p_b": p_h, 
                "p_w": p_h,
                "G": G,
                "d": d,
                "medium": medium,
                "verbose": verbose,
            },
            "cold": {
                "text": "for cold side",
                "T_b": T_c,
                "T_w": (T_h + T_c)/2,
                "p_b": p_c, 
                "p_w": p_c,
                "G": G,
                "d": d,
                "medium": medium,
                "verbose": verbose   
            }
        }

        (Nu, q1, q2, q3) = ([], [], [], [])

        for side in sides.values():
            (v1, v2, v3, v4) = self.Nu_cor(**side)
            Nu.append(v1)
            q1.append(v2)
            q2.append(v3)
            q3.append(v4)   

        if verbose:
            print(f"Averaged value: Nu_bar={mean(Nu)}, q1_bar={mean(q1)}, q2_bar={mean(q2)}, q3_bar={mean(q3)}")
        C1 = 1.503212939 * mean(q3) -1.18150071
        C2 = -1.618997638 * mean(q3) + 1.892017384

        return (C1, C2)
class FittingStage(Enum):
    train = 1,
    test = 2

class ErrorFunc(Enum):
    T = 1, # consider temperature only 
    Dp = 2, # consider pressure drop only
    T_Dp_both = 3

class ParamFitting(object):
    
    KEY_GNAME = "C1_fitting"

    def __init__(self, exp, cases, mapping, sim_ops, errfunc = ErrorFunc.T) -> None:
        super().__init__() 
        self.exp:Experiment = exp
        self.cases = cases
        self.mapping = mapping
        self.sim_ops = sim_ops
        self.errfunc = errfunc
        self.case_name = "" # current case 


    def run_fitting(self, pt = (1,1,1), max_steps=50, num_sumples=4, steplength=0.1):
        '''
            execute the fitting to find optimized C1
        '''
        self.stage = FittingStage.train.name

        (h_pt, h_eval) = self.param_random_local_search(
                            func_para=self.evaluate, 
                            # Cf_a,b,c_hot, Cf_a,b,c_cold
                            pt=pt, 
                            max_steps=max_steps,
                            num_samples=num_sumples, 
                            steplength=steplength)

        print('fitting done')

    def param_random_local_search(self, func_para,pt,max_steps,num_samples,steplength):
        '''
            Concurrently train C1 for hot/cold side due to simulation is a time consuming process
        '''
        # starting point evaluation
        (current_eval, tests) = func_para(pt)
        dim_x = len(pt)
        dim_y = len(current_eval)
        current_pt = pt

        exp_name = 'param_fitting_{:%Y-%m-%d-%H-%M-%S}'.format(datetime.datetime.now()) 
        # test dataset keep trace of fitting process
        ds_exp = TestDataSet.new_test_dataset(self.cases['cfg'], 'EmptyGroup', exp_name)
        ds_exp.add_view(self.mapping)

        for test in tests:
            gname = test.name
            ds_exp[gname] = TestGroup(gname, {'{:%H-%M-%S}'.format(datetime.datetime.now()): test})

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

                # set x, y independently
                d_C1 = rnd.uniform(-steplength_temp, steplength_temp) 
                d_C2 = rnd.uniform(-steplength_temp, steplength_temp) 
                d_C3 = rnd.uniform(-steplength_temp, steplength_temp) 

                new_pt = np.asarray([d_C1, d_C2, d_C3])
                temp_pt = deepcopy(keeper_pt) 
                new_pt += temp_pt
                
                # evaluate new point
                (new_eval, tests) = func_para(new_pt)

                y_old = [current_eval[x] for x in range(0, dim_y)]
                y_new = [new_eval[x] for x in range(0, dim_y)]

                # evaluate results respectively
                update = 0
                for k in range(0, dim_y):
                    if y_new[k] < y_old[k]:
                        update += 1 # increase update counter if find a better y in this dimension

                if update == dim_y: # find a better value for all dimensions out of the samples
                    current_pt = [new_pt[x] for x in range(0, dim_x)]
                    current_eval = [y_new[x] for x in range(0, dim_y)]                    

                    swap = 1
            
            # if swap happened
            if swap == 1:
                pt_history.append(current_pt)
                eval_history.append(current_eval)   
                
                for test in tests:
                    gname = test.name
                    # change test name to instance
                    test.name = '{:%H-%M-%S}'.format(datetime.datetime.now())
                    # merge it with previous result, keep trace of fitting
                    ds_exp[gname].append_test(test)
            
        self.exp.save_results(ds_exp) # save 'better' result and its parameters for furture training

        return pt_history,eval_history   

    def evaluate(self, x_new) -> Tuple[List[float], List[TestItem]]:
        init = False
        err_sum = []
        tests = []

        for case_name, data_ref in self.cases["data_ref"].items():

            cfg = deepcopy(self.get_cfg(case_name))

            # add to-be-searched coefficients 
            keys = ['Cf_C1', 'Cf_C2', 'Cf_C3']
            cfg.update(
                {
                    keys[i]:x_new[i]
                    for i in range(0, len(x_new))
                })          

            test = TestItem(case_name, cfg, {})

            self.exp.simulate(
                sim_ops=self.sim_ops,
                
                solution_dict=self.exp.gen_sol_dict(cfg), 
                test=test)

            err_cur = self.cal_func_err(test.post_data, data_ref)

            # test.name = '{:%H-%M-%S}'.format(datetime.datetime.now())

            for i in range(0, len(err_cur)):
                k = f'err_fun_{i}'
                test.result[k] = Variable(k, val = err_cur[i])

            if not init:
                err_sum = err_cur
                init = True
            else:
                for i in range(0, len(err_cur)):
                    err_sum[i] += err_cur[i]

            tests.append(test)
            
            # merge it with previous result
            # g.append_test(test) 

        return (err_sum, tests)

    def prepare_data(self, values):
        return values

    def cal_func_err(self, values, data_ref):
        '''
            calculate function errors against data referred
        '''
        # pre process on the values 
        values = self.prepare_data(values)

        T_hs = values['T_hot']
        dp_hs = values['dp_hot']
        T_cs = values['T_cold']
        dp_cs = values['dp_cold']        

        err_h_dp = 0.0
        err_c_dp = 0.0
        err_h_T = 0.0
        err_c_T = 0.0
        errFunc = self.errfunc

        if errFunc == ErrorFunc.Dp:
            err_h_dp = abs(dp_hs[-1] - data_ref['dp_hs'][-1]) 
            # cold side in reverse order - 0 is the last one 
            err_c_dp = abs(dp_cs[0] - data_ref['dp_cs'][0])
        elif errFunc == ErrorFunc.T_Dp_both:
            err_h_dp = self.cal_err(dp_hs[-1], data_ref['dp_hs'][-1])
            err_c_dp = self.cal_err(dp_cs[0], data_ref['dp_cs'][0])

        if errFunc == ErrorFunc.T:
            err_h_T = sum([abs(T_hs[i] - data_ref['T_hs'][i]) for i in range(0, len(T_hs))])
            err_c_T = sum([abs(T_cs[i] - data_ref['T_cs'][i]) for i in range(0, len(T_cs))])
        elif errFunc == ErrorFunc.T_Dp_both:
            err_h_T = sum([self.cal_err(T_hs[i], data_ref['T_hs'][i]) for i in range(0, len(T_hs))])
            err_c_T = sum([self.cal_err(T_cs[i], data_ref['T_cs'][i]) for i in range(0, len(T_cs))])

        return [err_h_T + err_h_dp + err_c_T + err_c_dp]

    def get_cfg(self, case_name) -> dict:
        '''
            return config
        '''
        cfg_table = self.cases['cfg']
        if not case_name in cfg_table.keys():
            raise KeyError(f"No case with name={case_name}")
        
        return {
            cfg_table['keys'][i] : cfg_table[case_name][i]
            for i in range(0, len(cfg_table['keys']))
        }

    def get_data_ref(self, case_name):
        '''
            return data referenced
        '''

        if case_name in self.cases["data_ref"].keys():
            raise KeyError(f"No case with name={case_name}")
        
        return self.cases["data_ref"][case_name]

    def cal_err(self, x, y) -> float:
        if abs(y) < 1e-10:
            return 0
        
        return ((x - y) / y) ** 2   