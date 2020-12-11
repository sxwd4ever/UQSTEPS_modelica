from collections import OrderedDict
from logging import exception
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

from os import path,sep
import os
import numpy as np

# from ex.CoolPropQuery import result

from model import TestDataSet, Variable, TestResult, TestItem, TestConstants
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

    def simulate(self, sim_ops, solution_dict, ds_test:TestDataSet, append_save=True):        
        
        print(f'loading model{self.model_name}...\n')        
        self.model = ModelicaSystem(sep.join([self.path_model,"package.mo"]), self.model_name,self.modelica_libs) 
        print(f'model{self.model_name}loaded, prepare to simulate\n')

        # options for simulation - steady state simulation, no iteration required, so set numberOfIntervals = 2 

        for group in ds_test.values():
            for (test_name, test) in group.items():
                sol_keys = list(solution_dict.keys())
                solutions = np.ones((len(sol_keys),1)) * -1 # default value
                
                try:
                    print(f' ready to run test={test_name} in group={group.name}\n')
                    self.model.setParameters(test.cfg.gen_params())
                    self.model.setSimulationOptions(sim_ops)
                    self.model.simulate()
                    solutions = self.model.getSolutions(sol_keys) 
                except:
                    pass

                test.result = TestResult.from_keyvalues(solution_dict, solutions) 
                
                self.post_process(test, ds_test)

                if(append_save): # append the result each time we finish the simulation
                    # use a copy for saving. avoid concurrent iteration on same object
                    self.save_results(deepcopy(ds_test))
        
        print('all simulation(s) done, results saved in Test Data Set')        

    def post_process(self, test:TestItem, ds_exp:TestDataSet):
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

    def load_results(self, exp_name, dir_name=None, ds_exp:TestDataSet = None) -> TestDataSet:
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

                cfg_dict = helper.read_dict((u, col_start))
                u = u + len(cfg_dict) + offset
                ds_exp.set_cfgs(gname, cfg_dict)

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

    def post_process(self, test:TestItem, ds_exp:TestDataSet):
        '''
            post process the test result after simulation
            data mapping and re-collection can be applied in this function
        '''
 
        cfg = test.cfg
        result = test.result

        N_seg = int(cfg['N_seg'])
        L_fp = cfg['L_fp']
        len_seg = L_fp / N_seg

        T_hot = []
        T_cold = []
        gamma_hot = []
        gamma_cold = []
        hc_hot = []
        hc_cold = []
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

            hc_hot.append(result['HE.gasFlow.heatTransfer.hc[%s]' % idx].val)
            hc_cold.append(result['HE.fluidFlow.heatTransfer.hc[%s]' % idx].val)

            if i <= N_seg - 2: # len(gamma) = len(T) - 1
                gamma_hot.append(result['HE.gasFlow.heatTransfer.gamma[%s]' % idx].val)
                gamma_cold.append(result['HE.fluidFlow.heatTransfer.gamma[%s]' % idx].val)
        
        T_cold.reverse()
        dp_cold.reverse()  
        hc_cold.reverse()          
        
        x_values = np.arange(0, len_seg * (N_seg) , len_seg) 
        x_values_gamma = np.arange(0, len_seg * (N_seg - 1) , len_seg) 

        test.set_post_data('x_hot', x_values)
        test.set_post_data('T_hot',np.array(T_hot))
        test.set_post_data('dp_hot',np.array(dp_hot))
        test.set_post_data('x_cold', x_values + len_seg)
        test.set_post_data('T_cold', np.array(T_cold))
        test.set_post_data('dp_cold',np.array(dp_cold))
        test.set_post_data('x_values_gamma',np.array(x_values_gamma))  
        test.set_post_data('hc_hot',np.array(hc_hot))
        test.set_post_data('hc_cold',np.array(hc_cold))      
        test.set_post_data('gamma_hot',np.array(gamma_hot))
        test.set_post_data('gamma_cold',np.array(gamma_cold))

    def gen_sol_dict(self, cfg:dict) -> dict:
        ports = [
            ('HE.fluidFlow.heatTransfer', 'cold'),
            ('HE.fluidFlow.heatTransfer', 'cold'),
            ('HE.gasFlow.heatTransfer', 'hot'),
            ('HE.gasFlow.heatTransfer', 'hot')]

        # (propName, unit)
        N_seg =  cfg['N_seg']
        props = [
            # ('fluid', '1'), # error in parsing this props, since returned name(string) are not digitial values
            ('T[%s]', 'K', N_seg), 
            ('dp[%s]', 'Pa', N_seg),
            ('Re[%s]', '[1]', N_seg),
            ('rho[%s]', 'kg/m3', N_seg),
            ('hc[%s]', 'W/m2-K', N_seg),
            ('gamma[%s]', 'W/m2-K', N_seg - 1)]

        sol_dict = OrderedDict()

        for (p, p_short) in ports:
            for (prop, unit, n) in props:
                for idx_seg in range(1, n + 1):
                    prop_name = f'{prop}' % idx_seg
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

