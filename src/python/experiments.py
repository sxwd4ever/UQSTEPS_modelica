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

from model import TestDataSet, TestConfig, TestResult, TestItem, TestConstants
from plotlib import PlotManager, DataSeries, AxisType
from physics import Temperature, Pressure, MDot
from utils import ExcelHelper

class Experiments(object):
        
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
                if(append_save): # append the result each time we finish the simulation
                    # use the copy for saving. avoid concurrent iteration on same object
                    self.save_results(deepcopy(ds_test))
        
        print('all simulation(s) done, results saved in Test Data Set')        
    


    

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
                # a sheet for each group

                sht = wbk.sheets.add(name = group.name[0:30]) # 30 - max lenth for a sheet name
                ex: ExcelHelper = ExcelHelper(sht)
                u, l = 1, TestConstants.DATA_FILE_COL_START
                title = {"Table 1": "Test Config"}

                # config
                (b, r) = ex.write_table(cfg_set, title=title, up=u, left=l, linespacing=True)

                title = {"Table 2": "Results"}
                (b, r) = ex.write_table(result_set, title=title, up=b, left=l, linespacing=True)
                
                # if view_set: # not empty
                #     idx = 3
                #     for (view_name, data_table) in view_set.items():                    
                #         title = {f"Table {idx}": view_name}
                #         (b, r) = ex.write_table(data_table, title=title, up=b, left=l, linespacing=True)
                #         idx += 1

            from pathlib import Path

            ex_file=Path(filename)       
            if ex_file.exists():
                shutil.move(filename, filename+".bak")

            wbk.save(filename)                    

        finally:
            xw_app.quit() # close the xlwings app finally                                   

    
    def load_result(self, exp_name) -> TestDataSet:
        '''
        Load simulation result from the latest, saved file
        '''
        pass

    def gen_plot_manager(self, test: TestItem, plot_cfg:dict) -> PlotManager:
        #update_cal_fields() 
        pass

    def run(self,  simulate = True):
        pass

    def find_latest_test(self, path_) -> str:
        pass
