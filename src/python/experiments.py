from logging import exception
from math import trunc
from os.path import join, sep
import sys
from OMPython import OMCSessionZMQ
from OMPython import ModelicaSystem
from datetime import datetime as dt
from numpy.core import test
import pandas as pd
import xlwings as xw

from os import path,sep
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

    def simulate(self, sim_ops, solution_names, ds_test:TestDataSet):        
        
        self.model = ModelicaSystem(sep.join([self.path_model,"package.mo"]), self.model_name,self.modelica_libs) 

        # options for simulation - steady state simulation, no iteration required, so set numberOfIntervals = 2 

        for group in ds_test:
            for test in group:
                self.model.setParameters(test.cfg.gen_params())
                self.model.setSimulationOptions(sim_ops)
                self.model.simulate()
                solutions = self.model.getSolutions(solution_names) 

                test.result = TestResult.from_keyvalues(solution_names, solutions) 
        
 

        # for test in ds_test:
        #     cfg = test.cfg
        #     print('prepare to simulate:' + cfg.full_name)
        #     mod.setParameters(test.gen_sim_param())

        #     mod.simulate()

        #     print('simulation {0} done, retrieving result'.format(cfg.full_name))
        #     result = test.result

        #     # collect data in solutions
        #     for sol_key in result:
        #         sol = mod.getSolutions(sol_key)                    
        #         if not sol is None:
        #             result.set_result(sol_key, sol)                        
        #         else:
        #             print("solution with key = {0} not exsits".format(sol_key))

        #     # result.update_cal_columns()

        # print('all simulation(s) done, save result next')   

        # self.save_results(ds_test)    

    def save_results(self, ds_test : TestDataSet):
        '''
        save the simulation result into files    
        '''
        pass                                  

    def load_result(self, test_name) -> TestDataSet:
        '''
        Load simulation result from the latest, saved file
        '''
        pass

    def gen_plot_manager(self, test: TestItem) -> PlotManager:
        #update_cal_fields() 
        pass

    def run(self,  simulate = True):
        pass

    def find_latest_test(self, path_) -> str:
        pass
