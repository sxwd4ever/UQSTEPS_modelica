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

        for group in ds_test:
            for test in group:
                sol_keys = list(solution_dict.keys())
                solutions = np.ones((len(sol_keys),1)) * -1 # default value
                
                try:
                    print(f' ready to run test={test.name} in group={group.name}\n')
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
        '''
        # create directory for current test batch
        name_batch = path.join(self.path_out, ds_test.name)

        xw_app: xw.App = xw.App(visible=False)
        wbs = {} # set of names of workbook

        try:

            if not os.path.exists(name_batch):
                os.mkdir(name_batch)        

            for group in ds_test:                
                gname = group.name
                filename = os.path.join(name_batch, gname + '.xlsx')
                
                if gname not in wbs:   

                    wbk = xw_app.books.add()
                    wbs[gname] = wbk.name
                else:
                    wbk = xw.Book(filename)

                for test in group:
                    cfg = test.cfg
                    sht = wbk.sheets.add(name = cfg.name[0:30]) # 30 - max lenth for a sheet name

                    ex: ExcelHelper = ExcelHelper(sht)
                    u, l = 1, TestConstants.DATA_FILE_COL_START
                    title = {"Table 1": "Test Config"}

                    # config
                    (b, r) = ex.write_table(cfg.to_dict(), title=title, up=u, left=l, linespacing=True)

                    title = {"Table 2": "Results"}
                    (b, r) = ex.write_table(test.result.to_dict(), title=title, up=b, left=l, linespacing=True)

                    wbk.save(filename)
                    # wbk.close()

        finally:
            xw_app.quit() # close the xlwings app finally                                     

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
