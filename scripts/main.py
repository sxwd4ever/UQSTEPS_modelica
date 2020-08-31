from os.path import pathsep
from OMPython import OMCSessionZMQ
from OMPython import ModelicaSystem
from datetime import datetime as dt
import pandas as pd
import xlwings as xw

import csv
import os
from os import path
import numpy as np

from model import TestDataSet, TestConfig
from plotlib import PlotManager, DataSeries, AxisType
from physics import Temperature, Pressure, MDot
from utils import ExcelHelper

class TestPCHEMeshram(object):
        
    def __init__(self, work_root):
        super().__init__()
        # common variable defination
        self.path_pics = "pics"
        self.path_out = "out"

        self.work_root = work_root
        self.model_path_root = path.join(work_root,"Steps")
        self.prepare_workspace()        

 
    def prepare_workspace(self):
        '''
        prepare the workspace for simulation, the directory structure is
        - Steps (currernt work directory)
        |--.vscode  - directory for vscode's configurations
        |--docs     - documents
        |--build    - workspace for modelica simulation
        |--scripts  - python scripts for test or parameter sweep
        |--Steps    - Directory for modelica models codes
        |.env       - configuration file for vs python extension
        |.gitignore - git ignore file
        ''' 

        import os
        import shutil      

        os.chdir(self.work_root)

        # setup working directory
        pwd = "build" 
        if not os.path.exists(pwd):    
            os.mkdir(pwd)

        os.chdir(pwd)        

        if os.path.exists(self.path_pics):
            shutil.rmtree(self.path_pics)
        
        # os.mkdir(self.path_pics)

        # copy pics for plot use
        shutil.copytree("../scripts/pics", self.path_pics)

        # directory for output
        if not os.path.exists(self.path_out):
            os.mkdir(self.path_out)        

        # copy cool prop lib
        libs = ["libCoolProp.a", 'libCoolProp.dll']
        for lib in libs:            
            if not os.path.exists(lib):
                paths = ["Resources", "Library"]
                lib_path = self.model_path_root
                for path in paths:
                    lib_path = pathsep.join([lib_path, path])

                shutil.copyfile(lib_path + lib, pathsep.join([".",lib])) # completa target name needed 

    def run_simulation(self, ds_test : TestDataSet):

        omc = OMCSessionZMQ()
        # path = inspect.getfile(omc.__class__)

        mod = ModelicaSystem(path.join(self.model_path_root,"package.mo"),"Steps.Test.TestPCHXMeshram",["Modelica 3.2.1"])
        # options for simulation - steady state simulation, no iteration required, so set numberOfIntervals = 2 
        mod.setSimulationOptions('stepSize  = 0.2')

        for test in ds_test:
            cfg = test.cfg

            mod.setParameters(test.gen_sim_param())

            mod.simulate()

            result = test.result

            # collect data in solutions
            for sol_key in result:
                sol = mod.getSolutions(sol_key)                    
                if not sol is None:
                    result.set_result(sol_key, sol)                        
                else:
                    print("solution with key = {0} not exsits".format(sol_key))

            # result.update_cal_columns()

        print('Simulation done, save result next')   

        self.save_results(ds_test)    

    def save_results(self, ds_test : TestDataSet):
        '''
        save the simulation result into files    
        '''
        # create directory for current test batch
        name_batch = path.join(self.path_out, ds_test.name)

        xw_app: xw.App = xw.App(visible=False)

        wbs = []

        try:

            if not os.path.exists(name_batch):
                os.mkdir(name_batch)        

            for test in ds_test:
                cfg = test.cfg
                gname = cfg.group_name

                if gname not in wbs:
                    wbk = xw_app.books.add()
                    wbs.append(gname)
                else:
                    wbk = xw_app.books[gname]

                sht = wbk.sheets.add(name = cfg.name)

                ex: ExcelHelper = ExcelHelper(sht)

                title = {"Table 1": "Test Config of " + cfg.name}

                # config
                ex.write_dict(cfg.to_dict(), title_dict=title, linespacing=True)

                title = {"Table 2": "Test Result of " + cfg.name}
                ex.write_dict(test.result.to_dict(), title_dict= title, linespacing= True)


            for k in wbs:
                wb = xw_app.books[k]
                wb.save(os.path.join(name_batch, k))

        finally:
            xw_app.quit() # close the xlwings app finally

            # result = test.result

            # file_name = path.join(name_batch, "result_{0}.csv".format(cfg.name))

            # with open(file_name, mode='w', newline='\n') as csv_file:            

            #     writer = csv.DictWriter(csv_file, delimiter=',', quotechar='', fieldnames=['name','value'])

            #     writer.writeheader()

            #     for k in result:
            #         v = result.get_values(k)
            #         val_str = '{0}'.format(v)[1:-1] # remove the surrounding square brackets
            #         val_str = val_str.replace('[','').replace(']','').replace(' ', ',')
            #         writer.writerow({'name': k, 'value': val_str})                                

    def load_result(self) -> TestDataSet:
        '''
        Load simulation result from the latest, saved file
        '''
        import glob
        import os

        list_of_files = glob.glob(self.path_out + '/*.csv') # get all files ends with csv

        if list_of_files == []:
            raise FileNotFoundError("No result file saved in {0}/out".format(self.work_root))

        latest_file = max(list_of_files, key = os.path.getctime)

        if latest_file is None:
            raise FileNotFoundError("No results file saved in " + self.path_out)
        
        result_dict = {}

        with open(latest_file, mode = 'r', newline='\n') as csv_file:           
            reader = csv.DictReader(csv_file, delimiter=',', quotechar='"')
            for line in reader:
                # print(line['name'], line['value'])
                value_str = line['value']
                values = [float(x) for x in value_str.split(',')]

                result_dict[line['name']] = values[0]

        return result_dict

    def gen_plot_manager(self, result_dict):
        #update_cal_fields()   
        N_seg = self.param_des.N_seg
        len_seg = self.param_des.len_seg

        zigzag = 0 # 0: Straight Channel, 1: Zigzag Channel        
        axis_x=[[0, 0.12], [0, 0.16]][zigzag]
        axis_T=[[400.0, 750.0],[400, 750.0]][zigzag]
        axis_dp=[[0, 80], [0, 80]][zigzag]

        T_hot = []
        T_cold = []
        # global pressure drop R. T. the inlet temperature 
        dp_hot = [] 
        dp_cold = []
        # started, fix pressure of each stream
        p_hot_start = result_dict['pchx.cell_hot[1].p']
        p_cold_start = result_dict['pchx.cell_cold[{0}].p'.format(N_seg)]

        for i in range(0, N_seg):
          T_hot.append(result_dict['pchx.cell_hot[{0}].T'.format(i+1)])
          T_cold.append(result_dict['pchx.cell_cold[{0}].T'.format(i+1)])
          # pa -> kPa
          dp_hot.append((p_hot_start - result_dict['pchx.cell_hot[{0}].p'.format(i+1)]) / 1e3)
          dp_cold.append((p_cold_start - result_dict['pchx.cell_cold[{0}].p'.format(i+1)]) / 1e3)
        
        x_values = np.arange(0, len_seg * (N_seg) , len_seg) 

        label_x = 'Z (m)'
        label_T = 'T (K)'
        label_dp = '$\Delta~p~(kPa)$ (K)'

        plot = PlotManager(title='Zigzag Channel @ High Temperature Range')
        plot.add(DataSeries(name = 'T_hot', x = x_values, y = np.array(T_hot), range_x=axis_x, range_y=axis_T, cs = 'r-s', label=[label_x, label_T]))
        plot.add(DataSeries(name = 'T_cold', x = x_values + len_seg, y = np.array(T_cold), range_x=axis_x, range_y=axis_T, cs = 'b--s', label=[label_x, label_T]))
        plot.add(DataSeries(name = 'dp_hot', x = x_values, y = np.array(dp_hot), range_x=axis_x, range_y=axis_dp, cs = 'r-^', label=[label_x, label_dp]), ax_type=AxisType.Secondary)
        plot.add(DataSeries(name = 'dp_cold', x = x_values + len_seg, y = np.array(dp_cold), range_x=axis_x, range_y=axis_dp, cs = 'b--^', label=[label_x, label_dp]), ax_type=AxisType.Secondary)

        return plot

    def run(self, ds_test : TestDataSet, simulate = True):

        if simulate:                          
            self.run_simulation(ds_test)            
            print('simulation done, ready to plot figures')
        else:
            # load the simulation result from latest, pre-saved file
            ds_test = self.load_result()

        zigzag = 1
        imgfile = [self.path_pics + "/Meshram_Fig_04.png", self.path_pics + "/Meshram_Fig_05.png"][zigzag]    
        
        self.gen_plot_manager(ds_test).draw(img_file=imgfile, dest_file= self.path_out + "/Meshram_Fig%db_compare.png"%(4 + zigzag))

        print('all done!') 

def main(work_root = []):
    # root path of modelica root
    if work_root == []:
        work_root = os.path.abspath(os.curdir)

    cfg_ref = TestConfig()

    ds_test = TestDataSet(name = dt.now().strftime("Test_%Y_%m_%d_%H_%M_%S"), cfg=cfg_ref)

    g = ds_test.new_para_group("zigzag_HT", ["Single"])
    g.add_para_seq("Re_des", [5000])
    g.submit() # only submited group will be valid for test

    g = ds_test.new_para_group("zigzag_HT_diff_Re", ["Re_des = 5000", "Re_des = 20000"])
    g.add_para_seq("Re_des", [5000, 20000])    

    g = ds_test.new_para_group("zigzag_HT_diff_mdot", [r"$\dot{m}_{off}$ = 10", r"${\dot{m}_{off}$ = 100"])
    g.add_para_seq("mdot_hot_odes", [MDot(10), MDot(100)])
    g.add_para_seq("mdot_cold_odes", [MDot(10), MDot(100)])
    # g.submit()

    g = ds_test.new_para_group("zigzag_HT_diff_pT_in", [r"$(p,T)_{hi}$ = 10 MPa, 450°", r"$(p,T)_{hi}$ = 12 MPa, 300°"])
    g.add_para_seq("p_hot_in", [Pressure.MPa(10), Pressure(12)])
    g.add_para_seq("T_hot_in", [Temperature.degC(730), Temperature.degC(300)])   
    # g.submit

    test = TestPCHEMeshram(work_root)

    test.run(ds_test=ds_test, simulate=True)    

###
if __name__ == "__main__":
    main()