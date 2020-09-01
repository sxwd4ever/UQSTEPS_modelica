from math import trunc
from os.path import join, sep
from OMPython import OMCSessionZMQ
from OMPython import ModelicaSystem
from datetime import datetime as dt
from numpy.core import test
import pandas as pd
import xlwings as xw

import csv
import os
from os import path,sep
import numpy as np
from ex.CoolPropQuery import result

from model import TestDataSet, TestConfig, TestResult, TestDataItem, TestConstants
from plotlib import PlotManager, DataSeries, AxisType
from physics import Temperature, Pressure, MDot
from utils import ExcelHelper

class TestPCHEMeshram(object):
        
    def __init__(self, work_root):
        super().__init__()
        # common variable defination
        self.path_root = work_root

        self.path_model = sep.join([self.path_root, "Steps"])
        self.path_workspace = sep.join([self.path_root, "build"])
        self.path_pics = sep.join([self.path_workspace, "pics"])
        self.path_out = sep.join([self.path_workspace, "out"])

        self.prepare_workspace()        

 
    def prepare_workspace(self):
        '''
        prepare the workspace for simulation, the directory structure is
        - Steps (currernt work directory) -> path_root
        |--.vscode              - directory for vscode's configurations
        |--docs                 - documents
        |--build                - -> path_workspace, workspace for modelica simulation, all
        |  |                        the relative paths in this program are relative to point
        |  |--out               - Final output of the Simulation
        |  |  |--Test_Batch 1   - Simulation output data for test 1
        |  |  ...               - Simulation output data for test ...
        |  |--pics              - Temp storage dir for referenced pic 
        |--scripts              - python scripts for test or parameter sweep
        |--Steps                - Directory for modelica models codes
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
        shutil.copytree(sep.join([self.path_root, 'scripts', 'pics']), self.path_pics)

        # directory for output        
        if not os.path.exists(self.path_out):
            os.mkdir(self.path_out)        

        # copy cool prop lib
        libs = ["libCoolProp.a", 'libCoolProp.dll']
        for lib in libs: 
            lib_path = sep.join([self.path_model, "Resources", "Library", lib])           
            if not os.path.exists(lib_path):
                shutil.copyfile(lib_path, sep.join([".",lib])) # completa target name needed 

    def simulate(self, ds_test : TestDataSet):

        omc = OMCSessionZMQ()
        # path = inspect.getfile(omc.__class__)

        mod = ModelicaSystem(sep.join([self.path_model,"package.mo"]),"Steps.Test.TestPCHXMeshram",["Modelica 3.2.1"])
        # options for simulation - steady state simulation, no iteration required, so set numberOfIntervals = 2 
        mod.setSimulationOptions('stepSize  = 0.2')

        for test in ds_test:
            cfg = test.cfg
            print('prepare to simulate=' + cfg.full_name)
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
        wbs = {} # set of names of workbook

        try:

            if not os.path.exists(name_batch):
                os.mkdir(name_batch)        

            for test in ds_test:
                cfg = test.cfg
                gname = cfg.group_name
                filename = os.path.join(name_batch, gname + '.xlsx')

                if gname not in wbs:
                    wbk = xw_app.books.add()                                     
                    wbs[gname] = wbk.name
                else:
                    wbk = xw.Book(filename)

                sht = wbk.sheets.add(name = cfg.name[0:30]) # 30 - max lenth for a sheet name

                ex: ExcelHelper = ExcelHelper(sht)
                u, l = 1, TestConstants.DATA_FILE_COL_START
                title = {"Table 1": "Test Config"}

                # config
                (b, r) = ex.write_table(cfg.to_dict(), title=title, up=u, left=l, linespacing=True)

                title = {"Table 2": "Results"}
                (b, r) = ex.write_table(test.result.to_dict(), title=title, up=b, left=l, linespacing=True)

                wbk.save(filename)
                wbk.close()

        finally:
            xw_app.quit() # close the xlwings app finally                                         

    def load_result(self, test_name) -> TestDataSet:
        '''
        Load simulation result from the latest, saved file
        '''
        import glob
        import os
        import xlwings as xw

        # find test_files 
        test_files = glob.glob(sep.join([self.path_out, test_name, '*.xlsx'])) # get all files ends with xlsx

        if test_files == []:
            raise FileNotFoundError("No result file saved in {0}/{1}".format(self.path_out,test_name))

        # filter temporary files 
        test_files = list(filter( lambda x: x.rfind('\\~$') < 0, test_files))

        app = xw.App(visible=False)
        ds_test:TestDataSet

        try:
            ds_test = TestDataSet(name = test_name,cfg=TestConfig())

            for test_file in test_files:
                    
                wbk:xw.Book = xw.Book(test_file,read_only=True)
                col_start = TestConstants.DATA_FILE_COL_START

                for sht in wbk.sheets:
                    if sht.name.startswith('Sheet'):
                        continue

                    helper = ExcelHelper(sht)
                    u = 1
                    offset = 2 # + 1 for linespacing, + 1 is start of next table 

                    dict_ = helper.read_dict((u, col_start))

                    cfg = TestConfig(name=sht.name, group_name=sht.book.name.replace('.xlsx',''))
                    cfg.from_dict(dict_)

                    u = u + len(dict_) + offset

                    dict_ = helper.read_dict((u, col_start))
                    (p_des, p_odes) = cfg.gen_test_param()

                    result = TestResult(sht.name, p_des)
                    result.from_dict(dict_)

                    test_item = TestDataItem(cfg)
                    test_item.result = result

                    ds_test.addTestDataItem(test_item)

                break

        finally:
            app.quit()

        return ds_test

    def gen_plot_manager(self, test: TestDataItem) -> PlotManager:
        #update_cal_fields() 
        cfg = test.cfg
        (p_des, p_odes) = cfg.gen_test_param()
        N_seg = int(p_des.N_seg)
        len_seg = p_des.len_seg

        result = test.result

        zigzag = 1 # 0: Straight Channel, 1: Zigzag Channel        
        # axis_x=[[0, 0.12], [0, 0.16]][zigzag]
        axis_x=[[0, 0.12], [0, 0.12]][zigzag]
        axis_T=[[400.0, 750.0],[400, 750.0]][zigzag]
        axis_dp=[[0, 80], [0, 80]][zigzag]

        T_hot = []
        T_cold = []
        # global pressure drop R. T. the inlet temperature 
        dp_hot = [] 
        dp_cold = []
        # started, fix pressure of each stream
        p_hot_start = result.get_values('pchx.cell_hot[1].p')[0]
        p_cold_start = result.get_values('pchx.cell_cold[{0}].p'.format(N_seg))[0]

        for i in range(0, N_seg):
            T_hot.append(result.get_values('pchx.cell_hot[{0}].T'.format(i+1))[0])
            T_cold.append(result.get_values('pchx.cell_cold[{0}].T'.format(i+1))[0])
            # pa -> kPa
            dp_hot.append((p_hot_start - result.get_values('pchx.cell_hot[{0}].p'.format(i+1))[0]) / 1e3)
            dp_cold.append((p_cold_start - result.get_values('pchx.cell_cold[{0}].p'.format(i+1))[0]) / 1e3)
        
        x_values = np.arange(0, len_seg * (N_seg) , len_seg) 

        label_x = 'Z (m)'
        label_T = 'T (K)'
        label_dp = '$\\Delta~p~(kPa)$ (K)'

        plot = PlotManager(title=cfg.full_name)
        plot.add(DataSeries(name = 'T_hot', x = x_values, y = np.array(T_hot), range_x=axis_x, range_y=axis_T, cs = 'r-s', label=[label_x, label_T]))
        plot.add(DataSeries(name = 'T_cold', x = x_values + len_seg, y = np.array(T_cold), range_x=axis_x, range_y=axis_T, cs = 'b--s', label=[label_x, label_T]))
        plot.add(DataSeries(name = 'dp_hot', x = x_values, y = np.array(dp_hot), range_x=axis_x, range_y=axis_dp, cs = 'r-^', label=[label_x, label_dp]), ax_type=AxisType.Secondary)
        plot.add(DataSeries(name = 'dp_cold', x = x_values + len_seg, y = np.array(dp_cold), range_x=axis_x, range_y=axis_dp, cs = 'b--^', label=[label_x, label_dp]), ax_type=AxisType.Secondary)

        return plot

    def run(self,  simulate = True):
        
        pk = TestConstants

        if simulate:
            # run simulation before loading results
            cfg_ref = TestConfig()

            ds_test = TestDataSet(cfg=cfg_ref)

            g = ds_test.new_para_group("zigzag_HT", ["Single"])
            g.add_para_seq(pk.Re_des, [5000])
            # g.submit() # only submited group will be valid for test

            g = ds_test.new_para_group("zigzag_HT_diff_Re", ["Re_des=5000", "Re_des=20000"])
            g.add_para_seq(pk.Re_des, [5000, 20000])    
            g.submit()

            g = ds_test.new_para_group("zigzag_HT_diff_mdot", [r"$\dot{m}_{off}$=50（-50%)"])
            g.add_para_seq(pk.mdot_hot_odes, [MDot(50)])
            g.add_para_seq(pk.mdot_cold_odes, [MDot(50)])
            g.submit()

            g = ds_test.new_para_group("zigzag_HT_diff_mdot", [r"$\dot{m}_{off}$=50（-50%)", r"${\dot{m}_{off}$=100(50%)"])
            g.add_para_seq(pk.mdot_hot_odes, [MDot(50), MDot(150)])
            g.add_para_seq(pk.mdot_cold_odes, [MDot(50), MDot(150)])            

            g = ds_test.new_para_group("zigzag_HT_diff_pT_in", [r"$(p,T)_{hi}$=10 MPa, 450°", r"$(p,T)_{hi}$=12 MPa, 300°"])
            g.add_para_seq(pk.p_hot_in, [Pressure.MPa(10), Pressure(12)])
            g.add_para_seq(pk.T_hot_in, [Temperature.degC(730), Temperature.degC(300)])   
            # g.submit

            self.simulate(ds_test)  
            print('Test data saved, ready to plot figures')

        # loading the newes result
        test_name = self.find_latest_test(self.path_out)                 
        # load the simulation result from latest, pre-saved file
        ds_test = self.load_result(test_name)

        zigzag = 1
        imgfile = ["Meshram_Fig_04.png", "Meshram_Fig_05.png"][zigzag]    
        imgfile = sep.join([self.path_pics, imgfile])

        for test in ds_test:
            destfile = sep.join([self.path_out, ds_test.name, "Meshram_Fig{0}b_compare_{1}.png".format(4 + zigzag, test.name)])
            self.gen_plot_manager(test).draw(img_file=imgfile, dest_file=destfile)

        print('all done!') 

    def find_latest_test(self, path_) -> str:
        lists = [x[0] for x in os.walk(path_)]
        dir_:str = max(lists, key = os.path.getctime)
        # -1 is the return if no substring found
        # do not use rindex() in case raising errors
        idx = dir_.rfind('\\') 
        return dir_[idx + 1:]

def main(work_root = []):
    # root path of modelica root
    if work_root == []:
        work_root = os.path.abspath(os.curdir)

    test = TestPCHEMeshram(work_root)

    ds_test = test.run(simulate=True)
###
if __name__ == "__main__":
    main()