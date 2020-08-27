from OMPython import OMCSessionZMQ
from OMPython import OMCSession
from OMPython import ModelicaSystem
from datetime import datetime as dt
from enum import IntEnum
from CoolProp.CoolProp import PropsSI

import inspect
import csv
import os
import numpy as np
import math

from pbmodel import DesignParam, OffDesignParam, StreamType
from plotlib import PlotManager, DataSeries, AxisType
from physics import Temperature, Pressure, MDot, Velocity, Density, Angle, Length
from pbmodel import ParamSweepSet, TestConfig

class TestPCHEMeshram(object):
        
    def __init__(self, work_root):
        super().__init__()
        # common variable defination
        self.path_pics = "pics"
        self.path_out = "out"

        self.work_root = work_root
        self.model_path_root = work_root + r"/Steps" 
        self.prepare_workspace()

        # params initialization
        self.param_des = None

        # Off design params
        self.param_odes = None
    
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
                shutil.copyfile(self.model_path_root + r"\Resources\Library\\" + lib, ".\\" + lib) # completa target name needed 

    def __gen_result_dict(self):
        
        # dict to store solutions under 
        result_dict = {}

        # for all the component's inlet/outlet that needs to be analyzed.
        ports_keys = [
            'source_hot.outlet', 
            'source_cold.outlet', 
            'sink_hot.inlet', 
            'sink_cold.inlet', 
            'pchx.inlet_hot',
            'pchx.outlet_hot',
            'pchx.inlet_cool',
            'pchx.outlet_cool'        
            ]
        
        # common variables for these port
        var_keys = ['h_outflow', 'p', 'm_flow']

        # fill in keys of result dict with above 
        for p_key in ports_keys:
            for val_key in var_keys:
                sol_key = p_key + '.' + val_key
                result_dict[sol_key] = []

        # for each HX cell in PCHE
        N_seg = 10 
        var_keys = ['T', 'p', 'h', 'u', 'k', 'rho', 'mu' ,'dp', 'G', 'Re', 'Nu', 'f', 'Q']
        node_keys = ['cell_cold', 'cell_hot']
        for var_key in var_keys:
            for node_key in node_keys:
                for i in range(1, N_seg + 1):        
                    result_dict["pchx.{node}[{idx}].{var}".format(idx = i, var = var_key, node = node_key)] = []

        # some special varible in the model, specify them explicitly like
        # input parameters
        result_dict['pchx.N_seg'] = []  
        result_dict['pchx.length_cell'] = [] 
        result_dict['pchx.phi'] = []  
        result_dict['pchx.d_c'] = []  
        result_dict['pchx.pitch'] = []  
        result_dict['pchx.kim_cor.a'] = []  
        result_dict['pchx.kim_cor.b'] = [] 
        result_dict['pchx.kim_cor.c'] = []  
        result_dict['pchx.kim_cor.d'] = []          

        # calculated initial parameters
        result_dict['pchx.d_h'] = []  
        result_dict['pchx.peri_c'] = [] 
        result_dict['pchx.t_wall'] = []  
        result_dict['pchx.N_ch'] = [] 
        result_dict['pchx.A_c'] = []  
        result_dict['pchx.A_flow'] = [] 
 
        result_dict['pchx.A_stack'] = []                  
        result_dict['pchx.length_ch'] = []     
        result_dict['pchx.Re_hot_start'] = [] 
        result_dict['pchx.Re_cold_start'] = []                           
        result_dict['mdot_hot'] = []  
        result_dict['mdot_cold'] = []  

        return result_dict

    def __cal_row(self, vals = [], eval_str = 'x + y'):

        real_eval_str = eval_str.replace('x','op1').replace('y','op2').replace('z', 'op3')

        r = []
        
        for i in range(0, len(vals[0])):
            op1,op2 = vals[0][i], vals[1][i]

            if(len(vals) == 3):
                op3 = vals[2][i]
            r.append(eval(real_eval_str))

        return r

    def __gen_cal_map(self):
        '''
        map containing all the calculated fields, which are calculated 
        according to the solutions of model
        '''
        # map for calculated values
        cal_map = []
        # recheck efficiency 
        cal_map.append({'key' : 'Q_in', 'vals':['pcm_heater.h_e', 'pcm_heater.h_i', 'pcm_heater.inlet.m_flow'], 'eval_str': 'z * (x - y)'})
        cal_map.append({'key' : 'W_turbine', 'vals':['turbine.h_i', 'turbine.h_ea', 'turbine.inlet.m_flow'], 'eval_str': 'z * (x - y)'})
        cal_map.append({'key' : 'W_comp', 'vals':['pump.h_ea', 'pump.h_i', 'pump.inlet.m_flow'], 'eval_str': 'z * (x - y)'})
        cal_map.append({'key' : 'W_comp_recom', 'vals':['recom_pump.h_ea', 'recom_pump.h_i', 'recom_pump.inlet.m_flow'], 'eval_str': 'z * (x - y)'})
        cal_map.append({'key' : 'W_net', 'vals':['W_turbine', 'W_comp', 'W_comp_recom'], 'eval_str': 'x - y - z'})
        cal_map.append({'key' : 'eta_all_recal', 'vals':['W_net', 'Q_in'], 'eval_str': 'x / y * 100'})
        cal_map.append({'key' : 'eta_all_recal', 'vals':['W_net', 'Q_in'], 'eval_str': 'x / y * 100'})

        # performance
        # LTR
        cal_map.append({'key' : 'LTR_cool_dt', 'vals':['recup_low.crec_out.T','recup_low.crec_in.T'], 'eval_str': 'x - y'})
        cal_map.append({'key' : 'LTR_hot_dt', 'vals':['recup_low.inlet.T','recup_low.outlet.T'], 'eval_str': 'x - y'})
        cal_map.append({'key' : 'LTR_cool_dQ', 'vals':['recup_low.h_cool_e','recup_low.h_cool_i', 'recup_low.crec_in.m_flow'], 'eval_str': '(x - y) * z'})
        cal_map.append({'key' : 'LTR_hot_dQ', 'vals':['recup_low.h_hot_i','recup_low.h_hot_e', 'recup_low.inlet.m_flow'], 'eval_str': '(x - y) * z'})

        # HTR
        cal_map.append({'key' : 'HTR_cool_dt', 'vals':['recup_high.crec_out.T','recup_high.crec_in.T'], 'eval_str': 'x - y'})
        cal_map.append({'key' : 'HTR_hot_dt', 'vals':['recup_high.inlet.T','recup_high.outlet.T'], 'eval_str': 'x - y'})
        cal_map.append({'key' : 'HTR_cool_dQ', 'vals':['recup_high.h_cool_e','recup_high.h_cool_i', 'recup_high.crec_in.m_flow'], 'eval_str': '(x - y) * z'})
        cal_map.append({'key' : 'HTR_hot_dQ', 'vals':['recup_high.h_hot_i','recup_high.h_hot_e', 'recup_high.inlet.m_flow'], 'eval_str': '(x - y) * z'})

        return cal_map

    def gen_model_param(self):
        # set up on design parameters
        para = self.param_des
        (hot, cold) = (StreamType.Hot, StreamType.Cold)
        para_dict = {}
        para_dict["N_ch"] = para.N_ch

        para_dict["T_cold_in"] = para.T[cold]
        para_dict["p_cold_in"] = para.p[cold]

        para_dict["T_hot_in"] = para.T[hot]
        para_dict["p_hot_in"] = para.p[hot]

        para_dict["phi"] = para.phi
        para_dict["pitch"] = para.pitch
        para_dict["d_c"] = para.d_c
        para_dict["length_cell"] = para.len_seg       

        # set up off design parameters
        para = self.param_odes
        para_dict["Re_hot_start"] = para.Re[hot]
        para_dict["Re_cold_start"] = para.Re[cold]
        
        para_dict["mdot_hot"] = para.mdot[hot]
        para_dict["mdot_cold"] = para.mdot[cold]

        params = []

        for k, v in para_dict.items():
            params.append("{0}={1}".format(k, str(v)))

        return params

    def run_simulation(self, cfg:TestConfig, param_set, result_dict):

        omc = OMCSessionZMQ()
        # path = inspect.getfile(omc.__class__)

        mod = ModelicaSystem(self.model_path_root + r"\package.mo","Steps.Test.TestPCHXMeshram",["Modelica 3.2.1"])
        # options for simulation - steady state simulation, no iteration required, so set numberOfIntervals = 2 
        mod.setSimulationOptions('stepSize  = 0.2')

        clones = param_set.gen_configs(cfg)

        for clone in clones:
            (self.param_des, self.param_odes) = clone.gen_test_param()
            pars = self.gen_model_param()

            mod.setParameters(pars)

            mod.simulate()

            # collect data in solutions
            for sol_key in result_dict.keys():
                l = result_dict[sol_key]
                sol = mod.getSolutions(sol_key)
                if not sol is None:
                    val = mod.getSolutions(sol_key)[0][0]
                    l.append(val)
                else:
                    print("solution with key = {0} not exsits".format(sol_key))

    def __update_cal_fields(self, result_dict):
        # fill calculated fields into the result_dict
        cal_map = self.__gen_cal_map()
        # fill in the calculated value by python's eval
        for item in cal_map:
            vals = item['vals']

            data = result_dict

            if len(vals) == 2:
                result_dict[item['key']] = __cal_row([data[vals[0]], data[vals[1]]], eval_str=item['eval_str'])
            elif len(vals) == 3:
                result_dict[item['key']] = __cal_row([data[vals[0]], data[vals[1]], data[vals[2]]], eval_str=item['eval_str'])      

    def save_results(self, result_dict, file_name = []):
        '''
        save the simulation result into files    
        '''

        if file_name == [] :
            file_name = self.path_out + "/result_%Y_%m_%d_%H_%M_%S.csv"

        with open(dt.now().strftime(file_name) ,mode='w', newline='\n') as csv_file:            

            writer = csv.DictWriter(csv_file, delimiter=',', quotechar='"', fieldnames=['name','value'])

            writer.writeheader()

            for k,v in result_dict.items():
                val_str = '{0}'.format(v)[1:-1] # remove the surrounding square brackets
                writer.writerow({'name': k, 'value': val_str})                                

    def load_result(self):
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

    def run(self, cfg: TestConfig, param_set: ParamSweepSet, simulate = True):

        if simulate:
            result_dict = self.__gen_result_dict()              

            self.run_simulation(cfg, param_set,result_dict)

            self.save_results(result_dict)
        
        # load the simulation result from latest, pre-saved file
        result_dict = self.load_result()

        zigzag = 1
        imgfile = [self.path_pics + "/Meshram_Fig_04.png", self.path_pics + "/Meshram_Fig_05.png"][zigzag]    
        
        self.gen_plot_manager(result_dict).draw(img_file=imgfile, dest_file= self.path_out + "/Meshram_Fig%db_compare.png"%(4 + zigzag))

        print('all done!') 

def main(work_root = []):
    # root path of modelica root
    if work_root == []:
        work_root = os.path.abspath(os.curdir)

    cfg = TestConfig()

    para_sweep_set = ParamSweepSet()

    g = para_sweep_set.new_group("group1", ["Re_des = 5000", "Re_des = 20000"])
    g.add_para_seq("Re_des", [5000, 20000])


    g = para_sweep_set.new_group("group2", [r"$\dot{m}_{off}$ = 10", r"${\dot{m}_{off}$ = 100"])
    g.add_para_seq("mdot_hot_odes", [MDot(10), MDot(100)])
    g.add_para_seq("mdot_cold_odes", [MDot(10), MDot(100)])

    g = para_sweep_set.new_group("group3", [r"$(p,T)_{hi}$ = 10 MPa, 450°", r"$(p,T)_{hi}$ = 12 MPa, 300°"])
    g.add_para_seq("p_hot_in", [Pressure.MPa(10), Pressure(12)])
    g.add_para_seq("T_hot_in", [Temperature.degC(730), Temperature.degC(300)])    

    test = TestPCHEMeshram(work_root)

    test.run( cfg, para_sweep_set, simulate=True)    

###
if __name__ == "__main__":
    main()