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

from plot_lib import PlotManager, DataSeries, AxisType

class StreamType(IntEnum):
    Hot = 0,
    Cold = 1

class DesignParam(object):
    '''
    Structure containing On Design parameters
    '''

    def __init__(self, d_c = 12e-3, p=[9e6, 20e6], T=[500, 300], mdot = [8.3, 8.3], Re = 2000, N_seg = 10, len_seg=12e-3, pitch = 12e-3, phi= math.pi / 4):
        # On design params initialization 
        self.d_c = d_c
        self.p = np.array(p)
        self.T = np.array(T)
        self.mdot = np.array(mdot)
        self.Re = np.array(Re)
        self.N_seg = N_seg
        self.len_seg = len_seg  
        self.pitch = pitch
        self.phi = phi

        # parameters determined in cal_on_design_params
        self.d_h = 0.0
        self.A_c = 0.0
        self.A_f = 0.0
        self.N_ch = 0.0
        self.mu = np.zeros(len(self.mdot))  
        self.G = np.zeros(len(self.mdot)) # mass flux kg/(m^2 * s) should be constant if Re, P, T, d_c are constant

        self.cal_design_params()

    def cal_geo_params(self):
        d_c = self.d_c
        A_c = math.pi * d_c * d_c / 8 
        peri_c = d_c * math.pi / 2 + d_c          
        d_h = 4 * A_c / peri_c

        return (A_c, d_h, peri_c)

    def area_flow(self):
        return self.A_c * self.N_ch

    def cal_design_params(self):
        
        #  0 = hot, 1 = cold, following Enum PipeType
        A_fx = np.zeros(len(StreamType))

        (self.A_c, self.d_h, peri_c) = self.cal_geo_params()

        for i in StreamType:
            self.mu[i] = PropsSI('V', 'P', self.p[i], 'T', self.T[i], "CO2")
            
        A_fx = np.divide(self.mdot, self.mu) * (self.d_h / self.Re)   

        self.A_f = max(A_fx)
        self.N_ch = math.ceil(self.A_f / self.A_c)
        self.G = self.mdot / self.A_f

class OffDesignParam(object):

    def __init__(self, param_des: DesignParam, mdot = [], G = []):

        if mdot == []:
            if G == []:
                raise ValueError('no mdot or G assigned for off design')
            else:
                mdot = np.array(G) * param_des.area_flow()     

        num_stream = len(StreamType)

        self.mdot = np.array(mdot)
        self.Re = np.zeros(num_stream)
        self.G = np.zeros(num_stream)
        self.param_des = param_des
        self.__update()

    def __update(self):
        p_des = self.param_des
        A_f = self.param_des.area_flow() 
        self.G = self.mdot / A_f
        self.Re = np.divide(self.mdot, p_des.mu) * p_des.d_h / A_f

class TestPCHEMeshram(object):
        
    def __init__(self, work_root, param_des = None, param_odes = []):
        super().__init__()
        # common variable defination
        self.path_pics = "pics"
        self.path_out = "out"

        self.work_root = work_root
        self.model_path_root = work_root + r"/Steps" 
        self.prepare_workspace()

        # params initialization
        if param_des == None:
            self.param_des = DesignParam()
        else:
            self.param_des = param_des

        # Off design params
        self.param_odes = param_odes
    
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

    def run_simulation(self, result_dict):

        omc = OMCSessionZMQ()
        # path = inspect.getfile(omc.__class__)

        mod = ModelicaSystem(self.model_path_root + r"\package.mo","Steps.Test.TestPCHXMeshram",["Modelica 3.2.1"])
        # options for simulation - steady state simulation, no iteration required, so set numberOfIntervals = 2 
        mod.setSimulationOptions('stepSize  = 0.2')

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

        plot = PlotManager()

        plot.add(DataSeries(name = 'T_hot', x = x_values, y = np.array(T_hot), range_x=axis_x, range_y=axis_T, cs = 'r-s'))
        plot.add(DataSeries(name = 'T_cold', x = x_values + len_seg, y = np.array(T_cold), range_x=axis_x, range_y=axis_T, cs = 'b-s'))
        plot.add(DataSeries(name = 'dp_hot', x = x_values, y = np.array(dp_hot), range_x=axis_x, range_y=axis_dp, cs = 'r-^'), ax_type=AxisType.Secondary)
        plot.add(DataSeries(name = 'dp_cold', x = x_values + len_seg, y = np.array(dp_cold), range_x=axis_x, range_y=axis_dp, cs = 'b-^'), ax_type=AxisType.Secondary)

        return plot

    def run(self, simulate = True):

        if simulate:
            result_dict = self.__gen_result_dict()              

            self.run_simulation(result_dict)

            self.save_results(result_dict)
        
        # load the simulation result from latest, pre-saved file
        result_dict = self.load_result()

        zigzag = 1
        imgfile = [self.path_pics + "/Meshram_Fig_04.png", self.path_pics + "/Meshram_Fig_05.png"][zigzag]    
        
        self.gen_plot_manager(result_dict).draw(img_file=imgfile, dest_file= self.path_out + "/Meshram_Fig%db_compare.png"%(4 + zigzag))

        print('all done!') 
class TestScenario(IntEnum):
    '''
    TestScenario Index for Meshram [2016]
    '''
    straight_low_T = 0,
    straight_high_T = 1,
    zigzag_low_T = 2,
    zigzag_high_T = 3

class Unit(object):
    '''
    class for unit conversion
    # IMPORTANT: SIunits will be used in this script : [T] = K, [p] = pa, [length] = m, [angle] = rad
    '''
    from pint import UnitRegistry
    ureg = UnitRegistry()

    @staticmethod
    def convert(val, src_unit, dest_unit):
        ureg = Unit.ureg
        ureg.Unit
        src = val * ureg[src_unit]
        return src.to(ureg[dest_unit]).magnitude

    @staticmethod
    def from_bar(val):
        return Unit.convert(val, 'bar', 'Pa')

    @staticmethod
    def from_degC(val):
        return Unit.convert(val, '°C','K')

    @staticmethod
    def from_deg(val): 
        '''
        deg to rad for angle
        '''
        return Unit.convert(val, '°', 'rad')

def main(work_root = []):
    # root path of modelica root
    if work_root == []:
        work_root = os.path.abspath(os.curdir)

    # index of test scenario
    scenario = TestScenario.zigzag_high_T    
  
    # configuration of different test scenario in meshram [2016] - table 3
    # arranged in same order as Enum TestScenario

    # IMPORTANT: SIunits will be used in this script : [T] = K, [p] = pa, [length] = m, [angle] = rad, [mass] = kg, [time] = s   

    T_cold_in = [400, 500, 400, 500][scenario] # unit K
    p_cold_in = Unit.from_bar(225)
    
    T_cold_out = [498.45, 615.48, 522.23, 639.15][scenario]
    p_cold_out = Unit.from_bar(225)

    T_hot_in = [630, 730, 630, 730][scenario]
    p_hot_in = Unit.from_bar(90)

    T_hot_out = [494.37, 601.83, 466.69, 576.69][scenario]
    p_hot_out = Unit.from_bar(90)
    
    mdot = np.array([10, 10])

    phi = [Unit.from_deg(0), Unit.from_deg(0), Unit.from_deg((180 - 108) /2), Unit.from_deg((180 - 108) /2)][scenario]

    param_des = DesignParam(d_c=2e-3, p=[p_hot_in, p_cold_in], T=[T_hot_in, T_cold_in], Re=18000, mdot=mdot, N_seg = 10, len_seg=12e-3, pitch = 12e-3, phi = phi)

    # index of array for hot/cold stream
    # 0 hot stream, 1 cold stream
    mdot_odes = np.array([10, 10])
    # array to store meshram's data     
    p = np.array([p_hot_in, p_cold_in])
    T = np.array([T_hot_in, T_cold_in])
    media_name = ['CO2', 'CO2']
    u_odes = np.array([7.564, 1.876]) # off design velocity for hot/cold inlet
    
    rho_odes  = np.zeros(len(u_odes)) 
    for stream in StreamType:
        rho_odes[stream] = PropsSI('D', 'P', p[stream], 'T', T[stream], media_name[stream])

    # param_odes = OffDesignParam(param_des=param_des, G = np.multiply(u_odes, rho_odes))
    param_odes = OffDesignParam(param_des=param_des, mdot = mdot_odes)

    test = TestPCHEMeshram(work_root, param_des=param_des, param_odes=param_odes)

    test.run(simulate=True)    

###
if __name__ == "__main__":
    main()