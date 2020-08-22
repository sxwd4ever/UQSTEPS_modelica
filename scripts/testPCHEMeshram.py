from OMPython import OMCSessionZMQ
from OMPython import OMCSession
from OMPython import ModelicaSystem
from datetime import datetime as dt

import inspect
import csv
import os
import numpy as np

class TestPCHEMeshram(object):
        
    def __init__(self, work_root):
        super().__init__()
        # common variable defination
        self.path_pics = "pics"
        self.path_out = "out"

        self.work_root = work_root
        self.model_path_root = work_root + r"/Steps" 
        
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
        var_keys = ['T', 'dp','Re', 'Nu', 'f']
        node_keys = ['cell_cold', 'cell_hot']
        for var_key in var_keys:
            for node_key in node_keys:
                for i in range(1, N_seg + 1):        
                    result_dict["pchx.{node}[{idx}].{var}".format(idx = i, var = var_key, node = node_key)] = []

        # some special varible in the model, specify them explicitly like
        # result_dict['eta_total'] = []  

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

    def run_simulation(self, result_dict):

        omc = OMCSessionZMQ()
        # path = inspect.getfile(omc.__class__)

        mod = ModelicaSystem(self.model_path_root + r"\package.mo","Steps.Test.TestPCHXMeshram",["Modelica 3.2.1"])
        # options for simulation - steady state simulation, no iteration required, so set numberOfIntervals = 2 
        mod.setSimulationOptions('stepSize  = 0.2')

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
                result_dict[item['key']] = cal_row([data[vals[0]], data[vals[1]]], eval_str=item['eval_str'])
            elif len(vals) == 3:
                result_dict[item['key']] = cal_row([data[vals[0]], data[vals[1]], data[vals[2]]], eval_str=item['eval_str'])      

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
                writer.writerow({'name': k, 'value': v})                                

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
                values = [float(x) for x in value_str[1:-1].split(',')]

                result_dict[line['name']] = values[0]

        return result_dict

    def draw_plot(self, x_values = [], y_values = []):
        '''
        Draw the result on figures of Meshram 2016 to compare
        '''
        from plot_lib import imgplot, plotanno, plot_xy_series

        zigzag = 0
        axisx=[[4000,26000],[4000,32000]][zigzag]
        axisT=[[400.0,750.0],[400, 750.0]][zigzag]
        axisz=[[0, 0.12], [0, 0.16]][zigzag]
        (Nu_plot, f_plot, T_plot, dp_plot)=(False, False, True, False) 

        imgfile = [self.path_pics + "/Meshram_Fig_04.png", self.path_pics + "/Meshram_Fig_05.jpg"][zigzag]            
        zticks=[["0.", "0.06", "0.12"],["0.", "0.08", "0.16"]][zigzag]
        Tticks=[["400", "575", "750"],["400", "575", "750"]][zigzag]
        dpyticks=[["0","2.50", "5.0"],["0", "40","80"]][zigzag]     
        if (T_plot):
            axis=[axisz,axisT]
            
            (fig, ax) = plot_xy_series(x=x_values, y=y_values, range_x=axisz, range_y=axisT, img_file=imgfile)#xticks=zticks, yticks=Tticks) 

            # (fig, ax)=imgplot(x_values, y_values, axis, "r", imgfile=imgfile, xticks=zticks, yticks=Tticks)
                 
            # for cold stream
            # imgplot(z, x.Tc+273, axis, "b--", ax=ax)
            plotanno(ax,xlabel="Z(m)", ylabel="T(oK)", title="Higher Temperature Range")
            fig.savefig( self.path_out + "/Meshram_Fig%db_T_compare.png"%(4+zigzag))
        if (dp_plot):
            axis=[axisz,axisp]
            z=np.linspace(0, x.length(), len(x.dpH))

    def run(self, simulate = True):

        if simulate:
            result_dict = self.__gen_result_dict()              

            self.run_simulation(result_dict)

            self.save_results(result_dict)
        
        # load the simulation result from latest, pre-saved file
        result_dict = self.load_result()

        #update_cal_fields()        

        # generate T series

        N_seg = 10
        len_seg = 12e-3
        T_hot = np.array([ float(i) for i in range(1, N_seg + 1)])
        x_values = np.arange(0, len_seg * (N_seg) , len_seg)        

        for i in range(1, N_seg + 1):
            T_hot[i - 1] = result_dict['pchx.cell_cold[{0}].T'.format(i)]

        self.draw_plot(x_values = x_values, y_values = T_hot)

        print('all done!') 


def main(work_root = []):
    # root path of modelica root
    if work_root == []:
        work_root = os.path.abspath(os.curdir)

    test = TestPCHEMeshram(work_root)

    test.run(simulate= False)    

###
if __name__ == "__main__":
    main()