from OMPython import OMCSessionZMQ
from OMPython import OMCSession
from OMPython import ModelicaSystem
from datetime import datetime as dt
import inspect

import csv

def gen_result_dict():
    
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

    # some special varible in the model, specify them explicitly
    # result_dict['eta_total'] = []
    return result_dict   

def cal_row( vals = [], eval_str = 'x + y'):

    real_eval_str = eval_str.replace('x','op1').replace('y','op2').replace('z', 'op3')

    r = []
    
    for i in range(0, len(vals[0])):
        op1,op2 = vals[0][i], vals[1][i]

        if(len(vals) == 3):
            op3 = vals[2][i]
        r.append(eval(real_eval_str))

    return r

def gen_cal_map():

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

def update_cal_fields(result_dict):
    # fill calculated fields into the result_dict
    cal_map = gen_cal_map()
    # fill in the calculated value by python's eval
    for item in cal_map:
        vals = item['vals']

        data = result_dict

        if len(vals) == 2:
            result_dict[item['key']] = cal_row([data[vals[0]], data[vals[1]]], eval_str=item['eval_str'])
        elif len(vals) == 3:
            result_dict[item['key']] = cal_row([data[vals[0]], data[vals[1]], data[vals[2]]], eval_str=item['eval_str'])      

def save_to_file(result_dict):
# save the simulation result into files
    with open(dt.now().strftime("result_%Y_%m_%d_%H_%M_%S.csv") ,mode='w') as csv_file:
        field_names = ['name']

        writer = csv.writer(csv_file, delimiter=',', quotechar='"', lineterminator="\n")

        writer.writerow(field_names)

        for k,v in result_dict.items():
            buf = [k]
            buf.extend(v)
            writer.writerow(buf)

def main():
    
    model_path = r"D:\Workspace\Steps\Steps"

    import os
    from shutil import copyfile
    # print(os.curdir)
    # os.chdir("..")
    pwd = "build" # setup working directory
    if not os.path.exists(pwd):    
        os.mkdir("build")
    os.chdir("build")

    libs = ["libCoolProp.a", 'libCoolProp.dll']
    for lib in libs:            
        if not os.path.exists(lib):
            copyfile(model_path + r"\Resources\Library\\" + lib, ".\\" + lib) # completa target name needed

    result_dict = gen_result_dict()

    omc = OMCSessionZMQ()
    path = inspect.getfile(omc.__class__)
    model_path = r"D:\Workspace\Steps\Steps"
    mod = ModelicaSystem(model_path + r"\package.mo","Steps.Test.TestPCHXMeshram",["Modelica 3.2.1"])
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

    #update_cal_fields()

    save_to_file(result_dict) 

    print('all done!')

###
if __name__ == "__main__":
    main()