within Steps.Model;

model DefaultPBConfiguration
  "Preset param set for Off-design power block test - DON'T Alter it directly"
  
  import Modelica.SIunits.Conversions.{from_degC, from_deg};
  import Modelica.SIunits.{Temperature, Pressure, SpecificEnthalpy};
  import Util = Utilities.Util;
  import Steps.Utilities.CoolProp.PropsSI; 
  import Steps.Components.PCHEBoundaryCondition;

  replaceable package PBMedia = Steps.Media.SCO2;   

  // efficiency of main compressor, bypass_compressor and turbine
  parameter Real eta_main_compressor = 0.89;
  
  parameter Real eta_bypass_compressor = 0.89;
  
  parameter Real eta_turbine = 0.89;
  
  // mass split ratio of splitter
  parameter Real splitter_split_ratio = mdot_bypass/mdot_main; 

  parameter Modelica.SIunits.Pressure p_pump_in = 8e6;
  parameter Modelica.SIunits.Pressure p_pump_out = 20e6;
  
  parameter Modelica.SIunits.MassFlowRate mdot_main = 51.91;
  parameter Modelica.SIunits.MassFlowRate mdot_pump = 31.31;
  parameter Modelica.SIunits.MassFlowRate mdot_bypass = mdot_main - mdot_pump;  
  
  parameter Modelica.SIunits.Temperature T_HTR_hot_out = from_degC(156.45);
  parameter Modelica.SIunits.Temperature T_HTR_cold_in = from_degC(151.45);
  
  parameter Modelica.SIunits.Temperature T_LTR_cold_out = from_degC(151.45);
  parameter Modelica.SIunits.Temperature T_LTR_hot_in =  from_degC(156.45);   
  
  parameter Modelica.SIunits.Temperature T_bypass_out =  from_degC(151.45);   
  
  // default configuration for gas 2 gas heat exchanger
  HEConfig cfg_heg2g;
    
   // **** Boundary Conditions as Start values for all components - end ****      
   
   // default sim parameters, 'slow' in solution finding with high accuracy, faster for PB's convergence 
   parameter SimParam sim_param_def(err=5e-4, delta_T_init = 5, N_iter = 20, step_rel=0.13, log_level = 4);
   // default sim parameters with verbose output for debug
   parameter SimParam sim_param_def_v(err=5e-4, delta_T_init = 5, N_iter = 20, step_rel=0.13, log_level = 0);
   // fast sim parameters, 'fast' in solution finding with lower accuracy, slower for PB's convergence 
   parameter SimParam sim_param_fast(err=1e-2, delta_T_init = 5, N_iter = 10, step_rel=0.3, log_level = 1);
end DefaultPBConfiguration;
