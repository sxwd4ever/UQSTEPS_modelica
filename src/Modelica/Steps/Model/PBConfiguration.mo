within Steps.Model;

model PBConfiguration
  "Preset param set for Off-design power block test - DON'T Alter it directly"
  
  import Modelica.SIunits.Conversions.{from_degC, from_deg};
  import Modelica.SIunits.{Temperature, Pressure, SpecificEnthalpy};
  import Util = Utilities.Util;
  import Steps.Utilities.CoolProp.PropsSI; 
  import Steps.Components.{PCHEBoundaryCondition, PCHEGeoParam};
  
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
  
  parameter PCHEGeoParam geo_HTR(
    // pitch length, m
    pitch = 12e-3,
    // pitch angle
    phi = from_deg((180 - 108) /2),
    // length of pche, m
    length = 2860e-3,
    // Diameter of semi_circular, m
    d_c = 2e-3,
    // number of channels
    N_ch = integer(94e3),
    // number of segments
    N_seg = 50);
      
  parameter PCHEGeoParam geo_LTR(
    // pitch length, m
    pitch = 12e-3,
    // pitch angle
    phi = from_deg((180 - 108) /2),
    // length of pche, m
    length = 3270e-3,
    // Diameter of semi_circular, m
    d_c = 2e-3,
    // number of channels
    N_ch = integer(125e3),
    // number of segments
    N_seg = 50);
  
  // test configuration for hot/cold/tube side of heatexchanger
  
  parameter EntityConfig cfg_hot(
    geo(V = 10, A_ex = 1708.2, N_seg = 7),
    thermo(gamma_he = 200)
  );
  
  parameter EntityConfig cfg_cold(
    geo(V = 2.234, A_ex = 225.073, N_seg = 7),
    thermo(gamma_he = 200)
  );
  
  parameter EntityConfig cfg_tube(
    geo(V = 0.573, A_ex = 252.286, N_seg = 6),
    thermo(rho_mcm = 7900 * 578.05, lambda = 20)
  );  

  // **** Boundary Conditions as Start values for all components - start ****  
  parameter PCHEBoundaryCondition bc_HTR(
    st_hot_in(p = p_pump_in, T = from_degC(578.22), h = PropsSI("H", "P",  bc_HTR.st_hot_in.p, "T", bc_HTR.st_hot_in.T, PBMedia.mediumName), mdot = mdot_main),    
    st_cold_in(p = p_pump_out, T = T_HTR_cold_in, h = PropsSI("H", "P", bc_HTR.st_cold_in.p, "T", bc_HTR.st_cold_in.T, PBMedia.mediumName), mdot = mdot_main),
    st_hot_out(p = p_pump_in, T = T_HTR_hot_out, h = PropsSI("H", "P", bc_HTR.st_hot_out.p, "T", bc_HTR.st_hot_out.T, PBMedia.mediumName), mdot = mdot_main),
    st_cold_out(p = p_pump_out, T = from_degC(533.5), h = PropsSI("H", "P", bc_HTR.st_cold_out.p, "T", bc_HTR.st_cold_out.T, PBMedia.mediumName), mdot = mdot_main));      
  
  // boundary condition for LTR test @ diff mdot
  parameter PCHEBoundaryCondition bc_LTR(
    st_hot_in(p = p_pump_in, T = T_LTR_hot_in, h = PropsSI("H", "P", bc_LTR.st_hot_in.p, "T", bc_LTR.st_hot_in.T, PBMedia.mediumName), mdot = mdot_main),    
    st_cold_in(p = p_pump_out, T = from_degC(62.229), h = PropsSI("H", "P", bc_LTR.st_cold_in.p, "T", bc_LTR.st_cold_in.T, PBMedia.mediumName), mdot = mdot_pump),
    st_hot_out(p = p_pump_in, T = from_degC(67.229), h = PropsSI("H", "P", bc_LTR.st_hot_out.p, "T", bc_LTR.st_hot_out.T, PBMedia.mediumName), mdot = mdot_main),
    st_cold_out(p = p_pump_out, T = T_LTR_cold_out, h = PropsSI("H", "P", bc_LTR.st_cold_out.p, "T", bc_LTR.st_cold_out.T, PBMedia.mediumName), mdot = mdot_pump)); 
       
  parameter ThermoState bc_cooler_out(p = bc_LTR.st_hot_out.p, T = from_degC(33), h = PropsSI("H", "P", bc_cooler_out.p, "T", bc_cooler_out.T, PBMedia.mediumName), mdot = mdot_pump);    
  parameter ThermoState bc_heater_out(p = bc_HTR.st_cold_out.p, T = from_degC(700), h = PropsSI("H", "P", bc_heater_out.p, "T", bc_heater_out.T, PBMedia.mediumName), mdot = mdot_main);    
  parameter ThermoState bc_bypass(p = bc_HTR.st_cold_in.p, T = T_bypass_out, h = PropsSI("H", "P", bc_bypass.p, "T", bc_bypass.T, PBMedia.mediumName), mdot = mdot_bypass);
  // **** Boundary Conditions as Start values for all components - end ****      
  
  // default sim parameters, 'slow' in solution finding with high accuracy, faster for PB's convergence 
  parameter SimParam sim_param_def(err=5e-4, delta_T_init = 5, N_iter = 20, step_rel=0.13, log_level = 4);
  // default sim parameters with verbose output for debug
  parameter SimParam sim_param_def_v(err=5e-4, delta_T_init = 5, N_iter = 20, step_rel=0.13, log_level = 0);
  // fast sim parameters, 'fast' in solution finding with lower accuracy, slower for PB's convergence 
  parameter SimParam sim_param_fast(err=1e-2, delta_T_init = 5, N_iter = 10, step_rel=0.3, log_level = 1);
end PBConfiguration;
