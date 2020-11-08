within Steps.Model;

model PBConfiguration
  "Preset param set for Off-design power block test - DON'T Alter it directly"
  
  import Modelica.SIunits.Conversions.{from_degC, from_deg};
  import Modelica.SIunits.{Temperature, Pressure, SpecificEnthalpy};
  import Util = Utilities.Util;
  import Steps.Utilities.CoolProp.PropsSI; 
  import Steps.Components.PCHEGeoParam;
  /*
  replaceable package medium_cold = Steps.Media.SCO2;   
  replaceable package medium_hot = Steps.Media.SCO2;  
  replaceable package medium_heater = SolarTherm.Media.Sodium.Sodium_pT;
  replaceable package medium_cooler = Modelica.Media.Water.WaterIF97_pT;
  */
  // efficiency of main compressor, bypass_compressor and turbine
  parameter Real eta_main_compressor = 0.89;  
  parameter Real eta_bypass_compressor = 0.89;  
  parameter Real eta_turbine = 0.89;
  
  // **** Boundary Conditions as Start values for recuperators - start ****
  // Core INPUT parameters for boundary conditions **** 
  parameter Modelica.SIunits.Pressure p_pump_in = 8e6;
  parameter Modelica.SIunits.Pressure p_pump_out = 20e6;
  parameter Modelica.SIunits.Pressure p_ATM = 1.01e6;
  parameter Modelica.SIunits.Pressure p_heater = 20e6;
  
  parameter Modelica.SIunits.MassFlowRate mdot_main = 51.91;
  parameter Modelica.SIunits.MassFlowRate mdot_pump = 31.31;
  parameter Modelica.SIunits.MassFlowRate mdot_heater = 42.73;
  parameter Modelica.SIunits.MassFlowRate mdot_cooler = 335.2;
  
  parameter Modelica.SIunits.Temperature T_amb = from_degC(15);
  parameter Modelica.SIunits.Temperature T_HTR_hot_in = from_degC(578.22);
  parameter Modelica.SIunits.Temperature T_HTR_cold_out = from_degC(533.5);
  
  parameter Modelica.SIunits.Temperature T_HTR_hot_out = from_degC(156.45);
  parameter Modelica.SIunits.Temperature T_HTR_cold_in = from_degC(151.45); 
  
  parameter Modelica.SIunits.Temperature T_LTR_cold_in = from_degC(62.229);
  parameter Modelica.SIunits.Temperature T_LTR_hot_out = from_degC(67.229);  
  
  parameter Modelica.SIunits.Temperature T_heater_hot_out = from_degC(600);
  parameter Modelica.SIunits.Temperature T_heater_hot_in = from_degC(800);
  // fixed pre-defined condition
  parameter Modelica.SIunits.Temperature T_heater_cold_out = from_degC(700);
  
  parameter Modelica.SIunits.Temperature T_cooler_cold_out = from_degC(30);
  parameter Modelica.SIunits.Temperature T_cooler_cold_in = from_degC(15); 
  // fixed pre-defined condition
  parameter Modelica.SIunits.Temperature T_cooler_hot_out = from_degC(33);
  
  // parameters based on input parameters.   
  parameter Real splitter_split_ratio = mdot_bypass/mdot_main "mass split ratio of splitter"; 
  parameter Modelica.SIunits.MassFlowRate mdot_bypass = mdot_main - mdot_pump "mdot in by pass path"; 
  
  parameter Modelica.SIunits.Temperature T_LTR_hot_in =  T_HTR_hot_out; 
  parameter Modelica.SIunits.Temperature T_LTR_cold_out = T_HTR_cold_in;     
  parameter Modelica.SIunits.Temperature T_bypass_out =  T_HTR_cold_in;    
  
  // DO NOT change following parameters - CHANGE input paramters instead  
  /*
  parameter HEBoundaryCondition bc_HTR(
    //st_hot_in = newThermoState_pT_CO2(p = p_pump_in, T = T_HTR_hot_in, mdot = mdot_main),    
    st_cold_in = newThermoState_pT_CO2(p = p_pump_out, T = T_HTR_cold_in, mdot = mdot_main),
    st_hot_out = newThermoState_pT_CO2(p = p_pump_in, T = T_HTR_hot_out, mdot = mdot_main),
    st_cold_out = newThermoState_pT_CO2(p = p_pump_out, T = T_HTR_cold_out, mdot = mdot_main));      


       
  parameter HEBoundaryCondition bc_heater(
    st_hot_in(newThermoState_pT_Sodium(p = p_heater, T = T_heater_hot_in, mdot = mdot_heater)),  
    st_cold_in = newThermoState_pT_CO2(p = bc_HTR.st_cold_out.p, T = bc_HTR.st_cold_out.T, mdot = mdot_main),
    st_hot_out = newThermoState_pT_Sodium(p = p_heater, T = T_heater_hot_out, mdot = mdot_heater),
    st_cold_out = newThermoState_pT_CO2(p = bc_HTR.st_cold_out.p, T = T_heater_cold_out, mdot = mdot_main));   
  */
  parameter HEBoundaryCondition bc_HTR(
    st_hot_in(p = p_pump_in, T = T_HTR_hot_in, h = specificEnthalpy_pT_CO2(bc_HTR.st_hot_in), mdot = mdot_main),    
    st_cold_in(p = p_pump_out, T = T_HTR_cold_in, h = specificEnthalpy_pT_CO2(bc_HTR.st_cold_in), mdot = mdot_main),
    st_hot_out(p = p_pump_in, T = T_HTR_hot_out, h = specificEnthalpy_pT_CO2(bc_HTR.st_hot_out), mdot = mdot_main),
    st_cold_out(p = p_pump_out, T = T_HTR_cold_out, h = specificEnthalpy_pT_CO2(bc_HTR.st_cold_out), mdot = mdot_main));      
  
  // boundary condition for LTR test @ diff mdot
  parameter HEBoundaryCondition bc_LTR(
    st_hot_in(p = p_pump_in, T = T_LTR_hot_in, h = specificEnthalpy_pT_CO2(bc_LTR.st_hot_in),mdot = mdot_main),    
    st_cold_in(p = p_pump_out, T = T_LTR_cold_in, h = specificEnthalpy_pT_CO2(bc_LTR.st_cold_in),mdot = mdot_pump),
    st_hot_out(p = p_pump_in, T = T_LTR_hot_out, h = specificEnthalpy_pT_CO2(bc_LTR.st_hot_out),mdot = mdot_main),
    st_cold_out(p = p_pump_out, T = T_LTR_cold_out, h = specificEnthalpy_pT_CO2(bc_LTR.st_cold_out),mdot = mdot_pump)); 
 
  parameter HEBoundaryCondition bc_cooler(
    st_hot_in(p = bc_LTR.st_hot_out.p, T = bc_LTR.st_hot_out.T, h = specificEnthalpy_pT_CO2(bc_cooler.st_hot_in), mdot = mdot_cooler),  
    st_cold_in(p = p_ATM, T = T_cooler_cold_in, h = specificEnthalpy_pT_Water(bc_cooler.st_cold_in), mdot = mdot_pump),
    st_hot_out(p = bc_LTR.st_hot_out.p, T = T_cooler_hot_out, h = specificEnthalpy_pT_CO2(bc_cooler.st_hot_out), mdot = mdot_cooler),
    st_cold_out(p = p_ATM, T = T_cooler_cold_out, h = specificEnthalpy_pT_Water(bc_cooler.st_cold_out),mdot = mdot_pump));
    
  parameter HEBoundaryCondition bc_heater(
    st_hot_in(p = p_heater, T = T_heater_hot_in, h = specificEnthalpy_pT_Sodium(bc_heater.st_hot_in), mdot = mdot_heater),  
    st_cold_in(p = bc_HTR.st_cold_out.p, h = specificEnthalpy_pT_CO2(bc_heater.st_cold_in), T = bc_HTR.st_cold_out.T, mdot = mdot_main),
    st_hot_out(p = p_heater, T = T_heater_hot_out, h = specificEnthalpy_pT_Sodium(bc_heater.st_hot_out), mdot = mdot_heater),
    st_cold_out(p = bc_HTR.st_cold_out.p, T = T_heater_cold_out, h = specificEnthalpy_pT_CO2(bc_heater.st_cold_out), mdot = mdot_main));   
 
  parameter ThermoState st_bypass(p = p_pump_out, T = T_bypass_out,  h = specificEnthalpy_pT_CO2(st_bypass), mdot = mdot_bypass);
  // **** Boundary Conditions as Start/Nominal values for recuperators - end ****      
  
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
  
  parameter EntityConfig cfg_HTR_hot(
    geo(V = 10, A_ex = 1708.2, N_seg = 7),
    thermo(gamma_he = 200)
  );
  
  parameter EntityConfig cfg_HTR_cold(
    geo(V = 2.234, A_ex = 225.073, N_seg = 7),
    thermo(gamma_he = 200)
  );
  
  parameter EntityConfig cfg_HTR_tube(
    geo(V = 0.573, A_ex = 252.286, N_seg = 6),
    thermo(rho_mcm = 7900 * 578.05, lambda = 20)
  ); 
  
  parameter EntityConfig cfg_mixer(
    geo(V = 3, A_ex = 1),
    thermo(gamma_he = 0) //steady state
  );   
  /*
  parameter EntityConfig cfg_heater(
    geo(V = 3, A_ex = 1),
    thermo(rho_mcm = 7900 * 578.05, lambda = 20, gamma_he = 200) //steady state
  ); 
  */
  // cfg for heater's hot/fluid side
  parameter EntityConfig cfg_heater_hot(
    geo(V = 2.234, A_ex = 225.073, N_seg = 7),
    thermo(gamma_he = 200 "4000")
  );
  
  // cfg for heater's cold/fluid side
  parameter EntityConfig cfg_heater_cold(
    geo(V = 1 "10", A_ex = 225.073 "1708.2", N_seg = 7),
    thermo(gamma_he = 200)
  );
  
  // cfg for heater's tube
  parameter EntityConfig cfg_heater_tube(
    geo(V = 0.573, A_ex = 252.286, N_seg = 6),
    thermo(rho_mcm = 7900 * 578.05, lambda = 20)
  );   
  
  // default sim parameters, 'slow' in solution finding with high accuracy, faster for PB's convergence 
  parameter SimParam sim_param_def(err=5e-4, delta_T_init = 5, N_iter = 20, step_rel=0.13, log_level = 4);
  // default sim parameters with verbose output for debug
  parameter SimParam sim_param_def_v(err=5e-4, delta_T_init = 5, N_iter = 20, step_rel=0.13, log_level = 0);
  // fast sim parameters, 'fast' in solution finding with lower accuracy, slower for PB's convergence 
  parameter SimParam sim_param_fast(err=1e-2, delta_T_init = 5, N_iter = 10, step_rel=0.3, log_level = 1);
  
  
  /** 
    * factory method to complete thermostate - fill in specific enthalpy h as determined by (p, T)
    * I can not replace the package for a function. 
    * So I use three inheirated functions to specify the medium explicitly
    */    
  function specificEnthalpy_pT
    input ThermoState state;
    output Modelica.SIunits.SpecificEnthalpy h;
    replaceable package Medium = Modelica.Media.Interfaces.PartialMedium;    
  algorithm
    h := Medium.specificEnthalpy(Medium.setState_pTX(state.p, state.T));
  end specificEnthalpy_pT;
  
  function specificEnthalpy_pT_CO2    
    extends specificEnthalpy_pT(redeclare package Medium = Steps.Media.CO2);
  end specificEnthalpy_pT_CO2;
  
  function specificEnthalpy_pT_Sodium  
    extends specificEnthalpy_pT(redeclare package Medium = SolarTherm.Media.Sodium.Sodium_pT);   
  end specificEnthalpy_pT_Sodium;
  
  function specificEnthalpy_pT_Water
    extends specificEnthalpy_pT(redeclare package Medium = ThermoPower.Water.StandardWater);          
  end specificEnthalpy_pT_Water;  
/*
  function completeThermoState_pT_CO2    
    extends completeThermoState_pT(redeclare package Medium = Steps.Media.CO2);
  end completeThermoState_pT_CO2;
  
  function completeThermoState_pT_Sodium  
    extends completeThermoState_pT(redeclare package Medium = SolarTherm.Media.Sodium.Sodium_pT);   
  end completeThermoState_pT_Sodium;
  
  function completeThermoState_pT_Water
    extends completeThermoState_pT(redeclare package Medium = ThermoPower.Water.StandardWater);          
  end completeThermoState_pT_Water;  
*/  
/*
  equation
  
    completeThermoState_pT_CO2(bc_HTR.st_hot_in);
    completeThermoState_pT_CO2(bc_HTR.st_cold_in);
    completeThermoState_pT_CO2(bc_HTR.st_hot_out);
    completeThermoState_pT_CO2(bc_HTR.st_cold_out);
    
    completeThermoState_pT_Sodium(bc_heater.st_hot_in);
    completeThermoState_pT_CO2(bc_heater.st_cold_in);
    completeThermoState_pT_Sodium(bc_heater.st_hot_out);
    completeThermoState_pT_CO2(bc_heater.st_cold_out);
   
    completeThermoState_pT_CO2(st_bypass);*/ 
    /*
    bc_HTR.st_hot_in := newThermoState_pT_CO2(p = p_pump_in, T = T_HTR_hot_in, mdot = mdot_main);    
    bc_HTR.st_cold_in := newThermoState_pT_CO2(p = p_pump_out, T = T_HTR_cold_in, mdot = mdot_main);
    bc_HTR.st_hot_out := newThermoState_pT_CO2(p = p_pump_in, T = T_HTR_hot_out, mdot = mdot_main);
    bc_HTR.st_cold_out := newThermoState_pT_CO2(p = p_pump_out, T = T_HTR_cold_out, mdot = mdot_main); 
    
    
    bc_heater.st_hot_in := newThermoState_pT_Sodium(p = p_heater, T = T_heater_hot_in, mdot = mdot_heater);
    bc_heater.st_cold_in := newThermoState_pT_CO2(p = bc_HTR.st_cold_out.p, T = bc_HTR.st_cold_out.T, mdot = mdot_main);
    bc_heater.st_hot_out := newThermoState_pT_Sodium(p = p_heater, T = T_heater_hot_out, mdot = mdot_heater);
    bc_heater.st_cold_out := newThermoState_pT_CO2(p = bc_HTR.st_cold_out.p, T = T_heater_cold_out, mdot = mdot_main);
    */
end PBConfiguration;
