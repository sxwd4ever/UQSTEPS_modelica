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
  constant Real pi=Modelica.Constants.pi;
  // efficiency of main compressor, bypass_compressor and turbine
  parameter Real eta_main_compressor = 0.89;  
  parameter Real eta_bypass_compressor = 0.89;  
  parameter Real eta_turbine = 0.89;
  
  // **** Boundary Conditions as Start values for recuperators - start ****
  // Core INPUT parameters for boundary conditions **** 
  // following values are calculated by sscar for 10 MW power block( mdot_main = 125, T_amb = 35oC, split_ratio = 0.675)
  parameter Modelica.SIunits.Pressure p_pump_in = 9e6;
  parameter Modelica.SIunits.Pressure p_pump_out = 20e6;
  parameter Modelica.SIunits.Pressure p_ATM = 101325;
  parameter Modelica.SIunits.Pressure p_heater = 20e6;
  
  parameter Modelica.SIunits.MassFlowRate mdot_main = 125;
  parameter Modelica.SIunits.MassFlowRate mdot_pump = 84.375;
  parameter Modelica.SIunits.MassFlowRate mdot_heater = 55;
  parameter Modelica.SIunits.MassFlowRate mdot_cooler = 40;
  
  parameter Modelica.SIunits.Temperature T_amb = from_degC(35);
  
  parameter Modelica.SIunits.Temperature T_HTR_hot_in = from_degC(636.95057734); 
  parameter Modelica.SIunits.Temperature T_HTR_cold_out = from_degC(605.011512655);
  
  parameter Modelica.SIunits.Temperature T_HTR_hot_out = from_degC(167.14450349);
  parameter Modelica.SIunits.Temperature T_HTR_cold_in = from_degC(162.14458875);     
  
  parameter Modelica.SIunits.Temperature T_LTR_cold_in = from_degC(85.84344872);
  parameter Modelica.SIunits.Temperature T_LTR_hot_out = from_degC(90.84344686);  
  
  parameter Modelica.SIunits.Temperature T_heater_hot_in = from_degC(800);
  parameter Modelica.SIunits.Temperature T_heater_hot_out = from_degC(600);
  
  // fixed pre-defined condition
  parameter Modelica.SIunits.Temperature T_heater_cold_out = from_degC(710);
  
  parameter Modelica.SIunits.Temperature T_cooler_cold_out = from_degC(61.7479462700001);
  parameter Modelica.SIunits.Temperature T_cooler_cold_in = T_amb; 
  // fixed pre-defined condition
  parameter Modelica.SIunits.Temperature T_cooler_hot_out = from_degC(45);  
  
  // default size of heat exchanger
  parameter Modelica.SIunits.Radius r_i = 10e-3 "10mm tube's internal radius";
  parameter Modelica.SIunits.Radius r_t = r_i + 2e-3 "tube's external radius";   
  parameter Modelica.SIunits.Radius r_o = r_t + 100e-3 "radius of external side of heat exchanger"; 
  parameter Modelica.SIunits.Length L = 10 "10m";
 
  /*
  // following values are calculated by IPESpro for 5 MW power block
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
  parameter Modelica.SIunits.Temperature T_cooler_cold_in = T_amb; 
  // fixed pre-defined condition
  parameter Modelica.SIunits.Temperature T_cooler_hot_out = from_degC(33);
  */
  // parameters based on input parameters.   
  parameter Real splitter_split_ratio = mdot_bypass/mdot_main "mass split ratio of splitter"; 
  parameter Modelica.SIunits.MassFlowRate mdot_bypass = mdot_main - mdot_pump "mdot in by pass path"; 
  
  parameter Modelica.SIunits.Temperature T_LTR_hot_in =  T_HTR_hot_out; 
  parameter Modelica.SIunits.Temperature T_LTR_cold_out = T_HTR_cold_in;     
  parameter Modelica.SIunits.Temperature T_bypass_out =  T_HTR_cold_in;    
  
  // DO NOT change following parameters - CHANGE input paramters instead  
  parameter HEBoundaryCondition bc_HTR(
    st_hot_in(p = p_pump_in, T = T_HTR_hot_in, h = specificEnthalpy_CO2(bc_HTR.st_hot_in), mdot = mdot_main),    
    st_cold_in(p = p_pump_out, T = T_HTR_cold_in, h = specificEnthalpy_CO2(bc_HTR.st_cold_in), mdot = mdot_main),
    st_hot_out(p = p_pump_in, T = T_HTR_hot_out, h = specificEnthalpy_CO2(bc_HTR.st_hot_out), mdot = mdot_main),
    st_cold_out(p = p_pump_out, T = T_HTR_cold_out, h = specificEnthalpy_CO2(bc_HTR.st_cold_out), mdot = mdot_main));      
  
  // boundary condition for LTR test @ diff mdot
  parameter HEBoundaryCondition bc_LTR(
    st_hot_in(p = p_pump_in, T = T_LTR_hot_in, h = specificEnthalpy_CO2(bc_LTR.st_hot_in),mdot = mdot_main),    
    st_cold_in(p = p_pump_out, T = T_LTR_cold_in, h = specificEnthalpy_CO2(bc_LTR.st_cold_in),mdot = mdot_pump),
    st_hot_out(p = p_pump_in, T = T_LTR_hot_out, h = specificEnthalpy_CO2(bc_LTR.st_hot_out),mdot = mdot_main),
    st_cold_out(p = p_pump_out, T = T_LTR_cold_out, h = specificEnthalpy_CO2(bc_LTR.st_cold_out),mdot = mdot_pump)); 
 
  parameter HEBoundaryCondition bc_cooler(
    st_hot_in(p = bc_LTR.st_hot_out.p, T = bc_LTR.st_hot_out.T, h = specificEnthalpy_CO2(bc_cooler.st_hot_in), mdot = mdot_pump),  
    st_cold_in(p = p_ATM, T = T_cooler_cold_in, h = specificEnthalpy_Water(bc_cooler.st_cold_in), mdot = mdot_cooler),
    st_hot_out(p = bc_LTR.st_hot_out.p, T = T_cooler_hot_out, h = specificEnthalpy_CO2(bc_cooler.st_hot_out), mdot = mdot_pump),
    st_cold_out(p = p_ATM, T = T_cooler_cold_out, h = specificEnthalpy_Water(bc_cooler.st_cold_out),mdot = mdot_cooler));
    
  parameter HEBoundaryCondition bc_heater(
    st_hot_in(p = p_heater, T = T_heater_hot_in, h = specificEnthalpy_Sodium(bc_heater.st_hot_in), mdot = mdot_heater),  
    st_cold_in(p = bc_HTR.st_cold_out.p, T = bc_HTR.st_cold_out.T, h = specificEnthalpy_CO2(bc_heater.st_cold_in), mdot = mdot_main),
    st_hot_out(p = p_heater, T = T_heater_hot_out, h = specificEnthalpy_Sodium(bc_heater.st_hot_out), mdot = mdot_heater),
    st_cold_out(p = bc_HTR.st_cold_out.p, T = T_heater_cold_out, h = specificEnthalpy_CO2(bc_heater.st_cold_out), mdot = mdot_main));   
 
  parameter ThermoState st_bypass(p = p_pump_out, T = T_bypass_out, h = specificEnthalpy_CO2(st_bypass), mdot = mdot_bypass);
  // **** Boundary Conditions as Start/Nominal values for recuperators - end ****      
  
  parameter PCHEGeoParam geo_HTR(
    // pitch length, m
    pitch = 12e-3,
    // pitch angle
    phi = from_deg((180 - 108) /2),
    // length of pche, m
    L = 1000e-3, // 2860e-3,
    // Diameter of semi_circular, m
    d = 2e-3,
    // number of channels
    N_ch = integer(94e3),
    // number of segments
    N_seg = 50);
  
  // test configuration for hot/cold/tube side of heatexchanger
  
  // In following calculation, V, A_ex are account for single tube/channel, not for total
  // check the Calculation in ThemoPower.PowerPlants.HRSG.Components.HE to understand the meaning of 
  // exsurface_G/F and extSurfaceTub, *Vol     
  parameter Real r_PCHE_HTR = 1e-3;
  parameter Real L_PCHE_HTR = 1000e-3;
  parameter EntityConfig cfg_PCHE_HTR_cold(
    geo(
    V = r_PCHE_HTR^2 * pi * L_PCHE_HTR, // * N_ch_LTR, 
    A_ex = pi * r_PCHE_HTR * L_PCHE_HTR, // * N_ch_LTR, exchange surface between fluid-tube
    L = L_PCHE_HTR, 
    d = r_PCHE_HTR * 2, 
    N_seg = 7, 
    N_ch = 332449), // calculate by steps the python code
    thermo(gamma_he = 3.96538615e6/cfg_HTR_cold.geo.A_ex "200")
  );
  
  parameter EntityConfig cfg_PCHE_HTR_tube(
    geo(
    V = cfg_PCHE_HTR_cold.geo.V * 2, //r_t_HTR^2 * pi * L_HTR * N_ch_HTR - cfg_HTR_cold.geo.V,
    A_ex = cfg_PCHE_HTR_cold.geo.A_ex, //  * N_ch_HTR, // assume thickness of tube approximately 0
    L = cfg_PCHE_HTR_cold.geo.L, 
    d = cfg_PCHE_HTR_cold.geo.d, 
    N_seg = 6, 
    N_ch = cfg_PCHE_HTR_cold.geo.N_ch),
    thermo(rho_mcm = 100, lambda = 20)
  ); 
  
  parameter EntityConfig cfg_PCHE_HTR_hot(
    geo(
    V = cfg_PCHE_HTR_cold.geo.V, // r_o_HTR^2 * pi * L_HTR * N_ch_HTR - cfg_HTR_tube.geo.V , 
    A_ex = cfg_PCHE_HTR_cold.geo.A_ex, // * N_ch_HTR, 
    L = cfg_PCHE_HTR_cold.geo.L, 
    d = cfg_PCHE_HTR_cold.geo.d, 
    N_seg = 7, 
    N_ch = cfg_PCHE_HTR_cold.geo.N_ch),
    thermo(gamma_he = 3.96538615e6/cfg_HTR_hot.geo.A_ex "200")
  );    
  
  parameter PCHEGeoParam geo_LTR(
    // pitch length, m
    pitch = 12e-3,
    // pitch angle
    phi = from_deg((180 - 108) /2),
    // length of pche, m
    L = 3270e-3,
    // Diameter of semi_circular, m
    d = 2e-3,
    // number of channels
    N_ch = integer(125e3),
    // number of segments
    N_seg = 50);
      
  parameter Real r_PCHE_LTR = 1e-3;
  parameter Real L_PCHE_LTR = 1000e-3;
  parameter EntityConfig cfg_PCHE_LTR_cold(
    geo(
    V = r_PCHE_LTR^2 * pi * L_PCHE_LTR, // * N_ch_LTR, 
    A_ex = pi * r_PCHE_LTR * L_PCHE_LTR, // * N_ch_LTR, exchange surface between fluid-tube
    L = L_PCHE_LTR, 
    d = r_PCHE_LTR * 2, 
    N_seg = 7, 
    N_ch = 422585), // calculate by steps the python code
    thermo(gamma_he = 1.938761018e6/cfg_LTR_cold.geo.A_ex "200")
  );
  
  parameter EntityConfig cfg_PCHE_LTR_tube(
    geo(
    V = cfg_PCHE_LTR_cold.geo.V * 2, //r_t_LTR^2 * pi * L_LTR * N_ch_LTR - cfg_LTR_cold.geo.V,
    A_ex = cfg_PCHE_LTR_cold.geo.A_ex, //  * N_ch_LTR, // assume thickness of tube approximately 0
    L = cfg_PCHE_LTR_cold.geo.L, 
    d = cfg_PCHE_LTR_cold.geo.d, 
    N_seg = 6, 
    N_ch = cfg_PCHE_LTR_cold.geo.N_ch),
    thermo(rho_mcm = 100, lambda = 20)
  ); 
  
  parameter EntityConfig cfg_PCHE_LTR_hot(
    geo(
    V = cfg_PCHE_LTR_cold.geo.V, // r_o_LTR^2 * pi * L_LTR * N_ch_LTR - cfg_LTR_tube.geo.V , 
    A_ex = cfg_PCHE_LTR_cold.geo.A_ex, // * N_ch_LTR, 
    L = cfg_PCHE_LTR_cold.geo.L, 
    d = cfg_PCHE_LTR_cold.geo.d, 
    N_seg = 7, 
    N_ch = cfg_PCHE_LTR_cold.geo.N_ch),
    thermo(gamma_he = 1.938761018e6/cfg_LTR_cold.geo.A_ex "200")
  );  
        
  // HTR's's size of heat exchanger gas - gas
  // N_ch_HTR groups of fluid(cold, inner)-tube-gas(hot, outter) tubes 
  parameter Modelica.SIunits.Radius r_i_HTR = 5e-3 "mm tube's internal radius";
  parameter Modelica.SIunits.Radius r_t_HTR = 10e-3 "tube's external radius";   
  parameter Modelica.SIunits.Radius r_o_HTR = 40e-3 "radius of external side of heat exchanger"; 
  parameter Integer N_ch_HTR = 400;
  parameter Modelica.SIunits.Length L_HTR = 1 "m";  
  
  // In following calculation, V, A_ex are account for single tube/channel, not for total
  // check the Calculation in ThemoPower.PowerPlants.HRSG.Components.HE to understand the meaning of 
  // exsurface_G/F and extSurfaceTub, *Vol     
  parameter EntityConfig cfg_HTR_cold(
    geo(
    V = r_i_HTR^2 * pi * L_HTR, // * N_ch_LTR, 
    A_ex = 2 * pi * r_i_HTR * L_HTR, // * N_ch_LTR, exchange surface between fluid-tube
    L = L_HTR, 
    d = 2 * r_i_HTR, 
    N_seg = 7, 
    N_ch = N_ch_HTR),
    thermo(gamma_he = 3.96538615e6/cfg_HTR_cold.geo.A_ex "200")
  );
  
  parameter EntityConfig cfg_HTR_tube(
    geo(
    V = r_t_HTR^2 * pi * L_HTR - cfg_HTR_cold.geo.V, //r_t_HTR^2 * pi * L_HTR * N_ch_HTR - cfg_HTR_cold.geo.V,
    A_ex = 2 * pi * r_t_HTR * L_HTR, //  * N_ch_HTR, // assume thickness of tube approximately 0
    L = L_HTR, 
    d = 2 * r_t_HTR, 
    N_seg = 6, 
    N_ch = N_ch_HTR),
    thermo(rho_mcm = 7900 * 578.05, lambda = 20)
  ); 
  
  parameter EntityConfig cfg_HTR_hot(
    geo(
    V = r_o_HTR^2 * pi * L_HTR - cfg_HTR_tube.geo.V, // r_o_HTR^2 * pi * L_HTR * N_ch_HTR - cfg_HTR_tube.geo.V , 
    A_ex = 2 * pi * r_o_HTR * L_HTR, // * N_ch_HTR, 
    L = L_HTR, 
    d = 2 * r_o_HTR, 
    N_seg = 7, 
    N_ch = N_ch_HTR),
    thermo(gamma_he = 3.96538615e6/cfg_HTR_hot.geo.A_ex "200")
  );  
  
  // LTR's's size of heat exchanger gas - gas
  // N_ch_LTR groups of fluid(cold, inner)-tube-gas(hot, outter) tubes 
  parameter Modelica.SIunits.Radius r_i_LTR = 100e-3 "mm tube's internal radius";
  parameter Modelica.SIunits.Radius r_t_LTR = 110e-3 "tube's external radius";   
  parameter Modelica.SIunits.Radius r_o_LTR = 170e-3 "radius of external side of one group"; 
  parameter Integer N_ch_LTR = 400;
  parameter Modelica.SIunits.Length L_LTR = 1 "m";  
  // In following calculation, V, A_ex are account for single tube/channel, not for total
  // check the Calculation in ThemoPower.PowerPlants.HRSG.Components.HE to understand the meaning of 
  // exsurface_G/F and extSurfaceTub, *Vol 
  parameter EntityConfig cfg_LTR_cold(
    geo(
    V = r_i_LTR^2 * pi * L_LTR, // * N_ch_LTR, 
    A_ex = 2 * pi * r_i_LTR * L_LTR, // * N_ch_LTR, exchange surface between fluid-tube
    L = L_LTR, 
    d = 2 * r_i_LTR, 
    N_seg = 7, 
    N_ch = N_ch_LTR),
    thermo(gamma_he = 1.938761018e6/cfg_LTR_cold.geo.A_ex "200")
  );
  
  parameter EntityConfig cfg_LTR_tube(
    geo(
    V = r_t_LTR^2 * pi * L_LTR - cfg_LTR_cold.geo.V,  //r_t_LTR^2 * pi * L_LTR * N_ch_LTR - cfg_LTR_cold.geo.V, 
    A_ex = 2 * pi * r_t_LTR * L_LTR, // * N_ch_LTR, // assume thickness of tube approximately 0
    L = L_LTR, 
    d = 2 * r_t_LTR, 
    N_seg = 6, 
    N_ch = N_ch_LTR),
    thermo(rho_mcm = 7900 * 578.05, lambda = 20)
  );   
  
  parameter EntityConfig cfg_LTR_hot(
    geo(
    V = r_o_LTR^2 * pi * L_LTR - cfg_LTR_tube.geo.V, //r_o_LTR^2 * pi * L_LTR * N_ch_LTR - cfg_LTR_tube.geo.V, 
    A_ex = 2 * pi * r_o_LTR * L_LTR, // * N_ch_LTR, 
    L = L_LTR, 
    d = 2 * r_o_LTR, 
    N_seg = 7, 
    N_ch = N_ch_LTR),
    thermo(gamma_he = 1.938761018e6/cfg_LTR_hot.geo.A_ex "200")
  );
  
  parameter EntityConfig cfg_mixer(
    geo(V = 3, A_ex = 1),
    thermo(gamma_he = 0) //steady state
  );
     
  // heater's size of heat exchanger
  // N_ch groups of fluid(hot, inner)-tube vs 1 gas(cold, outter) tube  
  parameter Modelica.SIunits.Radius r_i_h = 20e-3 "mm tube's internal radius";
  parameter Modelica.SIunits.Radius r_t_h = 30e-3 "tube's external radius";   
  parameter Modelica.SIunits.Radius r_o_h = 500e-3 "radius of external side of heat exchanger"; 
  parameter Integer N_ch_h = 100;
  parameter Modelica.SIunits.Length L_h = 1 "m"; 
  
  // cfg for heater's hot/fluid side
  parameter EntityConfig cfg_heater_hot(
    geo(
    V = r_i_h^2 * pi * L_h, 
    A_ex = 2 * pi * r_i_h * L_h, 
    L = L_h, 
    d = 2 * r_i_h, 
    N_seg = 7, 
    N_ch = 1),
    thermo(gamma_he = 200 "4000")
  );
    
  // cfg for heater's cold/fluid side
  parameter EntityConfig cfg_heater_cold(
    geo(
    V = r_o_h^2 * pi * L_h, 
    A_ex = 2 * pi * r_o_h * L_h, 
    L = L_h, 
    d = 2 * r_o_h, 
    N_seg = 7, 
    N_ch = N_ch_h),
    thermo(gamma_he = 200 "4000")
  );
  
  // cfg for heater's tube
  parameter EntityConfig cfg_heater_tube(
    geo(
    V = r_t_h^2 * pi * L_h, 
    A_ex = 2 * pi * r_t_h * L_h, 
    L = L_h, 
    d = 2 * r_t_h, 
    N_seg = 6, 
    N_ch = N_ch_h),
    thermo(rho_mcm = 7900 * 578.05, lambda = 20)
  );   
  
  // cooler's size of heat exchanger  
  // N_ch groups of fluid(cold, inner)-tube vs 1 gas(hot, outter) tube  
  parameter Modelica.SIunits.Radius r_i_c = 10e-3 "mm, tube's internal radius(single tube)";
  parameter Modelica.SIunits.Radius r_t_c = 15e-3 "tube's external radius(single tube)";   
  parameter Modelica.SIunits.Radius r_o_c = 500e-3 "radius of external side of heat exchanger"; 
  parameter Integer N_ch_c = 100;
  parameter Modelica.SIunits.Length L_c = 1 "m length of the tube/cooler";  
    
  // cfg for heater's cold/fluid side
  parameter EntityConfig cfg_cooler_cold(
    geo(
    V = r_i_c^2 * pi * L_c, 
    A_ex = 2 * pi * r_i_c * L_c, 
    L = L_c, 
    d = 2 * r_o_c, 
    N_ch = N_ch_c, 
    N_seg = 7),
    thermo(gamma_he = 200)
  );
  
  // cfg for heater's tube
  parameter EntityConfig cfg_cooler_tube(
    geo(
    V =  r_t_c^2 * pi * L_c, 
    A_ex = 2 * pi * r_t_c * L_c, 
    L = L_c, 
    d = 2 * r_t_c, 
    N_ch = N_ch_c, 
    N_seg = 6),
    thermo(rho_mcm = 7900 * 578.05, lambda = 20)
  );   
  
  // cfg for cooler's hot/fluid side
  parameter EntityConfig cfg_cooler_hot(
    geo(
    V = r_o_c^2 * pi * L_c, 
    A_ex = 2 * pi * r_o_c * L_c, 
    L = L_c, 
    d = 2 * r_i_c,
    N_ch = N_ch_c, 
    N_seg = 7),
    thermo(gamma_he = 200 "4000")
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
  function specificEnthalpy
    input ThermoState state;
    output Modelica.SIunits.SpecificEnthalpy h;
    replaceable package Medium = Modelica.Media.Interfaces.PartialMedium;    
  algorithm
    h := Medium.specificEnthalpy(Medium.setState_pTX(state.p, state.T));
  end specificEnthalpy;
  
  function specificEnthalpy_CO2    
    extends specificEnthalpy(redeclare package Medium = Steps.Media.CO2);
  end specificEnthalpy_CO2;
  
  function specificEnthalpy_Sodium  
    extends specificEnthalpy(redeclare package Medium = SolarTherm.Media.Sodium.Sodium_pT);   
  end specificEnthalpy_Sodium;
  
  function specificEnthalpy_Water
    extends specificEnthalpy(redeclare package Medium = ThermoPower.Water.StandardWater);          
  end specificEnthalpy_Water;  
  
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
