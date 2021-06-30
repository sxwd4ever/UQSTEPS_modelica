within Steps.Model;

model PBConfiguration
  "full set of params set for Off-design RCBC - DON'T Alter it directly, use override instead, see Model PBConfigs"
  
  import Modelica.SIunits.Conversions.{from_degC, from_deg};
  import Modelica.SIunits.{Temperature, Pressure, SpecificEnthalpy};
  import Util = Utilities.Util;
  import Steps.Utilities.CoolProp.PropsSI; 
  import Steps.Components.PCHEGeoParam;  

  replaceable package medium_HTR_cold = Steps.Media.SCO2 constrainedby Modelica.Media.Interfaces.PartialPureSubstance;
  replaceable package medium_HTR_hot  = Steps.Media.SCO2 constrainedby Modelica.Media.Interfaces.PartialPureSubstance;
  
  replaceable package medium_LTR_cold = Steps.Media.SCO2 constrainedby Modelica.Media.Interfaces.PartialPureSubstance;
  replaceable package medium_LTR_hot  = Steps.Media.SCO2 constrainedby Modelica.Media.Interfaces.PartialPureSubstance;

  replaceable package medium_cooler_cold = Modelica.Media.Water.WaterIF97_pT constrainedby Modelica.Media.Interfaces.PartialPureSubstance;
  replaceable package medium_cooler_hot  = Steps.Media.SCO2 constrainedby Modelica.Media.Interfaces.PartialPureSubstance;
  
  replaceable package medium_heater_hot  = Steps.Media.MoltenSalt.MoltenSalt_pT constrainedby Modelica.Media.Interfaces.PartialPureSubstance;
  replaceable package medium_heater_cold = Steps.Media.SCO2 constrainedby Modelica.Media.Interfaces.PartialPureSubstance;
  // replaceable package medium_heater_cold = Modelica
  //replaceable package medium_cooler = Modelica.Media.Water.WaterIF97_pT;
  
  constant Real pi = Modelica.Constants.pi;
  // efficiency of main compressor, bypass_compressor and turbine
  parameter Real eta_main_compressor   = 0.89;
  parameter Real eta_bypass_compressor = 0.89;
  parameter Real eta_turbine           = 0.89;
  
  // **** Boundary Conditions as Start values for recuperators - start ****
  // Core INPUT parameters for boundary conditions **** 
  // following values are calculated by sscar for 10 MW power block( mdot_main = 125, T_amb = 35oC, split_ratio = 0.675)
  parameter Modelica.SIunits.Pressure p_pump_in  = 9e6;
  parameter Modelica.SIunits.Pressure p_pump_out = 20e6;
  parameter Modelica.SIunits.Pressure p_ATM      = 101325;
  parameter Modelica.SIunits.Pressure p_heater   = 20e6;
  
  parameter Modelica.SIunits.MassFlowRate mdot_main   = 125;
  parameter Modelica.SIunits.MassFlowRate mdot_pump   = 84.375;
  parameter Modelica.SIunits.MassFlowRate mdot_heater = 55;
  parameter Modelica.SIunits.MassFlowRate mdot_cooler = 40;
  
  parameter Modelica.SIunits.Temperature T_amb = from_degC(35);
  
  parameter Modelica.SIunits.Temperature T_HTR_hot_in   = from_degC(636.95057734);
  parameter Modelica.SIunits.Temperature T_HTR_cold_out = from_degC(605.011512655);
  
  parameter Modelica.SIunits.Temperature T_HTR_hot_out = from_degC(167.14450349);
  parameter Modelica.SIunits.Temperature T_HTR_cold_in = from_degC(162.14458875);
  
  parameter Modelica.SIunits.Temperature T_LTR_cold_in = from_degC(85.84344872);
  parameter Modelica.SIunits.Temperature T_LTR_hot_out = from_degC(90.84344686);
  
  parameter Modelica.SIunits.Temperature T_heater_hot_in  = from_degC(800);
  parameter Modelica.SIunits.Temperature T_heater_hot_out = from_degC(600);
  
  // fixed pre-defined condition
  parameter Modelica.SIunits.Temperature T_heater_cold_out = from_degC(710);
  
  parameter Modelica.SIunits.Temperature T_cooler_cold_out = from_degC(61.7479462700001);
  parameter Modelica.SIunits.Temperature T_cooler_cold_in  = T_amb;
  // fixed pre-defined condition
  parameter Modelica.SIunits.Temperature T_cooler_hot_out = from_degC(45);
  
  parameter Integer N_seg = 10 "default number of discretized segments in one tube";
  
  /*
  // following values are calculated by IPESpro for 5 MW power block

  */
  // parameters based on input parameters.   
  parameter Real splitter_split_ratio                 = mdot_bypass/mdot_main "mass split ratio of splitter";
  parameter Modelica.SIunits.MassFlowRate mdot_bypass = mdot_main - mdot_pump "mdot in by pass path";
  
  parameter Modelica.SIunits.Temperature T_LTR_hot_in   = T_HTR_hot_out;
  parameter Modelica.SIunits.Temperature T_LTR_cold_out = T_HTR_cold_in;
  parameter Modelica.SIunits.Temperature T_bypass_out   = T_HTR_cold_in;
  
  // DO NOT change following parameters - CHANGE input paramters instead  
  parameter HEBoundaryCondition bc_HTR(
    st_hot_in(
    p    = p_pump_in,
    T    = T_HTR_hot_in,
    h    = medium_HTR_hot.specificEnthalpy(medium_HTR_hot.setState_pT(p = p_pump_in, T = T_HTR_hot_in)),
    mdot = mdot_main),
    st_cold_in(
    p    = p_pump_out,
    T    = T_HTR_cold_in,
    h    = medium_HTR_cold.specificEnthalpy(medium_HTR_cold.setState_pT(p = p_pump_out, T = T_HTR_cold_in)),
    mdot = mdot_pump),
    st_hot_out(
    p    = p_pump_in,
    T    = T_HTR_hot_out,
    h    = medium_HTR_hot.specificEnthalpy(medium_HTR_hot.setState_pT(p = p_pump_in, T = T_HTR_hot_out)),
    mdot = mdot_main),
    st_cold_out(
    p    = p_pump_out,
    T    = T_HTR_cold_out,
    h    = medium_HTR_cold.specificEnthalpy(medium_HTR_cold.setState_pT(p = p_pump_out, T = T_HTR_cold_out)),
    mdot = mdot_pump));

  
  // boundary condition for LTR test @ diff mdot
  parameter HEBoundaryCondition bc_LTR(
    st_hot_in(
    p    = p_pump_in,
    T    = T_LTR_hot_in,
    h    = medium_LTR_hot.specificEnthalpy(medium_LTR_hot.setState_pT(p = p_pump_in, T = T_LTR_hot_in)),
    mdot = mdot_main),
    st_cold_in(
    p    = p_pump_out,
    T    = T_LTR_cold_in,
    h    = medium_LTR_cold.specificEnthalpy(medium_LTR_cold.setState_pT(p = p_pump_out, T = T_LTR_cold_in)),
    mdot = mdot_pump),
    st_hot_out(
    p    = p_pump_in,
    T    = T_LTR_hot_out,
    h    = medium_LTR_hot.specificEnthalpy(medium_LTR_hot.setState_pT(p = p_pump_in, T = T_LTR_hot_out)),
    mdot = mdot_main),
    st_cold_out(
    p    = p_pump_out,
    T    = T_LTR_cold_out,
    h    = medium_LTR_cold.specificEnthalpy(medium_LTR_cold.setState_pT(p = p_pump_out, T = T_LTR_cold_out)),
    mdot = mdot_pump));
 
  parameter HEBoundaryCondition bc_cooler(
    st_hot_in(
    p    = bc_LTR.st_hot_out.p,
    T    = bc_LTR.st_hot_out.T,
    h    = medium_cooler_hot.specificEnthalpy(medium_cooler_hot.setState_pT(p = bc_LTR.st_hot_out.p, T = bc_LTR.st_hot_out.T)),
    mdot = mdot_pump),
    st_cold_in(
    p    = p_ATM,
    T    = T_cooler_cold_in,
    h    = medium_cooler_cold.specificEnthalpy(medium_cooler_cold.setState_pT(p = p_ATM, T = T_cooler_cold_in)),
    mdot = mdot_cooler),
    st_hot_out(
    p    = bc_LTR.st_hot_out.p,
    T    = T_cooler_hot_out,
    h    = medium_cooler_hot.specificEnthalpy(medium_cooler_hot.setState_pT(p = bc_LTR.st_hot_out.p, T = T_cooler_hot_out)),
    mdot = mdot_pump),
    st_cold_out(
    p    = p_ATM,
    T    = T_cooler_cold_out,
    h    = medium_cooler_cold.specificEnthalpy(medium_cooler_cold.setState_pT(p = p_ATM, T = T_cooler_cold_out)),
    mdot = mdot_cooler));
    
  parameter HEBoundaryCondition bc_heater(
    st_hot_in(
    p    = p_heater,
    T    = T_heater_hot_in,
    h    = medium_heater_hot.specificEnthalpy(medium_heater_hot.setState_pT(p_heater, T_heater_hot_in)),
    mdot = mdot_heater),
    st_cold_in(
    p    = bc_HTR.st_cold_out.p,
    T    = bc_HTR.st_cold_out.T,
    h    = medium_heater_cold.specificEnthalpy(medium_heater_cold.setState_pT(p = bc_HTR.st_cold_out.p, T = bc_HTR.st_cold_out.T)),
    mdot = mdot_main),
    st_hot_out(
    p    = p_heater,
    T    = T_heater_hot_out,
    h    = medium_heater_hot.specificEnthalpy(medium_heater_hot.setState_pT(p_heater, T_heater_hot_out)),
    mdot = mdot_heater),
    st_cold_out(
    p    = bc_HTR.st_cold_out.p,
    T    = T_heater_cold_out,
    h    = medium_heater_cold.specificEnthalpy(medium_heater_cold.setState_pT(p = bc_HTR.st_cold_out.p, T = T_heater_cold_out)),
    mdot = mdot_main));
 
  parameter ThermoState st_bypass(p = p_pump_out, T = T_bypass_out, h = medium_HTR_cold.specificEnthalpy(medium_HTR_cold.setState_pT(p = p_pump_out, T = T_bypass_out)), mdot = mdot_bypass);
  // **** Boundary Conditions as Start/Nominal values for recuperators - end ****      
  
  
  // HTR's's size of heat exchanger gas - gas
  // N_ch_HTR groups of fluid(cold, inner)-tube-gas(hot, outter) tubes 
  parameter Modelica.SIunits.Radius r_i_HTR = 5e-3 "mm tube's internal radius";
  parameter Modelica.SIunits.Radius r_t_HTR = 10e-3 "tube's external radius";
  parameter Modelica.SIunits.Radius r_o_HTR = 40e-3 "radius of external side of heat exchanger";
  parameter Modelica.SIunits.Length L_HTR   = 1 "m";
  parameter Integer N_ch_HTR                = 400;
  parameter Integer N_seg_HTR               = N_seg;
  
  // In following calculation, V, A_ex are account for single tube/channel, not for total
  // check the Calculation in ThemoPower.PowerPlants.HRSG.Components.HE to understand the meaning of 
  // exsurface_G/F and extSurfaceTub, *Vol     
  parameter EntityConfig cfg_HTR_cold(
    geo(
    V     = r_i_HTR^2 * pi * L_HTR,     // * N_ch_LTR, 
    A_ex  = 2 * pi * r_i_HTR * L_HTR,   // * N_ch_LTR, exchange surface between fluid-tube
    L     = L_HTR,
    d     = 2 * r_i_HTR,
    N_seg = N_seg_HTR,
    N_ch  = N_ch_HTR),
    thermo(UAnom = 3.96538615e6, gamma_he = cfg_HTR_cold.thermo.UAnom/(cfg_HTR_cold.geo.A_ex * N_ch_HTR)  "200")
  );
  
  parameter EntityConfig cfg_HTR_tube(
    geo(
    V     = r_t_HTR^2 * pi * L_HTR - cfg_HTR_cold.geo.V,   //r_t_HTR^2 * pi * L_HTR * N_ch_HTR - cfg_HTR_cold.geo.V,
    A_ex  = 2 * pi * r_t_HTR * L_HTR,                      //  * N_ch_HTR, // assume thickness of tube approximately 0
    L     = L_HTR,
    d     = 2 * r_t_HTR,
    N_seg = N_seg_HTR,
    N_ch  = N_ch_HTR),
    thermo(rho_mcm = 7900 * 578.05, lambda = 20)
  ); 
  
  parameter EntityConfig cfg_HTR_hot(
    geo(
    V     = r_o_HTR^2 * pi * L_HTR - cfg_HTR_tube.geo.V,   // r_o_HTR^2 * pi * L_HTR * N_ch_HTR - cfg_HTR_tube.geo.V , 
    A_ex  = 2 * pi * r_o_HTR * L_HTR,                      // * N_ch_HTR, 
    L     = L_HTR,
    d     = 2 * r_o_HTR,
    N_seg = N_seg_HTR,
    N_ch  = N_ch_HTR),
    thermo(UAnom = 3.96538615e6, gamma_he = cfg_HTR_hot.thermo.UAnom/(cfg_HTR_hot.geo.A_ex * N_ch_HTR) "200")
  );  
  
  // LTR's's size of heat exchanger gas - gas
  // N_ch_LTR groups of fluid(cold, inner)-tube-gas(hot, outter) tubes 
  parameter Modelica.SIunits.Radius r_i_LTR = 100e-3 "mm tube's internal radius";
  parameter Modelica.SIunits.Radius r_t_LTR = 110e-3 "tube's external radius";
  parameter Modelica.SIunits.Radius r_o_LTR = 170e-3 "radius of external side of one group";
  parameter Modelica.SIunits.Length L_LTR   = 1 "m";
  parameter Integer N_seg_LTR               = N_seg;
  parameter Integer N_ch_LTR                = 400;
  // In following calculation, V, A_ex are account for single tube/channel, not for total
  // check the Calculation in ThemoPower.PowerPlants.HRSG.Components.HE to understand the meaning of 
  // exsurface_G/F and extSurfaceTub, *Vol 
  parameter EntityConfig cfg_LTR_cold(
    geo(
    V     = r_i_LTR^2 * pi * L_LTR,     // * N_ch_LTR, 
    A_ex  = 2 * pi * r_i_LTR * L_LTR,   // * N_ch_LTR, exchange surface between fluid-tube
    L     = L_LTR,
    d     = 2 * r_i_LTR,
    N_seg = N_seg_LTR,
    N_ch  = N_ch_LTR),
    thermo(UAnom = 1.938761018e6, gamma_he = cfg_LTR_cold.thermo.UAnom/(cfg_LTR_cold.geo.A_ex * N_ch_LTR)  "200")
  );
  
  parameter EntityConfig cfg_LTR_tube(
    geo(
    V     = r_t_LTR^2 * pi * L_LTR - cfg_LTR_cold.geo.V,   //r_t_LTR^2 * pi * L_LTR * N_ch_LTR - cfg_LTR_cold.geo.V, 
    A_ex  = 2 * pi * r_t_LTR * L_LTR,                      // * N_ch_LTR, // assume thickness of tube approximately 0
    L     = L_LTR,
    d     = 2 * r_t_LTR,
    N_seg = N_seg_LTR,
    N_ch  = N_ch_LTR),
    thermo(rho_mcm = 7900 * 578.05, lambda = 20)
  );   
  
  parameter EntityConfig cfg_LTR_hot(
    geo(
    V     = r_o_LTR^2 * pi * L_LTR - cfg_LTR_tube.geo.V,   //r_o_LTR^2 * pi * L_LTR * N_ch_LTR - cfg_LTR_tube.geo.V, 
    A_ex  = 2 * pi * r_o_LTR * L_LTR,                      // * N_ch_LTR, 
    L     = L_LTR,
    d     = 2 * r_o_LTR,
    N_seg = N_seg_LTR,
    N_ch  = N_ch_LTR),
    thermo(UAnom = 1.938761018e6, gamma_he = cfg_LTR_hot.thermo.UAnom/(cfg_LTR_hot.geo.A_ex * N_ch_LTR) "200")
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
  parameter Modelica.SIunits.Length L_h   = 1 "m";
  parameter Integer N_ch_h                = 100;
  parameter Integer N_seg_heater          = N_seg;
  
  // cfg for heater's hot/fluid side
  parameter EntityConfig cfg_heater_hot(
    geo(
    V     = r_i_h^2 * pi * L_h,
    A_ex  = 2 * pi * r_i_h * L_h,
    L     = L_h,
    d     = 2 * r_i_h,
    N_seg = N_seg_heater,
    N_ch  = 1),
    thermo(UAnom = 1.938761018e6, gamma_he = cfg_heater_hot.thermo.UAnom/(cfg_heater_hot.geo.A_ex * N_ch_h) "200 4000")
  );
    
  // cfg for heater's cold/fluid side
  parameter EntityConfig cfg_heater_cold(
    geo(
    V     = r_o_h^2 * pi * L_h,
    A_ex  = 2 * pi * r_o_h * L_h,
    L     = L_h,
    d     = 2 * r_o_h,
    N_seg = N_seg_heater,
    N_ch  = N_ch_h),
    thermo(UAnom = 1.938761018e6, gamma_he = cfg_heater_cold.thermo.UAnom/(cfg_heater_cold.geo.A_ex * N_ch_h) "200 4000")
  );
  
  // cfg for heater's tube
  parameter EntityConfig cfg_heater_tube(
    geo(
    V     = r_t_h^2 * pi * L_h,
    A_ex  = 2 * pi * r_t_h * L_h,
    L     = L_h,
    d     = 2 * r_t_h,
    N_seg = N_seg_heater,
    N_ch  = N_ch_h),
    thermo(rho_mcm = 7900 * 578.05, lambda = 20)
  );   
  
  // cooler's size of heat exchanger  
  // N_ch groups of fluid(cold, inner)-tube vs 1 gas(hot, outter) tube  
  parameter Modelica.SIunits.Radius r_i_c = 10e-3 "mm, tube's internal radius(single tube)";
  parameter Modelica.SIunits.Radius r_t_c = 15e-3 "tube's external radius(single tube)";
  parameter Modelica.SIunits.Radius r_o_c = 500e-3 "radius of external side of heat exchanger";
  parameter Modelica.SIunits.Length L_c   = 1 "m length of the tube/cooler";
  parameter Integer N_ch_c                = 100;
  parameter Integer N_seg_cooler          = N_seg;
    
  // cfg for heater's cold/fluid side
  parameter EntityConfig cfg_cooler_cold(
    geo(
    V     = r_i_c^2 * pi * L_c,
    A_ex  = 2 * pi * r_i_c * L_c,
    L     = L_c,
    d     = 2 * r_o_c,
    N_ch  = N_ch_c,
    N_seg = N_seg_cooler),
    thermo(gamma_he = 200)
  );
  
  // cfg for heater's tube
  parameter EntityConfig cfg_cooler_tube(
    geo(
    V     = r_t_c^2 * pi * L_c,
    A_ex  = 2 * pi * r_t_c * L_c,
    L     = L_c,
    d     = 2 * r_t_c,
    N_ch  = N_ch_c,
    N_seg = N_seg_cooler),
    thermo(rho_mcm = 7900 * 578.05, lambda = 20)
  );   
  
  // cfg for cooler's hot/fluid side
  parameter EntityConfig cfg_cooler_hot(
    geo(
    V     = r_o_c^2 * pi * L_c,
    A_ex  = 2 * pi * r_o_c * L_c,
    L     = L_c,
    d     = 2 * r_i_c,
    N_ch  = N_ch_c,
    N_seg = N_seg_cooler),
    thermo(gamma_he = 200 "4000")
  );  
  
  // Thermal conductivity of LTR's metal wall 
  // material inconel_750
  parameter Real table_k_LTR_wall[:, :] = [149, 16.9; 316, 20.5; 538, 26.5; 649, 28.7; 760, 31.4; 871, 35.3];
  
  parameter Real table_k_HTR_wall [:, :] = table_k_LTR_wall;
  
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
   
  function specificEnthalpy
    input ThermoState state;
    output Modelica.SIunits.SpecificEnthalpy h;
    replaceable package Medium = Modelica.Media.Interfaces.PartialMedium;
  algorithm
    h : = Medium.specificEnthalpy(Medium.setState_pTX(state.p, state.T));
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
    bc_HTR.st_hot_in   : = newThermoState_pT_CO2(p = p_pump_in, T = T_HTR_hot_in, mdot = mdot_main);
    bc_HTR.st_cold_in  : = newThermoState_pT_CO2(p = p_pump_out, T = T_HTR_cold_in, mdot = mdot_main);
    bc_HTR.st_hot_out  : = newThermoState_pT_CO2(p = p_pump_in, T = T_HTR_hot_out, mdot = mdot_main);
    bc_HTR.st_cold_out : = newThermoState_pT_CO2(p = p_pump_out, T = T_HTR_cold_out, mdot = mdot_main);
    
    
    bc_heater.st_hot_in   : = newThermoState_pT_Sodium(p = p_heater, T = T_heater_hot_in, mdot = mdot_heater);
    bc_heater.st_cold_in  : = newThermoState_pT_CO2(p = bc_HTR.st_cold_out.p, T = bc_HTR.st_cold_out.T, mdot = mdot_main);
    bc_heater.st_hot_out  : = newThermoState_pT_Sodium(p = p_heater, T = T_heater_hot_out, mdot = mdot_heater);
    bc_heater.st_cold_out : = newThermoState_pT_CO2(p = bc_HTR.st_cold_out.p, T = T_heater_cold_out, mdot = mdot_main);
    */
end PBConfiguration;
