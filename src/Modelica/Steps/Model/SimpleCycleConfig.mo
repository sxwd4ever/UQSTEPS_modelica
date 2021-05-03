within Steps.Model;

model SimpleCycleConfig
  "Simple cycle Configuration - follows the layout of Viv's study - comp-heater-turb-cooler"
  
  import Modelica.SIunits.Conversions.{from_degC, from_deg};
  import Modelica.SIunits.{Temperature, Pressure, SpecificEnthalpy};
  import Util = Utilities.Util;
  import Steps.Utilities.CoolProp.PropsSI; 
  import Steps.Components.PCHEGeoParam;  

  replaceable package medium_heater = Steps.Media.MoltenSalt.MoltenSalt_pT constrainedby Modelica.Media.Interfaces.PartialPureSubstance;
  replaceable package medium_main = Steps.Media.SCO2 constrainedby Modelica.Media.Interfaces.PartialPureSubstance;
  // replaceable package medium_heater_cold = Modelica
  //replaceable package medium_cooler = Modelica.Media.Water.WaterIF97_pT;
  
  constant Real pi=Modelica.Constants.pi;
  // efficiency of main compressor, bypass_compressor and turbine
  parameter Real eta_main_compressor = 0.89;  
  parameter Real eta_turbine = 0.89;
  
  // **** Boundary Conditions as Start values for recuperators - start ****
  // Core INPUT parameters for boundary conditions **** 
  // following values are calculated by sscar for 10 MW power block( mdot_main = 125, T_amb = 35oC, split_ratio = 0.675)
  parameter Modelica.SIunits.Pressure p_source = 8.65e6;
  parameter Modelica.SIunits.Pressure p_comp_out = 12e6;
  parameter Modelica.SIunits.Pressure p_heater_hin = 4e6;
  parameter Modelica.SIunits.Pressure p_sink = 9.2e6;

  parameter Modelica.SIunits.MassFlowRate mdot_main = 10.5;
  parameter Modelica.SIunits.MassFlowRate mdot_heater = 10;  
  
  parameter Modelica.SIunits.Temperature T_amb = from_degC(35);  

  parameter Modelica.SIunits.Temperature T_source = 322;
  parameter Modelica.SIunits.Temperature T_comp_out = 348;
  parameter Modelica.SIunits.Temperature T_heater_hin = 595;
  parameter Modelica.SIunits.Temperature T_heater_hout = 478;
  parameter Modelica.SIunits.Temperature T_heater_cin = 350;
  parameter Modelica.SIunits.Temperature T_heater_cout = 571;
  parameter Modelica.SIunits.Temperature T_turb_in = 570;
  parameter Modelica.SIunits.Temperature T_sink = 545;
  
  parameter Integer N_seg = 10 "default number of discretized segments in one tube";  
  
  parameter Modelica.SIunits.Length pitch = 12.3e-3 "pitch length";
  parameter Real phi = 36 "pitch angle, degree";    
  
  // DO NOT change following parameters - CHANGE input paramters instead  
  
  // thermal state of points of cycles
  parameter Model.ThermoState st_source(
    p = p_source, 
    T = T_source, 
    h = medium_main.specificEnthalpy(medium_main.setState_pT(p = st_source.p, T = st_source.T)),
    mdot = mdot_main);

  parameter Model.ThermoState st_comp_in = st_source;
  
  parameter Model.ThermoState st_comp_out(
    p = p_comp_out, 
    T = T_comp_out, 
    h = medium_main.specificEnthalpy(medium_main.setState_pT(p = st_comp_out.p, T = st_comp_out.T)),
    mdot = mdot_main);       

  parameter Model.ThermoState st_turb_in(
    p = p_comp_out, 
    T = T_turb_in, 
    h = medium_main.specificEnthalpy(medium_main.setState_pT(p = st_turb_in.p, T = st_turb_in.T)),
    mdot = mdot_main);    
    
  parameter Model.ThermoState st_turb_out = st_sink;   

  parameter Model.ThermoState st_sink(
    p = p_sink, 
    T = T_sink, 
    h = medium_main.specificEnthalpy(medium_main.setState_pT(p = st_sink.p, T = st_sink.T)),
    mdot = mdot_main);     

  parameter Model.ThermoState st_heater_hin(
    p = p_heater_hin, 
    T = T_heater_hin, 
    h = medium_heater.specificEnthalpy(medium_heater.setState_pT(p = st_heater_hin.p, T = st_heater_hin.T)),
    mdot = mdot_heater);   

  parameter Model.ThermoState st_heater_hout(
    p = p_heater_hin, 
    T = T_heater_hout, 
    h = medium_heater.specificEnthalpy(medium_heater.setState_pT(p = st_heater_hout.p, T = st_heater_hout.T)),
    mdot = mdot_heater);         


  // boundary condition for heater
  parameter HEBoundaryCondition bc_heater(
    st_hot_in = st_heater_hin,
    st_hot_out = st_heater_hout,
    st_cold_in = st_comp_out,
    st_cold_out = st_turb_in);  
  
  // **** Boundary Conditions as Start/Nominal values for recuperators - end ****     
     
  // heater's size of heat exchanger
  // N_ch groups of fluid(hot, inner)-tube vs 1 gas(cold, outter) tube  
  
  // In following calculation, V, A_ex are account for single tube/channel, not for total
  // check the Calculation in ThemoPower.PowerPlants.HRSG.Components.HE to understand the meaning of 
  // exsurface_G/F and extSurfaceTub, *Vol     
  parameter Real r_h = 1e-3;  
  parameter Real p_h = (pi + 2) * r_h;
  parameter Real V_h = r_h^2 * pi * L_h / 2;
  parameter Real A_h = p_h * L_h;   

  parameter Modelica.SIunits.Length L_h = 1 "m"; 
  parameter Integer N_ch_h = 10000;  
  parameter Integer N_seg_h = N_seg;  

  // cfg for heater's hot/fluid side
  parameter EntityConfig cfg_heater_hot(
    geo(
      V = V_h, // * N_ch_h, 
      A_ex = A_h, // * N_ch_h, exchange surface between fluid-tube
      L = L_h, 
      d = r_h * 2,
      N_seg = N_seg_h, 
      N_ch = N_ch_h), // calculate by steps the python code
    thermo(gamma_he = 1.938761018e6/A_h "200")
  );
    
  // cfg for heater's cold/fluid side
  parameter EntityConfig cfg_heater_cold(
    geo(
      V = V_h, // * N_ch_h, 
      A_ex = A_h, // * N_ch_h, exchange surface between fluid-tube
      L = L_h, 
      d = r_h * 2,
      N_seg = N_seg_h, 
      N_ch = N_ch_h), // calculate by steps the python code
    thermo(gamma_he = 1.938761018e6/A_h "200")
  );
  
  // cfg for heater's tube
  parameter EntityConfig cfg_heater_tube(
    geo(
      V = V_h, // r_o_h^2 * pi * L_h * N_ch_h - cfg_h_tube.geo.V , 
      A_ex = A_h, // * N_ch_h, 
      L = L_h, 
      d = r_h * 2,
      N_seg = N_seg_h, 
      N_ch = N_ch_h),
    thermo(rho_mcm = rho_mcm, gamma_he = 1.938761018e6/A_h "200")
  );     
  
  // Thermal conductivity of LTR's metal wall 
  // material inconel_750
  parameter Real table_k_heater_wall[:, :] = [149, 16.9; 316, 20.5; 538, 26.5; 649, 28.7; 760, 31.4; 871, 35.3];  
  
    // default cp and rho for alloy X-750
  parameter Modelica.SIunits.Density rho_wall = 8280 "density of wall, kg/m3";
  parameter Modelica.SIunits.SpecificHeatCapacity cp_wall = 431 "cp of wall, J/kg-K";
  
  parameter Real rho_mcm = rho_wall * cp_wall "Metal heat capacity per unit volume [J/m^3.K]";
    
end SimpleCycleConfig;
