within Steps.Model;

model SimpleCycleConfig
  "Simple cycle Configuration - follows the layout of Viv's study - comp-heater-turb-cooler"
  
  import Modelica.SIunits.Conversions.{from_degC, from_deg};
  import Modelica.SIunits.{Temperature, Pressure, SpecificEnthalpy};
  import Util = Utilities.Util;
  import Steps.Utilities.CoolProp.PropsSI; 
  import Steps.Components.PCHEGeoParam;  

  replaceable package medium_heater = Steps.Media.MoltenSalt.MoltenSalt_pT constrainedby Modelica.Media.Interfaces.PartialPureSubstance;
  replaceable package medium_main   = Steps.Media.SCO2 constrainedby Modelica.Media.Interfaces.PartialPureSubstance;
  // replaceable package medium_heater_cold = Modelica
  //replaceable package medium_cooler = Modelica.Media.Water.WaterIF97_pT;
  
  constant Real pi = Modelica.Constants.pi;
  // efficiency of main compressor, bypass_compressor and turbine
  parameter Real eta_comp = 0.89;
  parameter Real eta_turb = 0.89;
  parameter Real Ns_comp  = 2100;
  parameter Real Ns_turb  = Ns_comp;
  
  // **** Boundary Conditions as Start values for recuperators - start ****
  // Core INPUT parameters for boundary conditions **** 
  // following values are calculated by sscar for 10 MW power block( mdot_main = 125, T_amb = 35oC, split_ratio = 0.675)
  parameter Modelica.SIunits.Pressure p_source     = 8.65e6;
  parameter Modelica.SIunits.Pressure p_comp_out   = 12e6;
  parameter Modelica.SIunits.Pressure p_heater_hin = 4e6;
  // parameter Modelica.SIunits.Pressure p_sink       = 9.2e6;
  parameter Modelica.SIunits.Pressure p_sink       = p_source; // try to equal the source and sink's pressure

  parameter Modelica.SIunits.MassFlowRate mdot_main   = 10.5;
  parameter Modelica.SIunits.MassFlowRate mdot_heater = 10;
  
  parameter Modelica.SIunits.Temperature T_amb         = from_degC(35);
  parameter Modelica.SIunits.Temperature T_source      = 322;
  parameter Modelica.SIunits.Temperature T_comp_out    = 348;
  parameter Modelica.SIunits.Temperature T_heater_hin  = 595;
  parameter Modelica.SIunits.Temperature T_heater_hout = 478;
  parameter Modelica.SIunits.Temperature T_heater_cin  = 350;
  parameter Modelica.SIunits.Temperature T_heater_cout = 571;
  parameter Modelica.SIunits.Temperature T_turb_in     = 570;
  parameter Modelica.SIunits.Temperature T_sink        = 545;
  
  parameter Modelica.SIunits.Length pitch = 12.3e-3 "pitch length";
  parameter Real phi                      = 36 "pitch angle, degree";
  
  // DO NOT change following parameters - CHANGE input paramters instead  
  parameter Model.ThermoState st_source(
    p    = p_source,
    T    = T_source,
    h    = medium_main.specificEnthalpy_pT(p = p_source, T = T_source),
    // medium_main.specificEnthalpy(medium_main.setState_pT(p = st_source.p, T = st_source.T)),
    mdot = mdot_main);  

  // parameter Model.ThermoState st_comp_in = st_source;
  
  parameter Model.ThermoState st_heater_cin(
    p    = p_comp_out,
    T    = T_comp_out,
    h    = medium_main.specificEnthalpy_pT(p = st_heater_cin.p, T = st_heater_cin.T),
    mdot = mdot_main);

  parameter Model.ThermoState st_heater_cout(
    p    = p_comp_out,
    T    = T_turb_in,
    h    = medium_main.specificEnthalpy_pT(p = st_heater_cout.p, T = st_heater_cout.T),
    mdot = mdot_main);
    
  // parameter Model.ThermoState st_turb_out = st_sink;

  parameter Model.ThermoState st_sink(
    p    = p_sink,
    T    = T_sink,
    h    = medium_main.specificEnthalpy_pT(p = st_sink.p, T = st_sink.T),
    mdot = mdot_main);

  parameter Model.ThermoState st_heater_hin(
    p    = p_heater_hin,
    T    = T_heater_hin,    
    h    = medium_heater.specificEnthalpy_pT(p = st_heater_hin.p, T = st_heater_hin.T),
    mdot = mdot_heater);

  parameter Model.ThermoState st_heater_hout(
    p    = p_heater_hin,
    T    = T_heater_hout,    
    h    = medium_heater.specificEnthalpy_pT(p = st_heater_hout.p, T = st_heater_hout.T),
    mdot = mdot_heater);  

  // geometry parameters
  // In following calculation, V, A_ex are account for single tube/channel, not for total
  // check the Calculation in ThemoPower.PowerPlants.HRSG.Components.HE to understand the meaning of 
  // exsurface_G/F and extSurfaceTub, *Vol     
  parameter Real r_h                    = 1e-3;
  parameter Modelica.SIunits.Length L_h = 1 "m";
  parameter Integer N_ch                = 10000;
  parameter Integer N_seg               = 10 "default number of discretized segments in one tube";
  parameter Real l_pitch_h              = 12.3e-3;
  parameter Real a_phi_h                = 0;

  // Thermal-hydraulic properties
  // default cp and rho for alloy X-750
  parameter Modelica.SIunits.Density rho_wall             = 8280 "density of wall, kg/m3";
  parameter Modelica.SIunits.SpecificHeatCapacity cp_wall = 431 "cp of wall, J/kg-K";
  // Thermal conductivity of LTR's metal wall 
  // material inconel_750
  parameter Real table_k_heater_wall[:, :] = [149, 16.9; 316, 20.5; 538, 26.5; 649, 28.7; 760, 31.4; 871, 35.3];
  parameter Real rho_mcm                   = rho_wall * cp_wall "Metal heat capacity per unit volume [J/m^3.K]";

 
  parameter Model.AreaGeometry ga_heater_flow = SetAreaGeometry_Circle(r = r_h);
  parameter Model.PathGeometry gp_heater_flow = SetPathGeometry(geo_area = ga_heater_flow, L = L_h, N_seg = N_seg);
  parameter Model.AreaGeometry ga_heater_wall = SetAreaGeometry_Wall(r = r_h, w = r_h, h = 2 * r_h);
  parameter Model.PathGeometry gp_heater_wall = SetPathGeometry(geo_area = ga_heater_wall, L = L_h, N_seg = N_seg);

  parameter Model.HeatExchangerConfig cfg_heater(
    // identical geometry for both sides
    cfg_hot(
      st_in    = st_heater_hin,
      st_out   = st_heater_hout,
      geo_area = ga_heater_flow,
      geo_path = gp_heater_flow,
      N_ch     = N_ch,
      u        = 0,
      l_pitch  = l_pitch_h,
      a_phi    = a_phi_h
    ),     
    cfg_cold(
      st_in    = st_heater_cin,
      st_out   = st_heater_cout,
      geo_area = ga_heater_flow,
      geo_path = gp_heater_flow,
      N_ch     = N_ch,
      u        = 0,
      l_pitch  = l_pitch_h,
      a_phi    = a_phi_h
    ),
    cfg_fluid  = cfg_cold,
    cfg_gas    = cfg_hot,
    cfg_wall(
      st_init  = st_heater_hin,
      geo_area = ga_heater_wall,
      geo_wall = gp_heater_wall,
      // table_k  = table_k_heater_wall,
      rho_mcm  = rho_wall * cp_wall,
      lambda   = 200
    )
  );
  
  parameter Model.TurbomachineryConfig cfg_turb(
    st_in  = st_heater_cout,
    st_out = st_sink,
    N      = Ns_turb,
    T_nom  = 0,
    eta    = eta_turb
  );
  
  parameter Model.TurbomachineryConfig cfg_comp(
    st_in  = st_source,
    st_out = st_heater_cin,
    N      = Ns_comp,
    T_nom  = 0,
    eta    = eta_comp
  );

  // **** Boundary Conditions as Start/Nominal values for recuperators - end ****     

end SimpleCycleConfig;
