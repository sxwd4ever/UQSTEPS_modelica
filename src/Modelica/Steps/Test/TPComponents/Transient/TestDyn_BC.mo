within Steps.Test.TPComponents.Transient;

model TestDyn_BC
  "Test boundary condition set-up components such as source and sink"  
  import Modelica.SIunits.Conversions.{from_degC, from_deg};
  import Modelica.SIunits.{Temperature, Pressure, SpecificEnthalpy};
  import Util = Utilities.Util;
  import Steps.Utilities.CoolProp.PropsSI;  
  import Steps.Components.{PCHEGeoParam};
  import Steps.Model.{PBConfiguration, SimParam, EntityConfig, EntityGeoParam, EntityThermoParam, ThermoState, HEBoundaryCondition} ;
  import Model.PBConfiguration;
  import ThermoPower.Choices.Init.Options;
  import ThermoPower.System;
  import ThermoPower.Gas;
  
  
  // package medium_main = Modelica.Media.IdealGases.SingleGases.CO2; //Steps.Media.CO2;
  package medium_main = Steps.Media.SCO2(
    // inputChoice = ExternalMedia.Common.InputChoice.pT,
    // substanceNames = {"CO2|debug=40"}
  );  
  // package medium_main = ExternalMedia.Examples.CO2CoolProp;
  // package medium_heater = Steps.Media.ThermiaOilD; // out of working range of this 10Mw high T loop
  package medium_heater = SolarTherm.Media.Sodium.Sodium_pT;
  package medium_cooler = ThermoPower.Water.StandardWater;//Modelica.Media.IdealGases.SingleGases.CO2;  
 
  // package medium_hot = ExternalMedia.Examples.CO2CoolProp;
  // package medium_cold = ExternalMedia.Examples.CO2CoolProp;    
  
  // package medium_hot = Modelica.Media.IdealGases.SingleGases.CO2;
  // package medium_cold = Modelica.Media.IdealGases.SingleGases.CO2;
  
  // geometry parameters
  constant  Real pi           = Modelica.Constants.pi;
  // parameter Integer N_ch      = integer(2400*1e1) "channel number";
  parameter Integer N_ch      = integer(2e6) "channel number";
  parameter Integer N_seg     = 10 "number of segments in one tube";
  parameter SI.Length D_ch    = 1.72e-3  "channel diameter, semi circular tube";
  parameter SI.Length r_ch    = D_ch / 2 "channel radiaus";
  parameter SI.Length L_fp    = 270 * 1e-3 * 2 "equivalent valid channel flow path length";
  parameter SI.Length L_pitch = 12e-3 "pitch length";
  parameter Real a_phi        = 0 "pitch angle °";
  parameter SI.Length H_ch    = 4.17e-3 "Height of the solid domain, containing one cold tube and one hot tube";
  parameter SI.Length W_ch    = 2.3e-3 "Width of the solid domain";
  parameter SI.Length T_wall  = 0.51e-3 "Wall thinckness";
  parameter SI.Length L_wall  = 420e-3 * 2 "Length of wall, not necessarily equals to length of flow path";
  parameter SI.Area A         = pi * r_ch ^2 / 2 "Area of cross section of semi circular tube";
  
  // Stainless 316, 316L, 317, 317L
  parameter Modelica.SIunits.Density rho_wall             = 8030 "density of wall, kg/m3";
  parameter Modelica.SIunits.SpecificHeatCapacity cp_wall = 485 "cp of wall, J/kg-K";
  // thermal conductivity (T in K) https://www.theworldmaterial.com/aisi-316-ss316-stainless-steel-properties-composition/
  // parameter Real table_k_metalwall[:,:] = [20, 12.1; 100, 16.3; 500, 21.5];
  parameter Real table_k_metalwall[:,:] = [293.15, 12.1; 373.15, 16.3; 773.15, 21.5];
  parameter Real Cf_C1                  = 1.626, Cf_C2 = 1, Cf_C3 = 1;
  // parameter Real Cf_C1_cold = 1, Cf_C2_cold = 1, Cf_C3_cold = 1;
  parameter Real use_rho_bar  = -1.0;
  parameter Real rho_bar_hot  = 1.0;
  parameter Real rho_bar_cold = 1.0;

  // input parameters of the power block
  parameter Modelica.SIunits.MassFlowRate mdot_main = 125 "kg/s, mass flow in the main path of PB, which follows the power demand";
  parameter Modelica.SIunits.MassFlowRate mdot_heater_hot = 90 "kg/s, mass flow rate of heater's hot fluid";
  parameter Real gamma = 0.4 "split ratio, mdot_bypass/mdot_main";
  
  parameter Modelica.SIunits.Temperature T_heater_hot = from_degC(800) "K, Temperature of heater's hot fluid";  
  parameter Modelica.SIunits.Temperature T_cooler_cold = from_degC(45) "K, Temperature of cooler's cold fluid";
  parameter Integer T_step = 3;
  
  // select the configuration of parameters
  parameter Model.RCBCycleConfig cfg(
    redeclare package medium_heater = medium_heater,
    redeclare package medium_main   = medium_main,
    // mdot_heater      = 40,
    // T_heater_hot_in  = from_degC(800),
    // T_heater_hot_out = from_degC(600),
    r_i_heater  = 1e-3,
    r_t_heater  = 2e-3, //cfg.r_i_heater + 10e-3,
    r_o_heater  = 3e-3,                      // agree with the final parameter Dhyd = 1 in HE, should be checked to see if it is capable of containing all fluid-metal tubes
    N_ch_heater = 10000,
    L_heater    = 1,
    N_ch_HTR    = 30000,
    L_HTR       = 2.5,
    r_i_HTR     = 1.5e-3,
    r_o_HTR     = 1.5e-3,
    N_ch_LTR    = 30000,
    L_LTR       = 2.5,
    r_i_LTR     = 1.5e-3,
    r_o_LTR     = 1.5e-3,
    Ns_turb     = 30000,
    // latest boundary conditions, following values are simulation results with sourcePressure.p = 200 bar and above geometry params
    p_comp_in  = 109.59e5,
    p_comp_out = 20e6,    
    p_heater   = 20e6,    
    T_HTR_hot_in      = from_degC(556.322),
    T_HTR_cold_out    = from_degC(515.234),
    T_HTR_hot_out     = from_degC(296.146),
    T_HTR_cold_in     = from_degC(266.181),
    T_LTR_cold_in     = from_degC(118.13),
    T_LTR_hot_out     = from_degC(154.909),
    T_heater_hot_in   = from_degC(800),
    T_heater_hot_out  = from_degC(707.918),
    T_heater_cold_out = from_degC(637.551),
    T_cooler_cold_out = from_degC(61.7479462700001),
    T_cooler_hot_out  = from_degC(68.494),
    
    mdot_main   = 128.774,
    mdot_comp   = 88.0661,
    mdot_heater = 40,
    mdot_cooler = 40.7188
  );
  
  // set the values of parameters accordingly
  parameter Model.HeatExchangerConfig cfg_heater = cfg.cfg_heater;
  
  parameter Model.ThermoState st_source  = cfg_heater.cfg_cold.st_in;
  parameter Model.ThermoState st_sink   = cfg_heater.cfg_cold.st_in;  

  parameter Integer N_seg_heater = cfg.cfg_heater.cfg_hot.geo_path.N_seg; 
  
  //Components
  // global init opition (system.initOpt) leads to order reduction error
  // use this flag to control the initialization of all components instead. 
  parameter Boolean SSInit = true "Steady-state initialization";  
  inner ThermoPower.System system(
    allowFlowReversal = false,
    initOpt           = ThermoPower.Choices.Init.Options.steadyState)
    annotation(
    Placement(transformation(extent = {{80, 80}, {100, 100}})));    
  
  ThermoPower.Gas.SourceMassFlow source(
    redeclare package Medium = medium_main, 
    T         = st_source.T,
    p0        = st_source.p,
    w0        = st_source.mdot,
    use_in_T  = false,
    use_in_w0 = true,
    gas(
      p(start = st_source.p, nominal = st_source.p), 
      T(start = st_source.T, nominal = st_source.T)))
  annotation(
    Placement(transformation(extent = {{-70, -10}, {-50, 10}}, rotation = 0))); 
    
  ThermoPower.Gas.SinkPressure sink(
    redeclare package Medium = medium_main, 
    p0 = st_sink.p,
    T  = st_sink.T)
  annotation(
    Placement(transformation(origin = {0, -80}, extent = {{-10, -10}, {10, 10}}, rotation = 270)));
  
  // value to accelerate the simulation
  constant Integer SEC_PER_HOUR = integer(60 * 60 / time_scaling_factor); 
  constant Real time_scaling_factor = 12; // 5 min = 1 hour
  // constant Real time_scaling_factor = 1; // 1 hour = 1 hour  

  // ramp input to simulate the ramp change in ASPEN+ simulation
  // ramp change for case I in Guan's report
  Steps.TPComponents.RampSeq_W ramp_mdot(    
    time_start = 3 * SEC_PER_HOUR,
    interval   = 1 * SEC_PER_HOUR,
    duration   = 0.15 * SEC_PER_HOUR,
    offset     = st_source.mdot
  );
  /*
  // variable for validation
  Modelica.SIunits.Power Q_out = (hE.gasIn.h_outflow - hE.gasOut.h_outflow) * hE.gasIn.m_flow; 
  Modelica.SIunits.Power Q_in = (hE.waterOut.h_outflow - hE.waterIn.h_outflow) * hE.waterIn.m_flow;
  Boolean isQMatch = abs(Q_out -Q_in) < 1e-3;
  */
initial equation
//hstart_F_Out = hE.waterOut.h_outflow;
equation
  connect(source.flange, sink.flange);
  connect(ramp_mdot.y, source.in_w0);

annotation(
    Diagram(graphics),
    experiment(StartTime = 0, StopTime = 2100, Tolerance = 1e-2, Interval = 30),
    __OpenModelica_commandLineOptions = "--matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection -d=initialization,NLSanalyticJacobian,aliasConflicts",    
    __OpenModelica_simulationFlags(lv = "LOG_DEBUG,LOG_NLS,LOG_NLS_V,LOG_STATS,LOG_INIT,LOG_STDOUT, -w", outputFormat = "mat", s = "dassl", nls = "homotopy"));
end TestDyn_BC;
