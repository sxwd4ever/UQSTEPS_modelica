within Steps.Test;

model TestTP_HE
  "Test for HE in ThermoPower"  
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
  
  package medium_hot = Steps.Media.CO2;
  //package medium_cold = Steps.Media.CO2;
  // package medium_hot = Steps.Media.CO2;
  // package medium_cold = Steps.Media.CO2;
  // package medium_heater = SolarTherm.Media.Sodium.Sodium_pT;
  // package medium_heater = Steps.Media.CO2;
  package medium_heater = SolarTherm.Media.Sodium.Sodium_pT; // ThermoPower.Water.StandardWater;//Modelica.Media.IdealGases.SingleGases.CO2;
  package medium_cold = Steps.Media.CO2; // Modelica.Media.IdealGases.SingleGases.CO2;
  package medium_cooler = ThermoPower.Water.StandardWater; 
  
  // parameter for C++ implementation of PCHE - based on Modelica impl's result
  parameter Model.PBConfiguration cfg_on_design( 
  p_pump_in = 8e6,
  p_pump_out = 20e6,  
  mdot_main = 51.51,
  mdot_pump = 31.31,
  mdot_heater = 20,
  T_HTR_cold_in = from_degC(141.3), 
  T_HTR_cold_out = from_degC(495.302),
  T_HTR_hot_out = from_degC(141.041),
  T_LTR_hot_out = from_degC(63.6726));    

  // for high temperature fluid - low temperature gas test
  parameter Model.PBConfiguration cfg_HFLG_test( 
  p_pump_in = 8e6,
  p_pump_out = 20e6,  
  mdot_main = 100,
  mdot_pump = 100,
  mdot_heater = 100,
  T_HTR_cold_in = from_degC(141.3), 
  T_HTR_cold_out = from_degC(495.302),
  T_HTR_hot_out = from_degC(141.041),
  T_LTR_hot_out = from_degC(63.6726),
  T_heater_hot_in = from_degC(600),
  T_heater_hot_out = from_degC(500),
  T_heater_cold_out = from_degC(550)  
  ); 
  
  parameter Model.PBConfiguration cfg_liquid_sodium_test( 
  p_pump_in = 20e6,
  p_pump_out = 20e6,  
  p_heater = 20e6, 
  mdot_main = 100,
  mdot_pump = 100,
  mdot_heater = 100,
  T_HTR_cold_in = from_degC(141.3), 
  T_HTR_cold_out = from_degC(400),
  T_HTR_hot_out = from_degC(141.041),
  T_LTR_hot_out = from_degC(63.6726),
  T_heater_hot_in = from_degC(800),
  T_heater_hot_out = from_degC(700),
  T_heater_cold_out = from_degC(500),
  cfg_heater_hot(thermo(gamma_he = 200))
  );  
  
  parameter Model.PBConfiguration cfg_hxb_test( 
  p_pump_in = 6.6e5,
  p_pump_out = 20e6,  
  mdot_main = 585.5,
  mdot_pump = 31.31,
  mdot_heater = 6.91,
  T_HTR_cold_in = from_degC(141.3), 
  T_HTR_cold_out = 519.15,
  T_HTR_hot_out = from_degC(141.041),
  T_LTR_hot_out = from_degC(63.6726),
  T_heater_hot_in = 435.75 - 100.0,
  T_heater_hot_out = 505.04,
  T_heater_cold_out = 517.44);    
    
  // select the configuration of parameters
  parameter Model.PBConfiguration cfg = cfg_on_design;
  
  // set the values of parameters accordingly
  parameter HEBoundaryCondition bc_heater = cfg.bc_heater;
  
  // use HTR's geo parameters as default 
  parameter EntityGeoParam geo_heater_hot = cfg.cfg_heater_hot.geo;
  parameter EntityGeoParam geo_heater_cold = cfg.cfg_heater_cold.geo;
  parameter EntityGeoParam geo_heater_tube = cfg.cfg_heater_tube.geo;
  
  parameter EntityThermoParam thermo_hot = cfg.cfg_heater_hot.thermo;
  parameter EntityThermoParam thermo_cold = cfg.cfg_heater_cold.thermo;
  parameter EntityThermoParam thermo_heater = cfg.cfg_heater_tube.thermo;
  
  ThermoPower.Water.SourceMassFlow source_heater_hot(
  redeclare package Medium = medium_heater,
    w0 = bc_heater.st_hot_in.mdot,
    p0 = bc_heater.st_hot_in.p,
    T = bc_heater.st_hot_in.T, 
    use_T = true) 
    annotation(
    Placement(transformation(origin = {0, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 270)));
    
  ThermoPower.Water.SinkPressure sink_heater_hot(
  redeclare package Medium = medium_heater,
    p0 = bc_heater.st_hot_out.p,
    T = bc_heater.st_hot_out.T,
    use_T = true) 
    annotation(
    Placement(transformation(origin = {0, -80}, extent = {{-10, -10}, {10, 10}}, rotation = 270)));
    
  ThermoPower.Gas.SourceMassFlow source_cold(
  redeclare package Medium = medium_cold,
    T = bc_heater.st_cold_in.T, 
    p0 = bc_heater.st_cold_in.p, 
    use_in_T = false, 
    w0 = bc_heater.st_cold_in.mdot) 
  annotation(
    Placement(transformation(origin = {0, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 270)));      
  
  ThermoPower.Gas.SinkPressure sink_cold(
  redeclare package Medium = medium_cold,
    p0 = bc_heater.st_cold_out.p, 
    T = bc_heater.st_cold_out.T) 
  annotation(
    Placement(transformation(origin = {0, -80}, extent = {{-10, -10}, {10, 10}}, rotation = 270)));  
  
  ThermoPower.PowerPlants.HRSG.Components.HE hE(
  //Components.HEG2G hE(
    redeclare package FluidMedium = medium_heater, 
    redeclare package FlueGasMedium = medium_cold, 
    fluidFlow(fixedMassFlowSimplified = true, hstartin = bc_heater.st_hot_in.h, hstartout=bc_heater.st_hot_out.h), // set the fluid flow as fixed mdot for simplarity
    N_G=geo_heater_cold.N_seg,
    N_F=geo_heater_hot.N_seg,
    Nw_G=geo_heater_tube.N_seg,
    gasNomFlowRate=bc_heater.st_cold_in.mdot,
    gasNomPressure=bc_heater.st_cold_in.p,
    fluidNomFlowRate=bc_heater.st_hot_in.mdot,
    fluidNomPressure=bc_heater.st_hot_in.p,
    exchSurface_G=geo_heater_cold.A_ex,
    exchSurface_F=geo_heater_hot.A_ex,
    extSurfaceTub=geo_heater_tube.A_ex,
    gasVol=geo_heater_cold.V,
    fluidVol=geo_heater_hot.V,
    metalVol=geo_heater_tube.V,
    SSInit=false,
    rhomcm=thermo_heater.rho_mcm,
    lambda=thermo_heater.lambda,
    Tstartbar_G=bc_heater.st_cold_in.T,
    Tstartbar_M=bc_heater.st_hot_in.T - 50,
    pstart_F = bc_heater.st_hot_in.p, 
    FluidPhaseStart=ThermoPower.Choices.FluidPhase.FluidPhases.Liquid,    
    redeclare replaceable model HeatTransfer_F = ThermoPower.Thermal.HeatTransferFV.IdealHeatTransfer, //ConstantHeatTransferCoefficientTwoGrids(gamma = thermo_hot.gamma_he),     
    redeclare replaceable model HeatTransfer_G = ThermoPower.Thermal.HeatTransferFV.IdealHeatTransfer, //ConstantHeatTransferCoefficient(gamma =  thermo_cold.gamma_he),     
    redeclare model HeatExchangerTopology = ThermoPower.Thermal.HeatExchangerTopologies.CounterCurrentFlow,
    metalTube(WallRes=false)) annotation(
    Placement(transformation(extent = {{-20, -20}, {20, 20}}, rotation = 0)));

  inner ThermoPower.System system(allowFlowReversal = false, initOpt=ThermoPower.Choices.Init.Options.noInit) annotation(
    Placement(transformation(extent = {{80, 80}, {100, 100}})));
  
  //ThermoPower.Water.SensT T_waterOut(redeclare package Medium = medium_heater) annotation(
  //  Placement(transformation(origin = {4, -50}, extent = {{-10, -10}, {10, 10}}, rotation = 270)));
  ThermoPower.Gas.SensT T_gasOut(redeclare package Medium = medium_cold) annotation(
    Placement(transformation(extent = {{30, -6}, {50, 14}}, rotation = 0)));
initial equation
//hstart_F_Out = hE.waterOut.h_outflow;
equation
  connect(source_heater_hot.flange, hE.waterIn) annotation(
    Line(points = {{-1.83697e-015, 50}, {-1.83697e-015, 20}, {0, 20}}, color = {0, 0, 255}, thickness = 0.5, smooth = Smooth.None));
  
  connect(hE.waterOut, sink_heater_hot.flange) annotation(
    Line(points = {{1.83697e-015, -70}, {1.83697e-015, -56}, {-8.88178e-016, -56}}, thickness = 0.5, color = {0, 0, 255}));
   
  // connect(hE.waterOut, T_waterOut.inlet) annotation(
    //Line(points = {{8.88178e-016, -44}, {8.88178e-016, -20}, {0, -20}}, thickness = 0.5, color = {0, 0, 255}));

  // connect(T_waterOut.outlet, sink_heater_hot.flange) annotation(
    // Line(points = {{1.83697e-015, -70}, {1.83697e-015, -56}, {-8.88178e-016, -56}}, thickness = 0.5, color = {0, 0, 255}));
  
  connect(source_cold.flange, hE.gasIn) annotation(
    Line(points = {{-50, 0}, {-20, 0}}, color = {159, 159, 223}, thickness = 0.5, smooth = Smooth.None));
  connect(hE.gasOut, T_gasOut.inlet) annotation(
    Line(points = {{34, 0}, {34, 0}, {20, 0}}, color = {159, 159, 223}, thickness = 0.5));
  connect(T_gasOut.outlet, sink_cold.flange) annotation(
    Line(points = {{46, 0}, {46, 0}, {60, 0}}, color = {159, 159, 223}, thickness = 0.5));

annotation(
    Diagram(graphics),
    experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-2, Interval = 1),
    __OpenModelica_commandLineOptions = "--matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection -d=initialization,NLSanalyticJacobian,aliasConflicts,bltdump",    
    __OpenModelica_simulationFlags(lv = "LOG_DEBUG,LOG_NLS,LOG_NLS_V,LOG_STATS,LOG_INIT,LOG_STDOUT, -w", outputFormat = "mat", s = "dassl", nls = "homotopy"));
end TestTP_HE;
