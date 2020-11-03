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
  package medium_cold = Steps.Media.CO2;
  // package medium_hot = Steps.Media.CO2;
  // package medium_cold = Steps.Media.CO2;
  // package medium_heater = SolarTherm.Media.Sodium.Sodium_pT;
  package medium_heater = Steps.Media.CO2;
  // package medium_heater = ThermoPower.Water.StandardWater;// Modelica.Media.IdealGases.SingleGases.CO2;

  
  // parameter for C++ implementation of PCHE - based on Modelica impl's result
  parameter Model.PBConfiguration cfg_on_design( 
  redeclare package medium_heater = medium_hot,
  p_pump_in = 8e6,
  p_pump_out = 20e6,  
  mdot_main = 51.51,
  mdot_pump = 31.31,
  T_HTR_cold_in = from_degC(141.3), 
  T_HTR_cold_out = from_degC(495.302),
  T_HTR_hot_out = from_degC(141.041),
  T_LTR_hot_out = from_degC(63.6726));    
  
  // select the configuration of parameters
  parameter Model.PBConfiguration cfg = cfg_on_design;
  
  // set the values of parameters accordingly
  parameter HEBoundaryCondition bc_HTR = cfg.bc_HTR;  
  parameter HEBoundaryCondition bc_LTR = cfg.bc_LTR;
  parameter HEBoundaryCondition bc_heater = cfg.bc_heater;
  
  parameter ThermoState st_bypass = cfg.st_bypass;
  
  parameter EntityGeoParam geo_hot = cfg.cfg_HTR_hot.geo;
  parameter EntityGeoParam geo_cold = cfg.cfg_HTR_cold.geo;
  parameter EntityGeoParam geo_tube = cfg.cfg_HTR_tube.geo;
  
  // use HTR's geo parameters as default 
  parameter EntityGeoParam geo_heater_hot = cfg.cfg_HTR_hot.geo;
  parameter EntityGeoParam geo_heater_cold = cfg.cfg_HTR_cold.geo;
  parameter EntityGeoParam geo_heater_tube = cfg.cfg_HTR_tube.geo;
  
  parameter EntityThermoParam thermo_hot = cfg.cfg_HTR_hot.thermo;
  parameter EntityThermoParam thermo_cold = cfg.cfg_HTR_cold.thermo;
  parameter EntityThermoParam thermo_tube = cfg.cfg_HTR_tube.thermo;  
  parameter EntityThermoParam thermo_mixer = cfg.cfg_mixer.thermo;
  parameter EntityThermoParam thermo_heater = cfg.cfg_heater.thermo;

  package FlueGasMedium = Modelica.Media.IdealGases.SingleGases.CO2;
  package FluidMedium = ThermoPower.Water.StandardWater;
  
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
    T = bc_HTR.st_cold_out.T, 
    p0 = bc_HTR.st_cold_out.p, 
    use_in_T = false, 
    w0 = bc_HTR.st_cold_out.mdot) 
  annotation(
    Placement(transformation(origin = {0, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 270)));      
  
  ThermoPower.Gas.SinkPressure sink_cold(
    redeclare package Medium = medium_cold, 
    p0 = bc_heater.st_cold_out.p, 
    T = bc_heater.st_cold_out.T) 
  annotation(
    Placement(transformation(origin = {0, -80}, extent = {{-10, -10}, {10, 10}}, rotation = 270)));  
  
  ThermoPower.PowerPlants.HRSG.Components.HE hE(
    redeclare package FluidMedium = medium_cold, 
    redeclare package FlueGasMedium = medium_heater, 
    fluidFlow(fixedMassFlowSimplified = true), // set the fluid flow as fixed mdot for simplarity
    N_G=geo_heater_hot.N_seg,
    N_F=geo_heater_cold.N_seg,
    Nw_G = geo_heater_tube.N_seg,
    gasNomFlowRate=bc_heater.st_hot_in.mdot,
    gasNomPressure=bc_heater.st_hot_in.p,
    fluidNomFlowRate=bc_heater.st_cold_in.mdot,
    fluidNomPressure=bc_heater.st_cold_in.p,
    exchSurface_G=geo_heater_hot.A_ex,
    exchSurface_F=geo_heater_cold.A_ex,
    extSurfaceTub=geo_heater_tube.A_ex,
    gasVol=geo_heater_hot.V,
    fluidVol=geo_heater_cold.V,
    metalVol=geo_heater_tube.V,
    SSInit=false,
    rhomcm=thermo_heater.rho_mcm,
    lambda=thermo_heater.lambda,
    Tstartbar_G=bc_heater.st_hot_in.T,
    Tstartbar_M=bc_heater.st_hot_in.T - 50,
    pstart_F = bc_heater.st_cold_in.p, 
    FluidPhaseStart=ThermoPower.Choices.FluidPhase.FluidPhases.Steam,    
    redeclare replaceable model HeatTransfer_F = ThermoPower.Thermal.HeatTransferFV.ConstantHeatTransferCoefficient(gamma = thermo_cold.gamma_he),     
    redeclare replaceable model HeatTransfer_G = ThermoPower.Thermal.HeatTransferFV.ConstantHeatTransferCoefficientTwoGrids(gamma = thermo_hot.gamma_he),     
    redeclare model HeatExchangerTopology = ThermoPower.Thermal.HeatExchangerTopologies.CounterCurrentFlow) annotation(
    Placement(transformation(extent = {{-20, -20}, {20, 20}}, rotation = 0)));

  inner ThermoPower.System system(allowFlowReversal = false) annotation(
    Placement(transformation(extent = {{80, 80}, {100, 100}})));
   /*
  ThermoPower.Water.SourceMassFlow source_heater_hot(redeclare package Medium = FluidMedium, w0 = fluidNomFlowRate, p0 = fluidNomPressure, h = hstart_F_In) annotation(
    Placement(transformation(origin = {0, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 270)));
  ThermoPower.Water.SinkPressure sink_heater_hot(redeclare package Medium = FluidMedium, p0 = fluidNomPressure, h = hstart_F_Out) annotation(
    Placement(transformation(origin = {0, -80}, extent = {{-10, -10}, {10, 10}}, rotation = 270)));

  ThermoPower.Gas.SinkPressure sink_cold(redeclare package Medium = FlueGasMedium, T = Tstart_G_Out, p0 = 10e6) annotation(
    Placement(transformation(extent = {{60, -10}, {80, 10}}, rotation = 0)));
  ThermoPower.Gas.SourceMassFlow source_cold(redeclare package Medium = FlueGasMedium, T = Tstart_G_In, p0 = 10e6, w0 = gasNomFlowRate) annotation(
    Placement(transformation(extent = {{-70, -10}, {-50, 10}}, rotation = 0)));
    */
  ThermoPower.Water.SensT T_waterOut(redeclare package Medium = FluidMedium) annotation(
    Placement(transformation(origin = {4, -50}, extent = {{-10, -10}, {10, 10}}, rotation = 270)));
  ThermoPower.Gas.SensT T_gasOut(redeclare package Medium = FlueGasMedium) annotation(
    Placement(transformation(extent = {{30, -6}, {50, 14}}, rotation = 0)));
initial equation
//hstart_F_Out = hE.waterOut.h_outflow;
equation
  connect(T_gasOut.inlet, hE.gasOut) annotation(
    Line(points = {{34, 0}, {34, 0}, {20, 0}}, color = {159, 159, 223}, thickness = 0.5));
  connect(T_gasOut.outlet, sink_cold.flange) annotation(
    Line(points = {{46, 0}, {46, 0}, {60, 0}}, color = {159, 159, 223}, thickness = 0.5));
  connect(sink_heater_hot.flange, T_waterOut.outlet) annotation(
    Line(points = {{1.83697e-015, -70}, {1.83697e-015, -56}, {-8.88178e-016, -56}}, thickness = 0.5, color = {0, 0, 255}));
  connect(T_waterOut.inlet, hE.waterOut) annotation(
    Line(points = {{8.88178e-016, -44}, {8.88178e-016, -20}, {0, -20}}, thickness = 0.5, color = {0, 0, 255}));
  connect(source_heater_hot.flange, hE.waterIn) annotation(
    Line(points = {{-1.83697e-015, 50}, {-1.83697e-015, 20}, {0, 20}}, color = {0, 0, 255}, thickness = 0.5, smooth = Smooth.None));
  connect(source_cold.flange, hE.gasIn) annotation(
    Line(points = {{-50, 0}, {-20, 0}}, color = {159, 159, 223}, thickness = 0.5, smooth = Smooth.None));

annotation(
    Diagram(graphics),
    experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-3, Interval = 1),
    __OpenModelica_commandLineOptions = "--matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection -d=initialization,NLSanalyticJacobian,aliasConflicts,bltdump",    
    __OpenModelica_simulationFlags(lv = "LOG_DEBUG,LOG_NLS,LOG_NLS_V,LOG_STATS,LOG_INIT,LOG_STDOUT, -w", outputFormat = "mat", s = "dassl", nls = "homotopy"));
end TestTP_HE;
