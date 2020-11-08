within Steps.Test;

model TestTP_Cooler
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
  package medium_cooler = ThermoPower.Water.StandardWater;//Modelica.Media.IdealGases.SingleGases.CO2;  
  
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
    
  // select the configuration of parameters
  parameter Model.PBConfiguration cfg = cfg_on_design;
  
  // set the values of parameters accordingly
  parameter HEBoundaryCondition bc = cfg.bc_cooler;
  
  // use HTR's geo parameters as default 
  parameter EntityGeoParam geo_hot = cfg.cfg_heater_hot.geo;
  parameter EntityGeoParam geo_cold = cfg.cfg_heater_cold.geo;
  parameter EntityGeoParam geo_tube = cfg.cfg_heater_tube.geo;
  
  parameter EntityThermoParam thermo_hot = cfg.cfg_heater_hot.thermo;
  parameter EntityThermoParam thermo_cold = cfg.cfg_heater_cold.thermo;
  parameter EntityThermoParam thermo_heater = cfg.cfg_heater_tube.thermo;
  
  ThermoPower.Gas.SourceMassFlow source_hot(
  redeclare package Medium = medium_hot,
    w0 = bc.st_hot_in.mdot,
    p0 = bc.st_hot_in.p,
    T = bc.st_hot_in.T) 
    annotation(
    Placement(transformation(origin = {0, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 270)));
    
  ThermoPower.Gas.SinkPressure sink_hot(
  redeclare package Medium = medium_hot,
    p0 = bc.st_hot_out.p,
    T = bc.st_hot_out.T) 
    annotation(
    Placement(transformation(origin = {0, -80}, extent = {{-10, -10}, {10, 10}}, rotation = 270)));
    
  ThermoPower.Water.SourceMassFlow source_cold(
  redeclare package Medium = medium_cooler,
    T = bc.st_cold_in.T, 
    p0 = bc.st_cold_in.p, 
    use_T = true,
    use_in_T = false, 
    w0 = bc.st_cold_in.mdot) 
  annotation(
    Placement(transformation(origin = {0, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 270)));      
  
  ThermoPower.Water.SinkPressure sink_cold(
  redeclare package Medium = medium_cooler,
    p0 = bc.st_cold_out.p, 
    T = bc.st_cold_out.T,
    use_T = true) 
  annotation(
    Placement(transformation(origin = {0, -80}, extent = {{-10, -10}, {10, 10}}, rotation = 270)));  
  
  ThermoPower.PowerPlants.HRSG.Components.HE hE(
  //Components.HEG2G hE(
    redeclare package FluidMedium = medium_cooler, 
    redeclare package FlueGasMedium = medium_hot, 
    fluidFlow(fixedMassFlowSimplified = true), // hstartin = bc.st_cold_in.h, hstartout=bc.st_cold_out.h), // set the fluid flow as fixed mdot for simplarity
    N_G=geo_hot.N_seg,
    N_F=geo_cold.N_seg,
    Nw_G=geo_tube.N_seg,
    gasNomFlowRate=bc.st_hot_in.mdot,
    gasNomPressure=bc.st_hot_in.p,
    fluidNomFlowRate=bc.st_cold_in.mdot,
    fluidNomPressure=bc.st_cold_in.p,
    exchSurface_G=geo_hot.A_ex,
    exchSurface_F=geo_cold.A_ex,
    extSurfaceTub=geo_tube.A_ex,
    gasVol=geo_hot.V,
    fluidVol=geo_cold.V,
    metalVol=geo_tube.V,
    SSInit=false,
    rhomcm=thermo_heater.rho_mcm,
    lambda=thermo_heater.lambda,
    Tstartbar_G=bc.st_hot_in.T,
    Tstartbar_M=bc.st_hot_in.T - 50,
    pstart_F = bc.st_cold_in.p, 
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
  ThermoPower.Gas.SensT T_gasOut(redeclare package Medium = medium_hot) annotation(
    Placement(transformation(extent = {{30, -6}, {50, 14}}, rotation = 0)));
initial equation
//hstart_F_Out = hE.waterOut.h_outflow;
equation
  connect(source_hot.flange, hE.gasIn) annotation(
    Line(points = {{-1.83697e-015, 50}, {-1.83697e-015, 20}, {0, 20}}, color = {0, 0, 255}, thickness = 0.5, smooth = Smooth.None));
  
  connect(hE.gasOut, sink_hot.flange) annotation(
    Line(points = {{1.83697e-015, -70}, {1.83697e-015, -56}, {-8.88178e-016, -56}}, thickness = 0.5, color = {0, 0, 255}));
   
  // connect(hE.waterOut, T_waterOut.inlet) annotation(
    //Line(points = {{8.88178e-016, -44}, {8.88178e-016, -20}, {0, -20}}, thickness = 0.5, color = {0, 0, 255}));

  // connect(T_waterOut.outlet, sink_heater_hot.flange) annotation(
    // Line(points = {{1.83697e-015, -70}, {1.83697e-015, -56}, {-8.88178e-016, -56}}, thickness = 0.5, color = {0, 0, 255}));
  
  connect(source_cold.flange, hE.waterIn) annotation(
    Line(points = {{-50, 0}, {-20, 0}}, color = {159, 159, 223}, thickness = 0.5, smooth = Smooth.None));
  connect(hE.waterOut, T_gasOut.inlet) annotation(
    Line(points = {{34, 0}, {34, 0}, {20, 0}}, color = {159, 159, 223}, thickness = 0.5));
  connect(T_gasOut.outlet, sink_cold.flange) annotation(
    Line(points = {{46, 0}, {46, 0}, {60, 0}}, color = {159, 159, 223}, thickness = 0.5));

annotation(
    Diagram(graphics),
    experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-2, Interval = 1),
    __OpenModelica_commandLineOptions = "--matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection -d=initialization,NLSanalyticJacobian,aliasConflicts,bltdump",    
    __OpenModelica_simulationFlags(lv = "LOG_DEBUG,LOG_NLS,LOG_NLS_V,LOG_STATS,LOG_INIT,LOG_STDOUT, -w", outputFormat = "mat", s = "dassl", nls = "homotopy"));
end TestTP_Cooler;
