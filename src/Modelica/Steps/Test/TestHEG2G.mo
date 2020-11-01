within Steps.Test;

model TestHEG2G
  "Test for HEG2G"  
    
  import Modelica.SIunits.Conversions.{from_degC, from_deg};
  import Modelica.SIunits.{Temperature, Pressure, SpecificEnthalpy};
  import Util = Utilities.Util;
  import Steps.Utilities.CoolProp.PropsSI;  
  import Steps.Components.{PCHEGeoParam, PCHEBoundaryCondition};
  import Steps.Model.{PBConfiguration, SimParam, EntityConfig, EntityGeoParam, EntityThermoParam, ThermoState} ;
  import Model.PBConfiguration;
  import ThermoPower.Choices.Init.Options;
  import ThermoPower.System;
  
  package medium_hot = Steps.Media.CO2;
  package medium_cold = Steps.Media.CO2;
  // package medium_cold = Modelica.Media.IdealGases.MixtureGases.CombustionAir;   
  // package medium_cold = Modelica.Media.IdealGases.SingleGases.O2;

  parameter Model.PBConfiguration cfg_def( 
  p_pump_in = 9e6,
  p_pump_out = 20e6,  
  mdot_main = 125,
  mdot_pump = cfg_def.mdot_main * 0.7,    
  bc_HTR.st_hot_in.T = 883,
  bc_HTR.st_hot_out.T = 643,
  bc_HTR.st_cold_in.T = 633,
  bc_HTR.st_cold_out.T = 843);  
  
  // parameter for C++ implementation of PCHE - based on Modelica impl's result
  parameter Model.PBConfiguration cfg_on_design( 
  p_pump_in = 8e6,
  p_pump_out = 20e6,  
  mdot_main = 51.51,
  mdot_pump = 31.31,
  T_HTR_cold_in = from_degC(141.3), 
  T_HTR_cold_out = from_degC(495.302),
  T_HTR_hot_out = from_degC(141.041),
  T_LTR_hot_out = from_degC(63.6726)); 
  
  parameter Model.PBConfiguration cfg_sscar_model_1_3( 
  p_pump_in = 8e6,
  p_pump_out = 20e6,  
  mdot_main = 125,
  mdot_pump = 68.75,    
  bc_HTR.st_hot_in.T = 883,
  bc_HTR.st_hot_out.T = 643,
  bc_HTR.st_cold_in.T = 633,
  bc_HTR.st_cold_out.T = 843);  
  
  // select the configuration of parameters
  parameter Model.PBConfiguration cfg = cfg_on_design;
  
  // set the values of parameters accordingly
  parameter PCHEBoundaryCondition bc_HTR = cfg.bc_HTR;  
  parameter PCHEBoundaryCondition bc_LTR = cfg.bc_LTR;
  parameter ThermoState st_bypass = cfg.bc_bypass;
  
  parameter EntityGeoParam geo_hot = cfg.cfg_HTR_hot.geo;
  parameter EntityGeoParam geo_cold = cfg.cfg_HTR_cold.geo;
  parameter EntityGeoParam geo_tube = cfg.cfg_HTR_tube.geo;
  parameter EntityGeoParam geo_mixer = cfg.cfg_mixer.geo;
  
  
  parameter EntityThermoParam thermo_hot = cfg.cfg_HTR_hot.thermo;
  parameter EntityThermoParam thermo_cold = cfg.cfg_HTR_cold.thermo;
  parameter EntityThermoParam thermo_tube = cfg.cfg_HTR_tube.thermo;  
  parameter EntityThermoParam thermo_mixer = cfg.cfg_mixer.thermo;

  //parameter Modelica.SIunits.SpecificEnthalpy hstart_F_In = medium_cold.specificEnthalpy_pT(fluidNomPressure, bc_HTR.st_cold_in.T) "Nominal specific enthalpy";
  //parameter Modelica.SIunits.SpecificEnthalpy hstart_F_Out = medium_cold.specificEnthalpy_pT(fluidNomPressure, bc_HTR.st_cold_out.T) "Nominal specific enthalpy";
  //Components
  inner ThermoPower.System system(allowFlowReversal = false, initOpt=ThermoPower.Choices.Init.Options.fixedState) annotation(
    Placement(transformation(extent = {{80, 80}, {100, 100}})));
    
  ThermoPower.Gas.SourceMassFlow sourceW_water(
    redeclare package Medium = medium_cold, 
    T = bc_LTR.st_cold_out.T, 
    p0 = bc_LTR.st_cold_out.p, 
    use_in_T = false, 
    w0 = bc_LTR.st_cold_out.mdot) 
  annotation(
    Placement(transformation(origin = {0, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 270)));
  
  ThermoPower.Gas.SourceMassFlow source_mixer_in(
    redeclare package Medium = medium_cold,
    T = st_bypass.T,
    p0 = st_bypass.p,
    use_in_T = false,
    w0 = st_bypass.mdot    
  );
  
  ThermoPower.Gas.SinkPressure sinkP_water(
    redeclare package Medium = medium_cold, 
    p0 = bc_HTR.st_cold_out.p, 
    T = bc_HTR.st_cold_out.T) 
  annotation(
    Placement(transformation(origin = {0, -80}, extent = {{-10, -10}, {10, 10}}, rotation = 270)));

  ThermoPower.Gas.SinkPressure sinkP_gas(
    redeclare package Medium = medium_hot,
    T = bc_HTR.st_hot_out.T, 
    p0 = bc_HTR.st_hot_out.p) 
  annotation(
    Placement(transformation(extent = {{60, -10}, {80, 10}}, rotation = 0)));
  
  ThermoPower.Gas.SourceMassFlow sourceW_gas(
    redeclare package Medium = medium_hot, 
    T = bc_HTR.st_hot_in.T, 
    p0 = bc_HTR.st_hot_in.p, 
    w0 = bc_HTR.st_hot_in.mdot) 
  annotation(
    Placement(transformation(extent = {{-70, -10}, {-50, 10}}, rotation = 0))); 
 
  Components.SSMixer mixer(
    redeclare package Medium = medium_cold,
    gamma=thermo_mixer.gamma_he,
    S=geo_mixer.A_ex,
    V=geo_mixer.V,
    pstart=bc_LTR.st_cold_out.p,
    Tstart=bc_LTR.st_cold_out.T,
    Tmstart=bc_LTR.st_cold_out.T - 50
  );  

  ThermoPower.Gas.SensT T_waterOut(
    redeclare package Medium = medium_cold) 
  annotation(
    Placement(transformation(origin = {4, -50}, extent = {{-10, -10}, {10, 10}}, rotation = 270)));
  
  ThermoPower.Gas.SensT T_gasOut(
    redeclare package Medium = medium_hot) 
  annotation(
    Placement(transformation(extent = {{30, -6}, {50, 14}}, rotation = 0)));
 
  Components.HEG2G HTR(
  redeclare package FluidMedium = medium_cold, 
  redeclare package FlueGasMedium = medium_hot, 
  redeclare replaceable model HeatTransfer_F = ThermoPower.Thermal.HeatTransferFV.ConstantHeatTransferCoefficient(gamma = thermo_cold.gamma_he), 
  redeclare replaceable model HeatTransfer_G = ThermoPower.Thermal.HeatTransferFV.ConstantHeatTransferCoefficientTwoGrids(gamma = thermo_hot.gamma_he), 
  redeclare model HeatExchangerTopology = ThermoPower.Thermal.HeatExchangerTopologies.CounterCurrentFlow,  
  N_F = geo_cold.N_seg, 
  N_G = geo_hot.N_seg,   
  SSInit = SSInit, 
  Tstartbar_G = bc_HTR.st_hot_in.T, 
  exchSurface_F = geo_cold.A_ex, 
  exchSurface_G = geo_hot.A_ex, 
  extSurfaceTub = geo_tube.A_ex, 
  fluidNomFlowRate = bc_HTR.st_cold_in.mdot, 
  fluidNomPressure = bc_HTR.st_cold_in.p, 
  fluidVol = geo_cold.V, 
  gasNomFlowRate = bc_HTR.st_hot_in.mdot, 
  gasNomPressure = bc_HTR.st_hot_in.p, 
  gasVol = geo_hot.V, 
  lambda = thermo_tube.lambda, 
  metalVol = geo_tube.V, 
  pstart_F = bc_HTR.st_cold_in.p, 
  rhomcm = thermo_tube.rho_mcm) annotation(
    Placement(transformation(extent = {{-20, -20}, {20, 20}}, rotation = 0)));
  //Start value
  // parameter Modelica.SIunits.Temperature Tstart_G = (bc_HTR.st_hot_in.T + bc_HTR.st_hot_out.T) / 2;
  // parameter Modelica.SIunits.Temperature Tstart_M = (bc_HTR.st_hot_in.T + bc_HTR.st_hot_out.T + bc_HTR.st_cold_in.T + bc_HTR.st_cold_out.T) / 4;

  parameter Boolean SSInit = true "Steady-state initialization";

initial equation
//hstart_F_Out = HTR.waterOut.h_outflow;
equation
/*
  // mixer alone
  connect(source_mixer_in.flange, mixer.in1);
  connect(sourceW_water.flange, mixer.in2);
  connect(mixer.out, sinkP_water.flange);
*/  

/*
  // HTR alone
  connect(sourceW_water.flange, HTR.waterIn) annotation(
    Line(points = {{-1.83697e-015, 50}, {-1.83697e-015, 20}, {0, 20}}, color = {0, 0, 255}, thickness = 0.5, smooth = Smooth.None));
  
  connect(T_waterOut.inlet, HTR.waterOut) annotation(
    Line(points = {{8.88178e-016, -44}, {8.88178e-016, -20}, {0, -20}}, thickness = 0.5, color = {0, 0, 255}));      
   
  connect(sinkP_water.flange, T_waterOut.outlet) annotation(
    Line(points = {{1.83697e-015, -70}, {1.83697e-015, -56}, {-8.88178e-016, -56}}, thickness = 0.5, color = {0, 0, 255}));
  
  connect(sourceW_gas.flange, HTR.gasIn) annotation(
        Line(points = {{-50, 0}, {-20, 0}}, color = {159, 159, 223}, thickness = 0.5, smooth = Smooth.None));
    
  connect(T_gasOut.inlet, HTR.gasOut) annotation(
    Line(points = {{34, 0}, {34, 0}, {20, 0}}, color = {159, 159, 223}, thickness = 0.5));
    
  connect(T_gasOut.outlet, sinkP_gas.flange) annotation(
    Line(points = {{46, 0}, {46, 0}, {60, 0}}, color = {159, 159, 223}, thickness = 0.5));    
*/


  //HTR + mixer
  connect(source_mixer_in.flange, mixer.in1);
  
  connect(sourceW_water.flange, mixer.in2);
  
  connect(mixer.out, HTR.waterIn);

  connect(T_gasOut.inlet, HTR.gasOut) annotation(
    Line(points = {{34, 0}, {34, 0}, {20, 0}}, color = {159, 159, 223}, thickness = 0.5));
  connect(T_gasOut.outlet, sinkP_gas.flange) annotation(
    Line(points = {{46, 0}, {46, 0}, {60, 0}}, color = {159, 159, 223}, thickness = 0.5));
  connect(sinkP_water.flange, T_waterOut.outlet) annotation(
    Line(points = {{1.83697e-015, -70}, {1.83697e-015, -56}, {-8.88178e-016, -56}}, thickness = 0.5, color = {0, 0, 255}));
  connect(T_waterOut.inlet, HTR.waterOut) annotation(
    Line(points = {{8.88178e-016, -44}, {8.88178e-016, -20}, {0, -20}}, thickness = 0.5, color = {0, 0, 255}));
  
  connect(sourceW_gas.flange, HTR.gasIn) annotation(
        Line(points = {{-50, 0}, {-20, 0}}, color = {159, 159, 223}, thickness = 0.5, smooth = Smooth.None));

annotation(
    Diagram(graphics),
    experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-3, Interval = 1),
    __OpenModelica_commandLineOptions = "--matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection -d=initialization,NLSanalyticJacobian,aliasConflicts,bltdump",    
    __OpenModelica_simulationFlags(lv = "LOG_DEBUG,LOG_NLS,LOG_NLS_V,LOG_STATS,LOG_INIT,LOG_STDOUT, -w", outputFormat = "mat", s = "dassl", nls = "homotopy"));
end TestHEG2G;
