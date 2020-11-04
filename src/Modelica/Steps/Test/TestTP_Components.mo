within Steps.Test;

model TestTP_Components
  "Test for combination of components in ThermoPower"  
    
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
  package medium_heater = Steps.Media.CO2;// SolarTherm.Media.Sodium;
  //package medium_heater = ThermoPower.Water.StandardWater;// Modelica.Media.IdealGases.SingleGases.CO2;

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

  //Components
  inner ThermoPower.System system(allowFlowReversal = false, initOpt=ThermoPower.Choices.Init.Options.noInit) annotation(
    Placement(transformation(extent = {{80, 80}, {100, 100}})));

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
    p0 = bc_heater.st_hot_out.mdot, 
    T = bc_heater.st_hot_out.T,
    use_T = true) 
    annotation(
    Placement(transformation(origin = {0, -80}, extent = {{-10, -10}, {10, 10}}, rotation = 270)));

  ThermoPower.PowerPlants.HRSG.Components.HE heater(
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
  
  ThermoPower.Gas.SourceMassFlow source_cold(
    redeclare package Medium = medium_cold, 
    T = bc_LTR.st_cold_in.T, 
    p0 = bc_LTR.st_cold_in.p, 
    use_in_T = false, 
    w0 = bc_LTR.st_cold_in.mdot) 
  annotation(
    Placement(transformation(origin = {0, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 270)));
  
  ThermoPower.Gas.SourceMassFlow source_mixer_in(
    redeclare package Medium = medium_cold,
    T = st_bypass.T,
    p0 = st_bypass.p,
    use_in_T = false,
    w0 = st_bypass.mdot    
  );
  /*
  ThermoPower.Gas.SinkPressure sink_cold(
    redeclare package Medium = medium_cold, 
    p0 = bc_heater.st_cold_out.p, 
    T = bc_heater.st_cold_out.T) 
  annotation(
    Placement(transformation(origin = {0, -80}, extent = {{-10, -10}, {10, 10}}, rotation = 270)));
 */
  ThermoPower.Gas.SinkPressure sink_hot(
    redeclare package Medium = medium_hot,
    T = bc_LTR.st_hot_out.T, 
    p0 = bc_LTR.st_hot_out.p) 
  annotation(
    Placement(transformation(extent = {{60, -10}, {80, 10}}, rotation = 0)));
  /*
  ThermoPower.Gas.SourceMassFlow source_hot(
    redeclare package Medium = medium_hot, 
    T = bc_HTR.st_hot_in.T, 
    p0 = bc_HTR.st_hot_in.p, 
    w0 = bc_HTR.st_hot_in.mdot,
    use_in_T = false) 
  annotation(
    Placement(transformation(extent = {{-70, -10}, {-50, 10}}, rotation = 0))); 
  */
  /*
  parameter EntityGeoParam geo_mixer = cfg.cfg_mixer.geo;  
  
  Components.SSMixer mixer(
  //Gas.Mixer mixer(
    redeclare package Medium = medium_cold,
    gamma=thermo_mixer.gamma_he,
    S=geo_mixer.A_ex,
    V=geo_mixer.V,
    pstart=bc_LTR.st_cold_out.p,
    Tstart=bc_LTR.st_cold_out.T,
    Tmstart=bc_LTR.st_cold_out.T - 50
  );
  */  
 
  // use FlowJoin to mix flow
  Gas.FlowJoin mixer(redeclare package Medium = medium_cold);  

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
  Tstartbar_F = bc_HTR.st_cold_in.T, 
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
  pstart_G = bc_HTR.st_hot_in.T,
  rhomcm = thermo_tube.rho_mcm,
  gasQuasiStatic = true,
  fluidQuasiStatic = true) annotation(
    Placement(transformation(extent = {{-20, -20}, {20, 20}}, rotation = 0)));

  Components.HEG2G LTR(
  redeclare package FluidMedium = medium_cold, 
  redeclare package FlueGasMedium = medium_hot, 
  redeclare replaceable model HeatTransfer_F = ThermoPower.Thermal.HeatTransferFV.ConstantHeatTransferCoefficient(gamma = thermo_cold.gamma_he), 
  redeclare replaceable model HeatTransfer_G = ThermoPower.Thermal.HeatTransferFV.ConstantHeatTransferCoefficientTwoGrids(gamma = thermo_hot.gamma_he), 
  redeclare model HeatExchangerTopology = ThermoPower.Thermal.HeatExchangerTopologies.CounterCurrentFlow,  
  N_F = geo_cold.N_seg, 
  N_G = geo_hot.N_seg,   
  SSInit = SSInit, 
  Tstartbar_G = bc_LTR.st_hot_in.T, 
  Tstartbar_F = bc_LTR.st_cold_in.T, 
  exchSurface_F = geo_cold.A_ex, 
  exchSurface_G = geo_hot.A_ex, 
  extSurfaceTub = geo_tube.A_ex, 
  fluidNomFlowRate = bc_LTR.st_cold_in.mdot, 
  fluidNomPressure = bc_LTR.st_cold_in.p, 
  fluidVol = geo_cold.V, 
  gasNomFlowRate = bc_LTR.st_hot_in.mdot, 
  gasNomPressure = bc_LTR.st_hot_in.p, 
  gasVol = geo_hot.V, 
  lambda = thermo_tube.lambda, 
  metalVol = geo_tube.V, 
  pstart_F = bc_LTR.st_cold_in.p, 
  pstart_G = bc_LTR.st_hot_in.T,
  rhomcm = thermo_tube.rho_mcm,
  gasQuasiStatic = true,
  fluidQuasiStatic = true) annotation(
    Placement(transformation(extent = {{-20, -20}, {20, 20}}, rotation = 0)));
/*
  Modelica.Mechanics.Rotational.Components.Inertia Inertia1(J = 10000) annotation(
    Placement(transformation(extent = {{10, -10}, {30, 10}}, rotation = 0)));
*/
  ThermoPower.Gas.Turbine Turbine1(
  redeclare package Medium = medium_hot, 
  tablePhic = tablePhic, 
  tableEta = tableEta, 
  pstart_in = bc_heater.st_cold_out.p, 
  pstart_out = bc_HTR.st_hot_in.p, 
  Tstart_in = bc_heater.st_cold_out.T, 
  Tstart_out = bc_HTR.st_hot_in.T, 
  Ndesign = 60000.0, 
  Tdes_in = bc_heater.st_cold_out.T, 
  Table = ThermoPower.Choices.TurboMachinery.TableTypes.matrix) 
  annotation(
    Placement(transformation(extent = {{-40, -20}, {0, 20}}, rotation = 0)));
 
  ThermoPower.Gas.SensT sens_turbine(redeclare package Medium = medium_hot) annotation(
    Placement(visible = true, transformation(origin = {20, 44}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));


  parameter Boolean SSInit = false "Steady-state initialization";

protected
  parameter Real tablePhic[5, 4] = [1, 37, 80, 100; 1.5, 7.10E-05, 7.10E-05, 7.10E-05; 2, 8.40E-05, 8.40E-05, 8.40E-05; 2.5, 8.70E-05, 8.70E-05, 8.70E-05; 3, 1.04E-04, 1.04E-04, 1.04E-04];
  parameter Real tableEta[5, 4] = [1, 37, 80, 100; 1.5, 0.57, 0.89, 0.81; 2, 0.46, 0.82, 0.88; 2.5, 0.41, 0.76, 0.85; 3, 0.38, 0.72, 0.82];

  
equation
  //HTR + mixer + LTR + Heater + Turbine
  // main stream, water/cold side  
  connect(source_mixer_in.flange, mixer.inlet1);
  
  connect(source_cold.flange, LTR.waterIn);
  
  connect(LTR.waterOut, mixer.inlet2);
  
  connect(mixer.outlet, HTR.waterIn);

  connect(HTR.waterOut, heater.waterIn);
  
  connect(heater.waterOut, T_waterOut.inlet) annotation(
    Line(points = {{8.88178e-016, -44}, {8.88178e-016, -20}, {0, -20}}, thickness = 0.5, color = {0, 0, 255}));

  connect(T_waterOut.outlet, Turbine1.inlet) annotation(
    Line(points = {{1.83697e-015, -70}, {1.83697e-015, -56}, {-8.88178e-016, -56}}, thickness = 0.5, color = {0, 0, 255}));
  
  // main stream, gas/hot side
  connect(Turbine1.outlet, sens_turbine.inlet) annotation(
   Line(points = {{-50, 0}, {-20, 0}}, color = {159, 159, 223}, thickness = 0.5, smooth = Smooth.None));
  
  //connect(Turbine1.shaft_b, Inertia1.flange_a);
  
  connect(sens_turbine.outlet, HTR.gasIn);
  
  connect(HTR.gasOut, LTR.gasIn);
  
  connect(LTR.gasOut, T_gasOut.inlet) annotation(
    Line(points = {{34, 0}, {34, 0}, {20, 0}}, color = {159, 159, 223}, thickness = 0.5));
    
  connect(T_gasOut.outlet, sink_hot.flange) annotation(
    Line(points = {{46, 0}, {46, 0}, {60, 0}}, color = {159, 159, 223}, thickness = 0.5));
   
  // hot stream for heater
  connect(source_heater_hot.flange, heater.gasIn);
  connect(heater.gasOut, sink_heater_hot.flange);
  

/*
  // mixer alone
  connect(source_mixer_in.flange, mixer.in1);
  connect(source_cold.flange, mixer.in2);
  connect(mixer.out, sink_cold.flange);
*/  

/*
  // HTR alone
  connect(source_cold.flange, HTR.waterIn) annotation(
    Line(points = {{-1.83697e-015, 50}, {-1.83697e-015, 20}, {0, 20}}, color = {0, 0, 255}, thickness = 0.5, smooth = Smooth.None));
  
  connect(T_waterOut.inlet, HTR.waterOut) annotation(
    Line(points = {{8.88178e-016, -44}, {8.88178e-016, -20}, {0, -20}}, thickness = 0.5, color = {0, 0, 255}));      
   
  connect(sink_cold.flange, T_waterOut.outlet) annotation(
    Line(points = {{1.83697e-015, -70}, {1.83697e-015, -56}, {-8.88178e-016, -56}}, thickness = 0.5, color = {0, 0, 255}));
  
  connect(source_hot.flange, HTR.gasIn) annotation(
        Line(points = {{-50, 0}, {-20, 0}}, color = {159, 159, 223}, thickness = 0.5, smooth = Smooth.None));
    
  connect(T_gasOut.inlet, HTR.gasOut) annotation(
    Line(points = {{34, 0}, {34, 0}, {20, 0}}, color = {159, 159, 223}, thickness = 0.5));
    
  connect(T_gasOut.outlet, sink_hot.flange) annotation(
    Line(points = {{46, 0}, {46, 0}, {60, 0}}, color = {159, 159, 223}, thickness = 0.5));    

  // LTR alone
  // water/cold side  
  connect(source_cold.flange, LTR.waterIn);
  
  connect(LTR.waterOut, sink_cold.flange) annotation(
    Line(points = {{1.83697e-015, -70}, {1.83697e-015, -56}, {-8.88178e-016, -56}}, thickness = 0.5, color = {0, 0, 255}));
  
  // gas/hot side
  connect(source_hot.flange, LTR.gasIn);
  
  connect(LTR.gasOut, sink_hot.flange) annotation(
    Line(points = {{46, 0}, {46, 0}, {60, 0}}, color = {159, 159, 223}, thickness = 0.5));
  */  
 
/*
  //HTR + mixer
  connect(source_mixer_in.flange, mixer.in1);
  
  connect(source_cold.flange, mixer.in2);
  
  connect(mixer.out, HTR.waterIn);

  connect(T_gasOut.inlet, HTR.gasOut) annotation(
    Line(points = {{34, 0}, {34, 0}, {20, 0}}, color = {159, 159, 223}, thickness = 0.5));
  connect(T_gasOut.outlet, sink_hot.flange) annotation(
    Line(points = {{46, 0}, {46, 0}, {60, 0}}, color = {159, 159, 223}, thickness = 0.5));
  connect(sink_cold.flange, T_waterOut.outlet) annotation(
    Line(points = {{1.83697e-015, -70}, {1.83697e-015, -56}, {-8.88178e-016, -56}}, thickness = 0.5, color = {0, 0, 255}));
  connect(T_waterOut.inlet, HTR.waterOut) annotation(
    Line(points = {{8.88178e-016, -44}, {8.88178e-016, -20}, {0, -20}}, thickness = 0.5, color = {0, 0, 255}));
  
  connect(source_hot.flange, HTR.gasIn) annotation(
   Line(points = {{-50, 0}, {-20, 0}}, color = {159, 159, 223}, thickness = 0.5, smooth = Smooth.None));
*/
/*
  // mixer + LTR
  // water/cold side  
  connect(source_mixer_in.flange, mixer.in1);
  
  connect(source_cold.flange, LTR.waterIn);
  
  connect(LTR.waterOut, mixer.in2);
  
  connect(mixer.out, sink_cold.flange) annotation(
    Line(points = {{1.83697e-015, -70}, {1.83697e-015, -56}, {-8.88178e-016, -56}}, thickness = 0.5, color = {0, 0, 255}));
  
  // gas/hot side
  connect(source_hot.flange, LTR.gasIn);
  
  connect(LTR.gasOut, sink_hot.flange) annotation(
    Line(points = {{46, 0}, {46, 0}, {60, 0}}, color = {159, 159, 223}, thickness = 0.5));
*/
/*
  //HTR + mixer + LTR 
  // water/cold side  
  connect(source_mixer_in.flange, mixer.inlet1);
  
  connect(source_cold.flange, LTR.waterIn);
  
  connect(LTR.waterOut, mixer.inlet2);
  
  connect(mixer.outlet, HTR.waterIn);

  connect(HTR.waterOut, T_waterOut.inlet) annotation(
    Line(points = {{8.88178e-016, -44}, {8.88178e-016, -20}, {0, -20}}, thickness = 0.5, color = {0, 0, 255}));


  connect(T_waterOut.outlet, sink_cold.flange) annotation(
    Line(points = {{1.83697e-015, -70}, {1.83697e-015, -56}, {-8.88178e-016, -56}}, thickness = 0.5, color = {0, 0, 255}));
  
  // gas/hot side
  connect(source_hot.flange, HTR.gasIn) annotation(
   Line(points = {{-50, 0}, {-20, 0}}, color = {159, 159, 223}, thickness = 0.5, smooth = Smooth.None));

  connect(HTR.gasOut, LTR.gasIn);
  
  connect(LTR.gasOut, T_gasOut.inlet) annotation(
    Line(points = {{34, 0}, {34, 0}, {20, 0}}, color = {159, 159, 223}, thickness = 0.5));
    
  connect(T_gasOut.outlet, sink_hot.flange) annotation(
    Line(points = {{46, 0}, {46, 0}, {60, 0}}, color = {159, 159, 223}, thickness = 0.5));
*/
/*
  //HTR + mixer + LTR + Heater
  // main stream, water/cold side  
  connect(source_mixer_in.flange, mixer.inlet1);
  
  connect(source_cold.flange, LTR.waterIn);
  
  connect(LTR.waterOut, mixer.inlet2);
  
  connect(mixer.outlet, HTR.waterIn);

  connect(HTR.waterOut, heater.waterIn);
  
  connect(heater.waterOut, T_waterOut.inlet) annotation(
    Line(points = {{8.88178e-016, -44}, {8.88178e-016, -20}, {0, -20}}, thickness = 0.5, color = {0, 0, 255}));

  connect(T_waterOut.outlet, sink_cold.flange) annotation(
    Line(points = {{1.83697e-015, -70}, {1.83697e-015, -56}, {-8.88178e-016, -56}}, thickness = 0.5, color = {0, 0, 255}));
  
  // main stream, gas/hot side
  connect(source_hot.flange, HTR.gasIn) annotation(
   Line(points = {{-50, 0}, {-20, 0}}, color = {159, 159, 223}, thickness = 0.5, smooth = Smooth.None));

  connect(HTR.gasOut, LTR.gasIn);
  
  connect(LTR.gasOut, T_gasOut.inlet) annotation(
    Line(points = {{34, 0}, {34, 0}, {20, 0}}, color = {159, 159, 223}, thickness = 0.5));
    
  connect(T_gasOut.outlet, sink_hot.flange) annotation(
    Line(points = {{46, 0}, {46, 0}, {60, 0}}, color = {159, 159, 223}, thickness = 0.5));
   
  // hot stream for heater
  connect(source_heater_hot.flange, heater.gasIn);
  connect(heater.gasOut, sink_heater_hot.flange);
*/



  
annotation(
    Diagram(graphics),
    experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-3, Interval = 1),
    __OpenModelica_commandLineOptions = "--matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection -d=initialization,NLSanalyticJacobian,aliasConflicts,bltdump",    
    __OpenModelica_simulationFlags(lv = "LOG_DEBUG,LOG_NLS,LOG_NLS_V,LOG_STATS,LOG_INIT,LOG_STDOUT, -w", outputFormat = "mat", s = "dassl", nls = "homotopy"));
end TestTP_Components;
