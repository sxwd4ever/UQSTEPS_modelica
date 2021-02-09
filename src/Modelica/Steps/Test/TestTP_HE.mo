within Steps.Test;

model TestTP_HE
  "Test against HE as a Heater in ThermoPower"  
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
  
  // package medium_hot = Steps.Media.CO2;
  package medium_heater_v1 = SolarTherm.Media.Sodium.Sodium_pT; // ThermoPower.Water.StandardWater; //Modelica.Media.IdealGases.SingleGases.CO2;
  package medium_heater = Media.MoltenSalt.MoltenSalt_pT;
  package medium_heater_v3 = Steps.Media.ThermiaOilD;//
  //package medium_cold = Steps.Media.CO2; // Modelica.Media.IdealGases.SingleGases.CO2;
  package medium_hot = Steps.Media.SCO2;//ExternalMedia.Examples.CO2CoolProp;
  package medium_cold = Steps.Media.SCO2;//ExternalMedia.Examples.CO2CoolProp;    
  
  parameter Model.PBConfiguration cfg_tune( 
  redeclare package medium_heater_hot = medium_heater,
  redeclare package medium_heater_cold = medium_hot,
  //mdot_heater = 40,
  //T_heater_hot_in = from_degC(800),
  //T_heater_hot_out = from_degC(600),
  r_i_h = 20e-3,
  r_t_h = cfg_tune.r_i_h + 10e-3,
  r_o_h = 1/2, // agree with the final parameter Dhyd = 1 in HE, should be checked to see if it is capable of containing all fluid-metal tubes
  N_ch_h = 100,
  L_h = 1); 
  
  parameter Model.PBConfiguration cfg_test( 
  redeclare package medium_heater_hot = medium_heater,
  redeclare package medium_heater_cold = medium_hot,  
  mdot_main = 93.75,
  mdot_heater = 55,
  T_heater_hot_in = from_degC(550),
  L_h = 1); 
  
  
  // select the configuration of parameters
  parameter Model.PBConfiguration cfg = cfg_tune;
  
  // set the values of parameters accordingly
  parameter HEBoundaryCondition bc_heater = cfg.bc_heater;
  
  ThermoPower.Water.SourceMassFlow source_heater_hot(
  redeclare package Medium = medium_heater,
    w0 = bc_heater.st_hot_in.mdot,
    p0 = bc_heater.st_hot_in.p,
    h = bc_heater.st_hot_in.h,
    T = bc_heater.st_hot_in.T) //,
    // use_T = true,
    // use_in_T = false) 
    annotation(
    Placement(transformation(origin = {0, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 270)));
    
  ThermoPower.Water.SinkPressure sink_heater_hot(
  redeclare package Medium = medium_heater,
    p0 = bc_heater.st_hot_out.p,
    T = bc_heater.st_hot_out.T,
    h = bc_heater.st_hot_out.h,
    use_T = true) 
    annotation(
    Placement(transformation(origin = {0, -80}, extent = {{-10, -10}, {10, 10}}, rotation = 270)));
    
  ThermoPower.Gas.SourceMassFlow source_cold(
  redeclare package Medium = medium_cold,
    T = bc_heater.st_cold_in.T, 
    p0 = bc_heater.st_cold_in.p, 
    use_in_T = false, 
    w0 = bc_heater.st_cold_in.mdot,
    gas(p(nominal = bc_heater.st_cold_in.p), 
    T(nominal=bc_heater.st_cold_in.T)))
  annotation(
    Placement(transformation(origin = {0, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 270)));      
  
  ThermoPower.Gas.SinkPressure sink_cold(
  redeclare package Medium = medium_cold,
    p0 = bc_heater.st_cold_out.p, 
    T = bc_heater.st_cold_out.T) 
  annotation(
    Placement(transformation(origin = {0, -80}, extent = {{-10, -10}, {10, 10}}, rotation = 270)));  
  
  TPComponents.HE hE(
  //Components.HEG2G hE(
    redeclare package FluidMedium = medium_heater, 
    redeclare package FlueGasMedium = medium_cold, 
    fluidFlow(fixedMassFlowSimplified = true, hstartin = bc_heater.st_hot_in.h, hstartout=bc_heater.st_hot_out.h), // set the fluid flow as fixed mdot for simplarity
    gasFlow(QuasiStatic = true, Tstartin = bc_heater.st_cold_in.T, Tstartout = bc_heater.st_cold_out.T),
    bc = bc_heater, 
    geo_hot = cfg.cfg_heater_hot.geo,
    geo_cold = cfg.cfg_heater_cold.geo,
    geo_tube = cfg.cfg_heater_tube.geo,  
    thermo_hot = cfg.cfg_heater_hot.thermo,
    thermo_cold = cfg.cfg_heater_cold.thermo,
    thermo_tube = cfg.cfg_heater_tube.thermo, 
    SSInit=true,
    FluidPhaseStart=ThermoPower.Choices.FluidPhase.FluidPhases.Liquid,    
    redeclare replaceable model HeatTransfer_F = ThermoPower.Thermal.HeatTransferFV.IdealHeatTransfer, //ConstantHeatTransferCoefficientTwoGrids(gamma = thermo_hot.gamma_he),     
    redeclare replaceable model HeatTransfer_G = ThermoPower.Thermal.HeatTransferFV.IdealHeatTransfer, //ConstantHeatTransferCoefficient(gamma =  thermo_cold.gamma_he),     
    redeclare model HeatExchangerTopology = ThermoPower.Thermal.HeatExchangerTopologies.CounterCurrentFlow,
    metalTube(WallRes=false, Tstartbar=bc_heater.st_hot_in.T - 50)) annotation(
    Placement(transformation(extent = {{-20, -20}, {20, 20}}, rotation = 0)));

  inner ThermoPower.System system(allowFlowReversal = false, initOpt=ThermoPower.Choices.Init.Options.noInit) annotation(
    Placement(transformation(extent = {{80, 80}, {100, 100}})));
  
  ThermoPower.Water.SensT T_waterIn(redeclare package Medium = medium_heater) annotation(
    Placement(transformation(origin = {4, -50}, extent = {{-10, -10}, {10, 10}}, rotation = 270)));
  ThermoPower.Water.SensT T_waterOut(redeclare package Medium = medium_heater) annotation(
    Placement(transformation(origin = {4, -50}, extent = {{-10, -10}, {10, 10}}, rotation = 270)));
  
  ThermoPower.Gas.SensT T_gasIn(redeclare package Medium = medium_cold) annotation(
    Placement(transformation(extent = {{30, -6}, {50, 14}}, rotation = 0)));
  ThermoPower.Gas.SensT T_gasOut(redeclare package Medium = medium_cold) annotation(
    Placement(transformation(extent = {{30, -6}, {50, 14}}, rotation = 0)));
  
  // variable for validation
  Modelica.SIunits.Power Q_out = (hE.gasIn.h_outflow - hE.gasOut.h_outflow) * hE.gasIn.m_flow; 
  Modelica.SIunits.Power Q_in = (hE.waterOut.h_outflow - hE.waterIn.h_outflow) * hE.waterIn.m_flow;
  Boolean isQMatch = abs(Q_out -Q_in) < 1e-3;
initial equation
//hstart_F_Out = hE.waterOut.h_outflow;
equation
  connect(source_heater_hot.flange, T_waterIn.inlet);
  connect(T_waterIn.outlet, hE.waterIn) annotation(
    Line(points = {{-1.83697e-015, 50}, {-1.83697e-015, 20}, {0, 20}}, color = {0, 0, 255}, thickness = 0.5, smooth = Smooth.None));
  connect(hE.waterOut, T_waterOut.inlet);
  connect(T_waterOut.outlet, sink_heater_hot.flange) annotation(
    Line(points = {{1.83697e-015, -70}, {1.83697e-015, -56}, {-8.88178e-016, -56}}, thickness = 0.5, color = {0, 0, 255}));
   
  connect(source_cold.flange, T_gasIn.inlet);
  connect(T_gasIn.outlet, hE.gasIn) annotation(
    Line(points = {{-50, 0}, {-20, 0}}, color = {159, 159, 223}, thickness = 0.5, smooth = Smooth.None));
  connect(hE.gasOut, T_gasOut.inlet);
  connect(T_gasOut.outlet, sink_cold.flange) annotation(
    Line(points = {{46, 0}, {46, 0}, {60, 0}}, color = {159, 159, 223}, thickness = 0.5));

annotation(
    Diagram(graphics),
    experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-2, Interval = 1),
    __OpenModelica_commandLineOptions = "--matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection -d=initialization,NLSanalyticJacobian,aliasConflicts,bltdump",    
    __OpenModelica_simulationFlags(lv = "LOG_DEBUG,LOG_NLS,LOG_NLS_V,LOG_STATS,LOG_INIT,LOG_STDOUT, -w", outputFormat = "mat", s = "dassl", nls = "homotopy"));
end TestTP_HE;
