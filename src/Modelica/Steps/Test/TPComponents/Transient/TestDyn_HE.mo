within Steps.Test.TPComponents.Transient;

model TestDyn_HE
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
  
  
  package medium_heater= SolarTherm.Media.Sodium.Sodium_pT;   
  // package medium_heater = Steps.Media.MoltenSalt.MoltenSalt_pT;
  //package medium_heater = ThermoPower.Water.StandardWater;
  // package medium_heater = Steps.Media.ThermiaOilD;//
  // package medium_cold = Steps.Media.CO2; // Modelica.Media.IdealGases.SingleGases.CO2;
  package medium_hot = Steps.Media.SCO2;//ExternalMedia.Examples.CO2CoolProp;
  package medium_cold = Steps.Media.SCO2;//ExternalMedia.Examples.CO2CoolProp;    
  
  // package medium_hot = ExternalMedia.Examples.CO2CoolProp;
  // package medium_cold = ExternalMedia.Examples.CO2CoolProp;    
  
  // package medium_hot = Modelica.Media.IdealGases.SingleGases.CO2;
  // package medium_cold = Modelica.Media.IdealGases.SingleGases.CO2;
  
  //gas
  parameter Modelica.SIunits.MassFlowRate gasNomFlowRate = 125 "Nominal mass flowrate";
  parameter Modelica.SIunits.Pressure gasNomPressure = 9e6 "Nominal pressure in the gas side inlet";
  parameter Modelica.SIunits.Temperature Tstart_G_In = 883 "Inlet gas temperature start value";
  parameter Modelica.SIunits.Temperature Tstart_G_Out = 643 "Outlet gas temperature start value";
  //fluid
  parameter Modelica.SIunits.MassFlowRate fluidNomFlowRate = 125 "Nominal flow rate through the fluid side";
  parameter Modelica.SIunits.Pressure fluidNomPressure = 9e6 "Nominal pressure in the fluid side inlet";
  parameter Modelica.SIunits.CoefficientOfHeatTransfer gamma_G = 200 "Constant heat transfer coefficient in the gas side";
  parameter Modelica.SIunits.CoefficientOfHeatTransfer gamma_F = 200 "Constant heat transfer coefficient in the fluid side";
  parameter Modelica.SIunits.Temperature Tstart_M_In = Tstart_F_In "Inlet metal wall temperature start value";
  parameter Modelica.SIunits.Temperature Tstart_M_Out = Tstart_F_Out "Outlet metal wall temperature start value";
  parameter Modelica.SIunits.Temperature Tstart_F_In = 633 "Inlet fluid temperature start value";
  parameter Modelica.SIunits.Temperature Tstart_F_Out = 843 "Outlet fluid temperature start value";

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
  parameter Model.PBConfiguration cfg = cfg_test;
  
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

  Steps.TPComponents.HE hE(
  //Components.HEG2G hE(
    redeclare package FluidMedium = medium_heater, 
    redeclare package FlueGasMedium = medium_cold, 
    fluidFlow(fixedMassFlowSimplified = true, hstartin = bc_heater.st_hot_in.h, hstartout=bc_heater.st_hot_out.h), // set the fluid flow as fixed mdot for simplarity
    gasFlow(Tstartin = bc_heater.st_cold_in.T, Tstartout = bc_heater.st_cold_out.T),
    bc = bc_heater, 
    geo_hot = cfg.cfg_heater_hot.geo,
    geo_cold = cfg.cfg_heater_cold.geo,
    geo_tube = cfg.cfg_heater_tube.geo,  
    thermo_hot = cfg.cfg_heater_hot.thermo,
    thermo_cold = cfg.cfg_heater_cold.thermo,
    thermo_tube = cfg.cfg_heater_tube.thermo, 
    SSInit=true,
    gasQuasiStatic = false,
    FluidPhaseStart=ThermoPower.Choices.FluidPhase.FluidPhases.Liquid,    
    // redeclare replaceable model HeatTransfer_F = ThermoPower.Thermal.HeatTransferFV.IdealHeatTransfer, 
    redeclare replaceable model HeatTransfer_F = ThermoPower.Thermal.HeatTransferFV.ConstantHeatTransferCoefficientTwoGrids(gamma = thermo_hot.gamma_he),     
    // redeclare replaceable model HeatTransfer_G = ThermoPower.Thermal.HeatTransferFV.IdealHeatTransfer,      
    redeclare replaceable model HeatTransfer_G = ThermoPower.Thermal.HeatTransferFV.ConstantHeatTransferCoefficient(gamma =  thermo_cold.gamma_he),     
    redeclare model HeatExchangerTopology = ThermoPower.Thermal.HeatExchangerTopologies.CounterCurrentFlow,
    metalTube(WallRes=false, Tstartbar=bc_heater.st_hot_in.T - 50)) annotation(
    Placement(transformation(extent = {{-20, -20}, {20, 20}}, rotation = 0)));

  /* 
  // ThermoPower.PowerPlants.HRSG.Components.HE hE(
  Steps.TPComponents.HE hE(
    redeclare package FluidMedium = medium_cold, 
    redeclare package FlueGasMedium = medium_hot, 
    redeclare replaceable model HeatTransfer_F = ThermoPower.Thermal.HeatTransferFV.ConstantHeatTransferCoefficient(gamma = gamma_F), 
    redeclare replaceable model HeatTransfer_G = ThermoPower.Thermal.HeatTransferFV.ConstantHeatTransferCoefficientTwoGrids(gamma = gamma_G), 
    redeclare model HeatExchangerTopology = ThermoPower.Thermal.HeatExchangerTopologies.CounterCurrentFlow, 
    N_F = 7, 
    N_G = 7,
    Nw_G = 6, 
    SSInit = true, 
    gasQuasiStatic = false,
    Tstartbar_G = 823.15, 
    Tstartbar_M = 773.15, 
    exchSurface_F = 225.073, 
    exchSurface_G = 1708.2, 
    extSurfaceTub = 252.286, 
    fluidNomFlowRate = fluidNomFlowRate, 
    fluidNomPressure = fluidNomPressure, 
    fluidVol = 2.234, 
    gasNomFlowRate = gasNomFlowRate, 
    gasNomPressure = gasNomPressure,     
    gasVol = 10, 
    lambda = 20, 
    metalVol = 0.573, 
    pstart_F = fluidNomPressure, 
    rhomcm = 7900 * 578.05) annotation(
            Placement(transformation(extent = {{-20, -20}, {20, 20}}, rotation = 0)));
 */
 
  inner ThermoPower.System system(allowFlowReversal = false, initOpt=ThermoPower.Choices.Init.Options.steadyState) annotation(
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
    __OpenModelica_commandLineOptions = "--matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection -d=initialization,NLSanalyticJacobian,aliasConflicts",    
    __OpenModelica_simulationFlags(lv = "LOG_DEBUG,LOG_NLS,LOG_NLS_V,LOG_STATS,LOG_INIT,LOG_STDOUT, -w", outputFormat = "mat", s = "dassl", nls = "homotopy"));
end TestDyn_HE;
