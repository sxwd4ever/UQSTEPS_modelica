within Steps.Test.TPComponents;

model TestTP_HEG2G
  "Test for Gas to Gas Heat Exchanger in ThermoPower"  
    
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
  // package medium_cold = Steps.Media.CO2;
  // package medium_hot = ExternalMedia.Examples.CO2CoolProp;
  // package medium_cold = ExternalMedia.Examples.CO2CoolProp;  
  package medium_hot = Steps.Media.SCO2;
  package medium_cold = Steps.Media.SCO2; 
  // package medium_hot = Modelica.Media.IdealGases.SingleGases.CO2;
  // package medium_cold = Modelica.Media.IdealGases.SingleGases.CO2;  
 
  // parameter for C++ implementation of PCHE - based on Modelica impl's result
  parameter Model.PBConfiguration cfg_default;
  
  parameter Model.PBConfiguration cfg_tune( 
  //mdot_main = 125,
  r_i_HTR = 5e-3,
  r_t_HTR = cfg_tune.r_i_HTR + 5e-3,
  r_o_HTR = cfg_tune.r_t_HTR + 30e-3,
  N_ch_HTR = 100,
  L_HTR = 1 "Don't modify this, since L in HE model is fixed as 1m. Modify Nt instead",
  r_i_LTR = 100e-3,
  r_t_LTR = cfg_tune.r_i_LTR + 10e-3,
  r_o_LTR = cfg_tune.r_t_LTR + 60e-3,
  N_ch_LTR = 400,
  L_LTR = 1 "Don't modify this, since L in HE model is fixed as 1m. Modify Nt instead");   
  
  // select the configuration of parameters
  parameter Model.PBConfiguration cfg = cfg_tune;  
/* 
  // set the values of parameters accordingly - For HTR test
  parameter HEBoundaryCondition bc_HE = cfg.bc_HTR; 
   
  parameter EntityGeoParam geo_hot = cfg.cfg_HTR_hot.geo;
  parameter EntityGeoParam geo_cold = cfg.cfg_HTR_cold.geo;
  parameter EntityGeoParam geo_tube = cfg.cfg_HTR_tube.geo;
  
  parameter EntityThermoParam thermo_hot = cfg.cfg_HTR_hot.thermo;
  parameter EntityThermoParam thermo_cold = cfg.cfg_HTR_cold.thermo;
  parameter EntityThermoParam thermo_tube = cfg.cfg_HTR_tube.thermo;  
*/

  // set the values of parameters accordingly - For LTR Test
  parameter HEBoundaryCondition bc_HE = cfg.bc_LTR; 
  
  parameter EntityGeoParam geo_hot = cfg.cfg_LTR_hot.geo;
  parameter EntityGeoParam geo_cold = cfg.cfg_LTR_cold.geo;
  parameter EntityGeoParam geo_tube = cfg.cfg_LTR_tube.geo;
  
  parameter EntityThermoParam thermo_hot = cfg.cfg_LTR_hot.thermo;
  parameter EntityThermoParam thermo_cold = cfg.cfg_LTR_cold.thermo;
  parameter EntityThermoParam thermo_tube = cfg.cfg_LTR_tube.thermo;  

  //Components
  inner ThermoPower.System system(allowFlowReversal = false, initOpt=ThermoPower.Choices.Init.Options.noInit) annotation(
    Placement(transformation(extent = {{80, 80}, {100, 100}})));  
  
  ThermoPower.Gas.SourceMassFlow source_cold(
    redeclare package Medium = medium_cold, 
    T = bc_HE.st_cold_in.T, 
    p0 = bc_HE.st_cold_in.p, 
    //use_in_T = false, 
    w0 = bc_HE.st_cold_in.mdot,
    gas(p(nominal = bc_HE.st_cold_in.p), 
    T(nominal=bc_HE.st_cold_in.T))) 
  annotation(
    Placement(transformation(origin = {0, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 270)));
  
  ThermoPower.Gas.SinkPressure sink_cold(
    redeclare package Medium = medium_cold, 
    p0 = bc_HE.st_cold_out.p, 
    T = bc_HE.st_cold_out.T) 
  annotation(
    Placement(transformation(origin = {0, -80}, extent = {{-10, -10}, {10, 10}}, rotation = 270)));
 
  ThermoPower.Gas.SinkPressure sink_hot(
    redeclare package Medium = medium_hot,
    T = bc_HE.st_hot_out.T, 
    p0 = bc_HE.st_hot_out.p) 
  annotation(
    Placement(transformation(extent = {{60, -10}, {80, 10}}, rotation = 0)));
  
  ThermoPower.Gas.SourceMassFlow source_hot(
    redeclare package Medium = medium_hot, 
    T = bc_HE.st_hot_in.T, 
    p0 = bc_HE.st_hot_in.p, 
    w0 = bc_HE.st_hot_in.mdot,
    //use_in_T = false,
    gas(p(nominal = bc_HE.st_hot_in.p), 
    T(nominal=bc_HE.st_hot_in.T))) 
  annotation(
    Placement(transformation(extent = {{-70, -10}, {-50, 10}}, rotation = 0))); 
  /*  
  ThermoPower.Gas.SensT T_waterIn(redeclare package Medium = medium_cold) annotation(
    Placement(transformation(origin = {4, -50}, extent = {{-10, -10}, {10, 10}}, rotation = 270)));
  ThermoPower.Gas.SensT T_waterOut(redeclare package Medium = medium_cold) annotation(
    Placement(transformation(origin = {4, -50}, extent = {{-10, -10}, {10, 10}}, rotation = 270)));
  
  ThermoPower.Gas.SensT T_gasIn(redeclare package Medium = medium_hot);
  ThermoPower.Gas.SensT T_gasOut(redeclare package Medium = medium_hot);
 */
  TPComponents.HEG2G HE(
  redeclare package FluidMedium = medium_cold, 
  redeclare package FlueGasMedium = medium_hot, 
  redeclare replaceable model HeatTransfer_F =
   ThermoPower.Thermal.HeatTransferFV.FlowDependentThermalConductance(
    UAnom = thermo_cold.gamma_he * geo_cold.A_ex * HE.Nt,    
    alpha = 0.8
   ), 
   //ThermoPower.Thermal.HeatTransferFV.ConstantHeatTransferCoefficient(gamma = thermo_cold.gamma_he), 
   //ThermoPower.Thermal.HeatTransferFV.IdealHeatTransfer, 
   //ThermoPower.Thermal.HeatTransferFV.ConstantHeatTransferCoefficient(gamma = thermo_cold.gamma_he), 
  redeclare replaceable model HeatTransfer_G =
   ThermoPower.Thermal.HeatTransferFV.FlowDependentThermalConductance(
    UAnom = thermo_hot.gamma_he * geo_hot.A_ex * HE.Nt,
    alpha = 0.8
   ),
   //ThermoPower.Thermal.HeatTransferFV.ConstantHeatTransferCoefficientTwoGrids(gamma = thermo_hot.gamma_he), 
   // ThermoPower.Thermal.HeatTransferFV.IdealHeatTransfer, 
   //ThermoPower.Thermal.HeatTransferFV.ConstantHeatTransferCoefficientTwoGrids(gamma = thermo_hot.gamma_he), 
  redeclare model HeatExchangerTopology = ThermoPower.Thermal.HeatExchangerTopologies.CounterCurrentFlow,  
  gasFlow(QuasiStatic = true, Tstartin = bc_HE.st_hot_in.T, Tstartout = bc_HE.st_hot_out.T),
  fluidFlow(Tstartin = bc_HE.st_cold_in.T, Tstartout = bc_HE.st_cold_out.T),
  N_F = geo_cold.N_seg, 
  N_G = geo_hot.N_seg,  
  Nt = geo_hot.N_ch,  
  Tstartbar_G = bc_HE.st_hot_in.T, 
  Tstartbar_F = bc_HE.st_cold_in.T, 
  exchSurface_F = geo_cold.A_ex, 
  exchSurface_G = geo_hot.A_ex, 
  extSurfaceTub = geo_tube.A_ex, 
  fluidNomFlowRate = bc_HE.st_cold_in.mdot, 
  fluidNomPressure = bc_HE.st_cold_in.p, 
  fluidVol = geo_cold.V, 
  gasNomFlowRate = bc_HE.st_hot_in.mdot, 
  gasNomPressure = bc_HE.st_hot_in.p, 
  gasVol = geo_hot.V, 
  lambda = thermo_tube.lambda, 
  metalVol = geo_tube.V, 
  pstart_F = bc_HE.st_cold_in.p, 
  pstart_G = bc_HE.st_hot_in.p,
  rhomcm = thermo_tube.rho_mcm,
  SSInit = false,
  gasQuasiStatic = true,
  fluidQuasiStatic = true,
  metalTube(WallRes=false)) annotation(
    Placement(transformation(extent = {{-20, -20}, {20, 20}}, rotation = 0)));

  // variable for validation
  Modelica.SIunits.Power Q_out = (HE.gasIn.h_outflow - HE.gasOut.h_outflow) * HE.gasIn.m_flow; 
  Modelica.SIunits.Power Q_in = (HE.waterOut.h_outflow - HE.waterIn.h_outflow) * HE.waterIn.m_flow;
  Boolean isQMatch = abs(Q_out -Q_in) < 1e-3;  
equation

  // HE alone
  connect(source_cold.flange, HE.waterIn) annotation(
    Line(points = {{-1.83697e-015, 50}, {-1.83697e-015, 20}, {0, 20}}, color = {0, 0, 255}, thickness = 0.5, smooth = Smooth.None));
  connect(HE.waterOut, sink_cold.flange) annotation(
    Line(points = {{1.83697e-015, -70}, {1.83697e-015, -56}, {-8.88178e-016, -56}}, thickness = 0.5, color = {0, 0, 255}));
  
  connect(source_hot.flange, HE.gasIn) annotation(
        Line(points = {{-50, 0}, {-20, 0}}, color = {159, 159, 223}, thickness = 0.5, smooth = Smooth.None)); 
  connect(HE.gasOut, sink_hot.flange) annotation(
    Line(points = {{46, 0}, {46, 0}, {60, 0}}, color = {159, 159, 223}, thickness = 0.5));    
   
  /*
  connect(source_cold.flange, T_waterIn.inlet);
  connect(T_waterIn.outlet, HE.waterIn) annotation(
    Line(points = {{-1.83697e-015, 50}, {-1.83697e-015, 20}, {0, 20}}, color = {0, 0, 255}, thickness = 0.5, smooth = Smooth.None));
  connect(HE.waterOut, T_waterOut.inlet) annotation(
    Line(points = {{8.88178e-016, -44}, {8.88178e-016, -20}, {0, -20}}, thickness = 0.5, color = {0, 0, 255}));      
  connect(T_waterOut.outlet, sink_cold.flange) annotation(
    Line(points = {{1.83697e-015, -70}, {1.83697e-015, -56}, {-8.88178e-016, -56}}, thickness = 0.5, color = {0, 0, 255}));
  
  connect(source_hot.flange, T_gasIn.inlet);
  connect(T_gasIn.outlet, HE.gasIn) annotation(
        Line(points = {{-50, 0}, {-20, 0}}, color = {159, 159, 223}, thickness = 0.5, smooth = Smooth.None)); 
  connect(HE.gasOut, T_gasOut.inlet) annotation(
    Line(points = {{34, 0}, {34, 0}, {20, 0}}, color = {159, 159, 223}, thickness = 0.5));
  connect(T_gasOut.outlet, sink_hot.flange) annotation(
    Line(points = {{46, 0}, {46, 0}, {60, 0}}, color = {159, 159, 223}, thickness = 0.5));    
  */
annotation(
    Diagram(graphics),
    experiment(StartTime = 0, StopTime = 10, Tolerance = 1e-3, Interval = 2),
    __OpenModelica_commandLineOptions = "--matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection -d=initialization,NLSanalyticJacobian,aliasConflicts,bltdump",    
    __OpenModelica_simulationFlags(lv = "LOG_DEBUG,LOG_NLS,LOG_NLS_V,LOG_STATS,LOG_INIT,LOG_STDOUT, -w", outputFormat = "mat", s = "dassl", nls = "homotopy"));
end TestTP_HEG2G;
