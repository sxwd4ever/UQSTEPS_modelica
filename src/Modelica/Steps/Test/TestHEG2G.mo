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
  package FlueGasMedium = Steps.Media.SCO2;
  package FluidMedium = Steps.Media.SCO2;
    
  // parameter for C++ implementation of PCHE - based on Modelica impl's result
  parameter Model.PBConfiguration cfg_def( 
  p_pump_in = 9e6,
  p_pump_out = 20e6,  
  mdot_main = 125,
  mdot_pump = cfg_def.mdot_main * 0.7,    
  bc_HTR.st_hot_in.T = 883,
  bc_HTR.st_hot_out.T = 643,
  bc_HTR.st_cold_in.T = 633,
  bc_HTR.st_cold_out.T = 843);  
  
  // select the configuration of parameters
  parameter Model.PBConfiguration cfg = cfg_def;
  
  // set the values of parameters accordingly
  parameter PCHEBoundaryCondition bc = cfg.bc_HTR;  
  
  parameter EntityGeoParam geo_hot = cfg.cfg_hot.geo;
  parameter EntityGeoParam geo_cold = cfg.cfg_cold.geo;
  parameter EntityGeoParam geo_tube = cfg.cfg_tube.geo;
  
  parameter EntityThermoParam thermo_hot = cfg.cfg_hot.thermo;
  parameter EntityThermoParam thermo_cold = cfg.cfg_cold.thermo;
  parameter EntityThermoParam thermo_tube = cfg.cfg_tube.thermo;  
  
  parameter Modelica.SIunits.Temperature Tstart_M_In = bc.st_cold_in.T "Inlet metal wall temperature start value";
  parameter Modelica.SIunits.Temperature Tstart_M_Out = bc.st_cold_out.T "Outlet metal wall temperature start value";

  //parameter Modelica.SIunits.SpecificEnthalpy hstart_F_In = FluidMedium.specificEnthalpy_pT(fluidNomPressure, bc.st_cold_in.T) "Nominal specific enthalpy";
  //parameter Modelica.SIunits.SpecificEnthalpy hstart_F_Out = FluidMedium.specificEnthalpy_pT(fluidNomPressure, bc.st_cold_out.T) "Nominal specific enthalpy";
  //Components
  inner ThermoPower.System system(allowFlowReversal = false) annotation(
    Placement(transformation(extent = {{80, 80}, {100, 100}})));
  ThermoPower.Gas.SourceMassFlow sourceW_water(redeclare package Medium = FluidMedium, T = bc.st_cold_in.T, p0 = bc.st_cold_in.p, use_in_T = false, w0 = bc.st_cold_in.mdot) annotation(
    Placement(transformation(origin = {0, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 270)));
  ThermoPower.Gas.SinkPressure sinkP_water(redeclare package Medium = FluidMedium, p0 = bc.st_cold_in.p, T = bc.st_cold_out.T) annotation(
    Placement(transformation(origin = {0, -80}, extent = {{-10, -10}, {10, 10}}, rotation = 270)));
  ThermoPower.Gas.SinkPressure sinkP_gas(redeclare package Medium = FlueGasMedium, T = bc.st_hot_out.T, p0 = bc.st_hot_in.p) annotation(
    Placement(transformation(extent = {{60, -10}, {80, 10}}, rotation = 0)));
  ThermoPower.Gas.SourceMassFlow sourceW_gas(redeclare package Medium = FlueGasMedium, T = bc.st_hot_in.T, p0 = bc.st_hot_in.p, w0 = bc.st_hot_in.mdot) annotation(
    Placement(transformation(extent = {{-70, -10}, {-50, 10}}, rotation = 0)));
  ThermoPower.Gas.SensT T_waterOut(redeclare package Medium = FluidMedium) annotation(
    Placement(transformation(origin = {4, -50}, extent = {{-10, -10}, {10, 10}}, rotation = 270)));
  ThermoPower.Gas.SensT T_gasOut(redeclare package Medium = FlueGasMedium) annotation(
    Placement(transformation(extent = {{30, -6}, {50, 14}}, rotation = 0)));
  Components.HEG2G hE(
  redeclare package FluidMedium = FluidMedium, 
  redeclare package FlueGasMedium = FlueGasMedium, 
  redeclare replaceable model HeatTransfer_F = ThermoPower.Thermal.HeatTransferFV.ConstantHeatTransferCoefficient(gamma = thermo_cold.gamma_he), 
  redeclare replaceable model HeatTransfer_G = ThermoPower.Thermal.HeatTransferFV.ConstantHeatTransferCoefficientTwoGrids(gamma = thermo_hot.gamma_he), 
  redeclare model HeatExchangerTopology = ThermoPower.Thermal.HeatExchangerTopologies.CounterCurrentFlow,  
  N_F = geo_cold.N_seg, 
  N_G = geo_hot.N_seg,   
  SSInit = SSInit, 
  Tstartbar_G = bc.st_hot_in.T, 
  exchSurface_F = geo_cold.A_ex, 
  exchSurface_G = geo_hot.A_ex, 
  extSurfaceTub = geo_tube.A_ex, 
  fluidNomFlowRate = bc.st_cold_in.mdot, 
  fluidNomPressure = bc.st_cold_in.p, 
  fluidVol = geo_cold.V, 
  gasNomFlowRate = bc.st_hot_in.mdot, 
  gasNomPressure = bc.st_hot_in.p, 
  gasVol = geo_hot.V, 
  lambda = thermo_tube.lambda, 
  metalVol = geo_tube.V, 
  pstart_F = bc.st_cold_in.p, 
  rhomcm = thermo_tube.rho_mcm) annotation(
    Placement(transformation(extent = {{-20, -20}, {20, 20}}, rotation = 0)));
  //Start value
  // parameter Modelica.SIunits.Temperature Tstart_G = (bc.st_hot_in.T + bc.st_hot_out.T) / 2;
  // parameter Modelica.SIunits.Temperature Tstart_M = (bc.st_hot_in.T + bc.st_hot_out.T + bc.st_cold_in.T + bc.st_cold_out.T) / 4;
  parameter Boolean SSInit = true "Steady-state initialization";
  
initial equation
//hstart_F_Out = hE.waterOut.h_outflow;
equation
  connect(T_gasOut.inlet, hE.gasOut) annotation(
    Line(points = {{34, 0}, {34, 0}, {20, 0}}, color = {159, 159, 223}, thickness = 0.5));
  connect(T_gasOut.outlet, sinkP_gas.flange) annotation(
    Line(points = {{46, 0}, {46, 0}, {60, 0}}, color = {159, 159, 223}, thickness = 0.5));
  connect(sinkP_water.flange, T_waterOut.outlet) annotation(
    Line(points = {{1.83697e-015, -70}, {1.83697e-015, -56}, {-8.88178e-016, -56}}, thickness = 0.5, color = {0, 0, 255}));
  connect(T_waterOut.inlet, hE.waterOut) annotation(
    Line(points = {{8.88178e-016, -44}, {8.88178e-016, -20}, {0, -20}}, thickness = 0.5, color = {0, 0, 255}));
  connect(sourceW_water.flange, hE.waterIn) annotation(
    Line(points = {{-1.83697e-015, 50}, {-1.83697e-015, 20}, {0, 20}}, color = {0, 0, 255}, thickness = 0.5, smooth = Smooth.None));
  connect(sourceW_gas.flange, hE.gasIn) annotation(
        Line(points = {{-50, 0}, {-20, 0}}, color = {159, 159, 223}, thickness = 0.5, smooth = Smooth.None));

annotation(
    Diagram(graphics),
    experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-3, Interval = 1),
    __OpenModelica_commandLineOptions = "--matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection -d=initialization,NLSanalyticJacobian,aliasConflicts,bltdump");  
end TestHEG2G;
