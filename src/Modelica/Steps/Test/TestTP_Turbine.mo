within Steps.Test;

model TestTP_Turbine
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
/*  
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
*/
//extends Modelica.Icons.Example;
  //package Medium = Modelica.Media.IdealGases.MixtureGases.CombustionAir;
  package Medium = Modelica.Media.IdealGases.SingleGases.CO2;
  //package Medium = Media.CO2;//Modelica.Media.IdealGases.SingleGases.CO2;
  ThermoPower.Gas.SourcePressure SourceP1(redeclare package Medium = Medium, T = 700 + 273.15, p0 = 20.0e6) annotation(
    Placement(transformation(extent = {{-80, 6}, {-60, 26}}, rotation = 0)));
  Modelica.Mechanics.Rotational.Components.Inertia Inertia1(J = 10000) annotation(
    Placement(transformation(extent = {{10, -10}, {30, 10}}, rotation = 0)));
  ThermoPower.Gas.Turbine Turbine1(redeclare package Medium = Medium, tablePhic = tablePhic, tableEta = tableEta, pstart_in = 20.0e6, pstart_out = 8.0e6, Tstart_in = 700 + 273.15, Tstart_out = 616 + 273.15, Ndesign = 60000.0, Tdes_in = 700 + 273.15, Table = ThermoPower.Choices.TurboMachinery.TableTypes.matrix) annotation(
    Placement(transformation(extent = {{-40, -20}, {0, 20}}, rotation = 0)));
  ThermoPower.Gas.SinkPressure SinkP1(redeclare package Medium = Medium, p0 = 8.0e6, T = 616 + 273.15) annotation(
    Placement(visible = true, transformation(extent = {{50, 6}, {70, 26}}, rotation = 0)));
  inner ThermoPower.System system annotation(
    Placement(transformation(extent = {{80, 80}, {100, 100}})));
  ThermoPower.Gas.SensT sensT1(redeclare package Medium = Medium) annotation(
    Placement(visible = true, transformation(origin = {20, 44}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
protected
  parameter Real tablePhic[5, 4] = [1, 37, 80, 100; 1.5, 7.10E-05, 7.10E-05, 7.10E-05; 2, 8.40E-05, 8.40E-05, 8.40E-05; 2.5, 8.70E-05, 8.70E-05, 8.70E-05; 3, 1.04E-04, 1.04E-04, 1.04E-04];
  parameter Real tableEta[5, 4] = [1, 37, 80, 100; 1.5, 0.57, 0.89, 0.81; 2, 0.46, 0.82, 0.88; 2.5, 0.41, 0.76, 0.85; 3, 0.38, 0.72, 0.82];
equation
  connect(SourceP1.flange, Turbine1.inlet) annotation(
    Line(points = {{-60, 16}, {-36, 16}}, color = {159, 159, 223}, thickness = 0.5));
  connect(Turbine1.outlet, sensT1.inlet) annotation(
    Line(points = {{-4, 16}, {6, 16}, {6, 40}, {14, 40}, {14, 40}}, color = {159, 159, 223}));
  connect(sensT1.outlet, SinkP1.flange) annotation(
    Line(points = {{26, 40}, {36, 40}, {36, 16}, {50, 16}}, color = {159, 159, 223}));
initial equation
  Inertia1.w = 60000.0;
equation
  connect(Turbine1.shaft_b, Inertia1.flange_a) annotation(
    Line(points = {{-8, 0}, {-4, 0}, {-4, 0}, {10, 0}}, color = {0, 0, 0}, thickness = 0.5));
  annotation(
    Diagram(graphics),
    experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-3, Interval = 1),
    __OpenModelica_commandLineOptions = "--matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection -d=initialization,NLSanalyticJacobian,aliasConflicts,bltdump",    
    __OpenModelica_simulationFlags(lv = "LOG_DEBUG,LOG_NLS,LOG_NLS_V,LOG_STATS,LOG_INIT,LOG_STDOUT, -w", outputFormat = "mat", s = "dassl", nls = "homotopy"));

end TestTP_Turbine;
