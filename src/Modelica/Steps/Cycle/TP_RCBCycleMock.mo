within Steps.Cycle;

model TP_RCBCycleMock
  " to accelerate the development and integration with Python script"
  
import Modelica.SIunits.Conversions.{from_degC,from_deg};
  import Modelica.SIunits.{Temperature,Pressure,SpecificEnthalpy};
  import Util = Utilities.Util;
  import Steps.Utilities.CoolProp.PropsSI;
  import Steps.Components.{PCHEGeoParam};
  import Steps.Model.{PBConfiguration,SimParam,EntityConfig,EntityGeoParam,EntityThermoParam,ThermoState,HEBoundaryCondition};
  import Model.PBConfiguration;
  import ThermoPower.Choices.Init.Options;
  import ThermoPower.System;
  import ThermoPower.Gas;
  
  // input parameters of the power block
  parameter Modelica.SIunits.MassFlowRate mdot_main = 125 "kg/s, mass flow in the main path of PB, which follows the power demand";
  parameter Modelica.SIunits.MassFlowRate mdot_heater_hot = 55 "kg/s, mass flow rate of heater's hot fluid";
  parameter Modelica.SIunits.Temperature T_heater_hot = from_degC(800) "K, Temperature of heater's hot fluid";  
  parameter Modelica.SIunits.Temperature T_cooler_cold = from_degC(45) "K, Temperature of cooler's cold fluid";  
  parameter Modelica.SIunits.AbsolutePressure p_main = from_bar(80);
  parameter Modelica.SIunits.AbsolutePressure p_water = from_bar(10);

  replaceable package medium_main = Steps.Media.SCO2;
  replaceable package medium_water = Modelica.Media.Water.StandardWater;
  
  inner ThermoPower.System system(allowFlowReversal = false, initOpt=ThermoPower.Choices.Init.Options.noInit) annotation(
    Placement(transformation(extent = {{100, 80}, {120, 100}})));
    
  ThermoPower.Gas.SourceMassFlow source(
  redeclare package Medium = medium_main,
    w0 = mdot_main,
    p0 = p_main,
    T = T_heater_hot,
    gas(T(nominal = source.T), p(nominal=source.p0))) annotation(
    Placement(visible = true, transformation(origin = {-91, -17}, extent = {{-5, -5}, {5, 5}}, rotation = -90)));
    
  ThermoPower.Gas.SinkPressure sink(
  redeclare package Medium = medium_main,
    p0 = p_main,
    T = T_heater_hot) annotation(
    Placement(visible = true, transformation(origin = {-92, -4}, extent = {{-4, -4}, {4, 4}}, rotation = 270)));  

  ThermoPower.PowerPlants.HRSG.Components.StateReader_gas r01(redeclare package Medium = medium_main) annotation(
    Placement(visible = true, transformation(origin = {66, -56}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
  ThermoPower.PowerPlants.HRSG.Components.StateReader_gas r02(redeclare package Medium = medium_main) annotation(
    Placement(visible = true, transformation(origin = {96, 44}, extent = {{-6, -6}, {6, 6}}, rotation = 180)));
  ThermoPower.PowerPlants.HRSG.Components.StateReader_gas r03(redeclare package Medium = medium_main) annotation(
    Placement(visible = true, transformation(origin = {-22, 50}, extent = {{-6, -6}, {6, 6}}, rotation = -90)));
  ThermoPower.PowerPlants.HRSG.Components.StateReader_gas r04(redeclare package Medium = medium_main) annotation(
    Placement(visible = true, transformation(origin = {-32, 20}, extent = {{-6, -6}, {6, 6}}, rotation = 180)));
  ThermoPower.PowerPlants.HRSG.Components.StateReader_gas r05(redeclare package Medium = medium_main) annotation(
    Placement(visible = true, transformation(origin = {-92, 18}, extent = {{-6, -6}, {6, 6}}, rotation = -90)));
  ThermoPower.PowerPlants.HRSG.Components.StateReader_gas r05a(redeclare package Medium = medium_main) annotation(
    Placement(visible = true, transformation(origin = {-92, -32}, extent = {{-6, -6}, {6, 6}}, rotation = -90)));
  ThermoPower.PowerPlants.HRSG.Components.StateReader_gas r06(redeclare package Medium = medium_main) annotation(
    Placement(visible = true, transformation(origin = {-46, -16}, extent = {{-6, -6}, {6, 6}}, rotation = 90)));
  ThermoPower.PowerPlants.HRSG.Components.StateReader_gas r07(redeclare package Medium = medium_main) annotation(
    Placement(visible = true, transformation(origin = {-4, 34}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
  ThermoPower.PowerPlants.HRSG.Components.StateReader_gas r08(redeclare package Medium = medium_main) annotation(
    Placement(visible = true, transformation(origin = {32, 18}, extent = {{-6, -6}, {6, 6}}, rotation = -90)));
  ThermoPower.PowerPlants.HRSG.Components.StateReader_gas r08a(redeclare package Medium = medium_main) annotation(
    Placement(visible = true, transformation(origin = {20, -24}, extent = {{-6, -6}, {6, 6}}, rotation = -90)));
  ThermoPower.PowerPlants.HRSG.Components.StateReader_gas r08b(redeclare package Medium = medium_main) annotation(
    Placement(visible = true, transformation(origin = {34, 4}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
  ThermoPower.PowerPlants.HRSG.Components.StateReader_gas r09(redeclare package Medium = medium_main) annotation(
    Placement(visible = true, transformation(origin = {50, 62}, extent = {{-6, -6}, {6, 6}}, rotation = 180)));
  ThermoPower.PowerPlants.HRSG.Components.StateReader_gas r10(redeclare package Medium = medium_main) annotation(
    Placement(visible = true, transformation(origin = {4, 58}, extent = {{-6, -6}, {6, 6}}, rotation = 180)));

  ThermoPower.Water.SourceMassFlow source_water(
  redeclare package Medium = medium_main,
    w0 = mdot_heater_hot,
    p0 = p_water,
    T = T_cooler_cold);
    
  ThermoPower.Water.SinkPressure sink_water(
  redeclare package Medium = medium_main,
    p0 = p_water,
    T = T_cooler_cold);


  ThermoPower.PowerPlants.HRSG.Components.StateReader_water rh1(redeclare package Medium = medium_water) annotation(
    Placement(visible = true, transformation(origin = {-66, 22}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
  ThermoPower.PowerPlants.HRSG.Components.StateReader_water rh2(redeclare package Medium = medium_water) annotation(
    Placement(visible = true, transformation(origin = {-64, 64}, extent = {{-6, -6}, {6, 6}}, rotation = 180)));
  ThermoPower.PowerPlants.HRSG.Components.StateReader_water rc1(redeclare package Medium = medium_water) annotation(
    Placement(visible = true, transformation(origin = {44, -36}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
  ThermoPower.PowerPlants.HRSG.Components.StateReader_water rc2(redeclare package Medium = medium_water) annotation(
    Placement(visible = true, transformation(origin = {44, -74}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));


equation
  // open loop with compressor + recompressor
  // main path - source -> sink
  connect(source.flange, r01.inlet); 
  connect(r01.outlet, r02.inlet);
  connect(r02.outlet, r03.inlet);
  connect(r03.outlet, r04.inlet);
  connect(r04.outlet, r05.inlet);
  connect(r05.outlet, r05a.inlet);
  connect(r05a.outlet, r06.inlet);
  connect(r06.outlet, r07.inlet);
  connect(r07.outlet, r08.inlet);
  connect(r08.outlet, r08a.inlet);
  connect(r08a.outlet, r08b.inlet);
  connect(r08b.outlet, r09.inlet);
  connect(r09.outlet, r10.inlet);
  connect(r10.outlet, sink.flange);

  connect(source_water.flange, rh1.inlet); 
  connect(rh1.outlet, rh2.inlet);
  connect(rh2.outlet, rc1.inlet);
  connect(rc1.outlet, rc2.inlet);  
  connect(rc2.outlet, sink_water.flange);

  annotation(
    Diagram(coordinateSystem(extent = {{-100, -100}, {120, 100}})),
    experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-3, Interval = 1),
    __OpenModelica_commandLineOptions = "--matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection -d=initialization,NLSanalyticJacobian,bltdump",
    __OpenModelica_simulationFlags(lv = "LOG_DEBUG,LOG_NLS,LOG_NLS_V,LOG_STATS,LOG_INIT,LOG_STDOUT, -w", outputFormat = "mat", s = "dassl", nls = "homotopy"));
end TP_RCBCycleMock;
