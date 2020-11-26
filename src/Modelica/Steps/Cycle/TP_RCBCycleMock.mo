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
  parameter Modelica.SIunits.Temperature T_heater_hot = from_degC(700) "K, Temperature of heater's hot fluid";  
  parameter Modelica.SIunits.Temperature T_cooler_cold = from_degC(35) "K, Temperature of cooler's cold fluid";  
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

  Components.GasStateReader r01(redeclare package Medium = medium_main) annotation(
    Placement(visible = true, transformation(origin = {66, -56}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
  Components.GasStateReader r02(redeclare package Medium = medium_main) annotation(
    Placement(visible = true, transformation(origin = {96, 44}, extent = {{-6, -6}, {6, 6}}, rotation = 180)));
  Components.GasStateReader r03(redeclare package Medium = medium_main) annotation(
    Placement(visible = true, transformation(origin = {-22, 50}, extent = {{-6, -6}, {6, 6}}, rotation = -90)));
  Components.GasStateReader r04(redeclare package Medium = medium_main) annotation(
    Placement(visible = true, transformation(origin = {-32, 20}, extent = {{-6, -6}, {6, 6}}, rotation = 180)));
  Components.GasStateReader r05(redeclare package Medium = medium_main) annotation(
    Placement(visible = true, transformation(origin = {-92, 18}, extent = {{-6, -6}, {6, 6}}, rotation = -90)));
  Components.GasStateReader r05a(redeclare package Medium = medium_main) annotation(
    Placement(visible = true, transformation(origin = {-92, -32}, extent = {{-6, -6}, {6, 6}}, rotation = -90)));
  Components.GasStateReader r06(redeclare package Medium = medium_main) annotation(
    Placement(visible = true, transformation(origin = {-46, -16}, extent = {{-6, -6}, {6, 6}}, rotation = 90)));
  Components.GasStateReader r07(redeclare package Medium = medium_main) annotation(
    Placement(visible = true, transformation(origin = {-4, 34}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
  Components.GasStateReader r08(redeclare package Medium = medium_main) annotation(
    Placement(visible = true, transformation(origin = {32, 18}, extent = {{-6, -6}, {6, 6}}, rotation = -90)));
  Components.GasStateReader r08a(redeclare package Medium = medium_main) annotation(
    Placement(visible = true, transformation(origin = {20, -24}, extent = {{-6, -6}, {6, 6}}, rotation = -90)));
  Components.GasStateReader r08b(redeclare package Medium = medium_main) annotation(
    Placement(visible = true, transformation(origin = {34, 4}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
  Components.GasStateReader r09(redeclare package Medium = medium_main) annotation(
    Placement(visible = true, transformation(origin = {50, 62}, extent = {{-6, -6}, {6, 6}}, rotation = 180)));
  Components.GasStateReader r10(redeclare package Medium = medium_main) annotation(
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


  Components.WaterStateReader rh1(redeclare package Medium = medium_water) annotation(
    Placement(visible = true, transformation(origin = {-66, 22}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
  Components.WaterStateReader rh2(redeclare package Medium = medium_water) annotation(
    Placement(visible = true, transformation(origin = {-64, 64}, extent = {{-6, -6}, {6, 6}}, rotation = 180)));
  Components.WaterStateReader rc1(redeclare package Medium = medium_water) annotation(
    Placement(visible = true, transformation(origin = {44, -36}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
  Components.WaterStateReader rc2(redeclare package Medium = medium_water) annotation(
    Placement(visible = true, transformation(origin = {44, -74}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));

  // (Mock) calculated variables
  // for turbine
  Modelica.SIunits.Power W_turb = (r05a.h - r06.h) * r05a.w / 1e6 "W->MW, net power for turbine";
  Modelica.SIunits.Efficiency eta_turb = 0.85;
  // main compressor
  Modelica.SIunits.Power W_MC = (r02.h - r01.h) * r01.w / 1e6 "W->MW, net power for main compressor";
  Modelica.SIunits.Efficiency eta_MC = 0.86;
  // re compressor
  Modelica.SIunits.Power W_RC = (r09.h - r08b.h) * r08b.w / 1e6 "W->MW, net power for recompressor";
  Modelica.SIunits.Efficiency eta_RC = 0.87;
  
  // power block
  Modelica.SIunits.Power W_net = W_turb - W_MC - W_RC "net power generated";
  Modelica.SIunits.Power Q_heater = (rh1.h - rh2.h) * rh1.w / 1e6 + 1 "heat input for the power block";
  Modelica.SIunits.Efficiency eta_pb = W_net / Q_heater * 100 "power block efficiency";
  Real SR = r08b.w / r08.w * 100 "split ratio";
  
  // heat transfer coefficient for HTR and LTR
  // HTR
  Modelica.SIunits.Power Q_HTR = (r04.h - r03.h) * r03.w / 1e6 "W->MW, heat input for HTR";
  Modelica.SIunits.TemperatureDifference dT1_HTR = (r06.T - r04.T) + 1;
  Modelica.SIunits.TemperatureDifference dT2_HTR = (r07.T - r03.T);
  Real T_ltmd_HTR = (dT2_HTR - dT1_HTR) / Modelica.Math.log(dT2_HTR / dT1_HTR + 0.5);
  Real UA_HTR = Q_HTR / T_ltmd_HTR;
  // LTR
  Modelica.SIunits.Power Q_LTR = (r10.h - r02.h) * r02.w / 1e6 "W->MW";
  Modelica.SIunits.TemperatureDifference dT1_LTR = (r07.T - r10.T) + 1;
  Modelica.SIunits.TemperatureDifference dT2_LTR = (r08.T - r02.T);
  Real T_ltmd_LTR = (dT2_LTR - dT1_LTR) / Modelica.Math.log(dT2_LTR / dT1_LTR + 0.5);
  Real UA_LTR = Q_LTR / T_ltmd_LTR;  
  
  // Liquid Na exit temperature
  Modelica.SIunits.Temperature T_heater_hot_out = rh2.T;
  
  
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
