within Steps.Test.TPComponents;

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
*/
  
  // parameter for C++ implementation of PCHE - based on Modelica impl's result    
  parameter Model.RCBCycleConfig cfg(
    redeclare package medium_main = Medium
  );
  /*
  (
    mdot_main = 100,
    T_heater_cold_out = from_degC(700),
    Ns_comp = 30000
  );
  */
  // set the values of parameters accordingly
  parameter Model.TurbomachineryConfig cfg_turb = cfg.cfg_turb;
  parameter Model.ThermoState st_source         = cfg_turb.st_in;
  parameter Model.ThermoState st_sink           = cfg_turb.st_out;
  
  // package Medium = Media.CO2;
  package Medium = Steps.Media.SCO2(
    inputChoice = ExternalMedia.Common.InputChoice.pT,
    substanceNames = {"CO2|debug=40"}    
  );//ExternalMedia.Examples.CO2CoolProp;
  
  // package Medium = ExternalMedia.Examples.CO2CoolProp;

  ThermoPower.Gas.SourceMassFlow SourceP1(
    redeclare package Medium = Medium, 
    T        = st_source.T,
    p0       = st_source.p,
    use_in_T = false,
    w0       = st_source.mdot,
    gas(
      p(nominal = st_source.p), 
      T(nominal = st_source.T),
      h(nominal = st_source.h))) 
  annotation(
    Placement(transformation(origin = {0, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 270)));

 /* 
  ThermoPower.Gas.SourcePressure SourceP1(
    redeclare package Medium = Medium, 
    T = bc_heater.st_cold_out.T, 
    p0 = bc_heater.st_cold_out.p,
    //h = bc_heater.st_cold_out.h, 
    use_in_T = false, 
    //w0 = bc_heater.st_cold_out.mdot,    
    gas(p(nominal = bc_heater.st_cold_out.p), 
    T(nominal=bc_heater.st_cold_out.T))) 
  annotation(
    Placement(transformation(origin = {0, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 270)));
  */ 
  
  ThermoPower.Gas.Turbine Turbine1(
    redeclare package Medium = Medium, 
    fileName                   = Modelica.Utilities.Files.loadResource("modelica://Steps/Resources/Data/turbine_map_10MW.txt"),
    tablePhic                  = fill(0.0, 14, 12),                                                                               //tablePhic, 
    tableEta                   = fill(0.0, 14, 12),                                                                               //tableEta, 
    pstart_in                  = cfg_turb.st_in.p,
    pstart_out                 = cfg_turb.st_out.p,
    Tstart_in                  = cfg_turb.st_in.T,
    Tstart_out                 = cfg_turb.st_out.T,
    Ndesign                    = cfg_turb.N,
    Tdes_in                    = cfg_turb.st_in.T,
    Table                      = ThermoPower.Choices.TurboMachinery.TableTypes.file,
    //explicitIsentropicEnthalpy = false,
    gas_in(
      p(nominal = Turbine1.pstart_in), 
      T(nominal = Turbine1.Tstart_in)),
    gas_iso(
      p(nominal = Turbine1.pstart_out), 
      T(nominal = Turbine1.Tstart_out)))    
    annotation(
      Placement(transformation(extent = {{-40, -20}, {0, 20}}, rotation = 0)));
  
  ThermoPower.Gas.SinkPressure SinkP1(
  redeclare package Medium = Medium, 
  p0 = st_sink.p,
  T  = st_sink.T,
  gas(
    p(nominal = st_sink.p), 
    T(nominal = st_sink.T)))
  // h = bc_HTR.st_hot_in.h) 
  annotation(
    Placement(visible = true, transformation(extent = {{50, 6}, {70, 26}}, rotation = 0)));

  Modelica.Mechanics.Rotational.Sources.Speed speed1 annotation(
    Placement(visible = true, transformation(origin = {84, 0}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.Constant const1(k = cfg_turb.N) annotation(
    Placement(visible = true, transformation(origin = {130, 0}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));  

  inner ThermoPower.System system(allowFlowReversal = false, initOpt=ThermoPower.Choices.Init.Options.noInit) annotation(
    Placement(transformation(extent = {{80, 80}, {100, 100}})));

  //ThermoPower.Gas.SensT T_gasIn(redeclare package Medium = Medium);
  //ThermoPower.Gas.SensT T_gasOut(redeclare package Medium = Medium);
  Modelica.SIunits.Power W_turb = (Turbine1.gas_in.h - Turbine1.hout) * Turbine1.inlet.m_flow / 1e6 "W->MW, net power for turbine";
  
  Modelica.SIunits.Efficiency eta_turb = Turbine1.eta * 100;
  
protected

  parameter Real tablePhic[5, 4] = [1, 37, 80, 100; 1.5, 7.10E-05, 7.10E-05, 7.10E-05; 2, 8.40E-05, 8.40E-05, 8.40E-05; 2.5, 8.70E-05, 8.70E-05, 8.70E-05; 3, 1.04E-04, 1.04E-04, 1.04E-04];
  parameter Real tableEta[5, 4] = [1, 37, 80, 100; 1.5, 0.57, 0.89, 0.81; 2, 0.46, 0.82, 0.88; 2.5, 0.41, 0.76, 0.85; 3, 0.38, 0.72, 0.82];

equation
  /*
  connect(SourceP1.flange, T_gasIn.inlet);
  connect(T_gasIn.outlet, Turbine1.inlet) annotation(
    Line(points = {{-60, 16}, {-36, 16}}, color = {159, 159, 223}, thickness = 0.5));

  connect(Turbine1.outlet, T_gasOut.inlet) annotation(
    Line(points = {{-4, 16}, {6, 16}, {6, 40}, {14, 40}, {14, 40}}, color = {159, 159, 223}));
  connect(T_gasOut.outlet, SinkP1.flange) annotation(
    Line(points = {{26, 40}, {36, 40}, {36, 16}, {50, 16}}, color = {159, 159, 223}));
  */
  
  connect(SourceP1.flange, Turbine1.inlet) annotation(
    Line(points = {{-60, 16}, {-36, 16}}, color = {159, 159, 223}, thickness = 0.5));

  connect(Turbine1.outlet, SinkP1.flange) annotation(
    Line(points = {{26, 40}, {36, 40}, {36, 16}, {50, 16}}, color = {159, 159, 223}));

  connect(Turbine1.shaft_b, speed1.flange) annotation(
    Line(points = {{30, 0}, {74, 0}, {74, 0}, {74, 0}}));
  connect(speed1.w_ref, const1.y) annotation(
    Line(points = {{96, 0}, {120, 0}, {120, 0}, {118, 0}}, color = {0, 0, 127}));

  annotation(
    Diagram(graphics),
    experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-3, Interval = 1),
    __OpenModelica_commandLineOptions = "--matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection -d=initialization,NLSanalyticJacobian,aliasConflicts,bltdump",    
    __OpenModelica_simulationFlags(lv = "LOG_DEBUG,LOG_NLS,LOG_NLS_V,LOG_STATS,LOG_INIT,LOG_STDOUT, -w", outputFormat = "mat", s = "dassl", nls = "homotopy"));

end TestTP_Turbine;
