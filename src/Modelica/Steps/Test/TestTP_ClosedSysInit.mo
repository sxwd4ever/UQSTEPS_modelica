within Steps.Test;

model TestTP_ClosedSysInit "Test for ClosedSystemInit in ThermoPower"
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
  parameter Model.PBConfiguration cfg_default;
  // select the configuration of parameters
  parameter Model.PBConfiguration cfg = cfg_default;
  // set the values of parameters accordingly
  parameter HEBoundaryCondition bc_HTR = cfg.bc_HTR;
  parameter HEBoundaryCondition bc_heater = cfg.bc_heater;
  // package Medium = Media.CO2;
  package Medium = Steps.Media.SCO2;
  //ExternalMedia.Examples.CO2CoolProp;
  ThermoPower.Gas.Turbine Turbine1(redeclare package Medium = Medium, fileName = Modelica.Utilities.Files.loadResource("modelica://Steps/Resources/Data/turbine_map.txt"), tablePhic = fill(0.0, 14, 12), tableEta = fill(0.0, 14, 12), pstart_in = bc_heater.st_cold_out.p, pstart_out = bc_HTR.st_hot_in.p, Tstart_in = bc_heater.st_cold_out.T, Tstart_out = bc_HTR.st_hot_in.T, Ndesign = 60000.0, Tdes_in = bc_heater.st_cold_out.T, Table = ThermoPower.Choices.TurboMachinery.TableTypes.file, explicitIsentropicEnthalpy = false, gas_in(p(nominal = Turbine1.pstart_in), T(nominal = Turbine1.Tstart_in)), gas_iso(p(nominal = Turbine1.pstart_out), T(nominal = Turbine1.Tstart_out))) annotation(
    Placement(visible = true, transformation(extent = {{-62, -24}, {-22, 16}}, rotation = 0)));
  //tablePhic,
  //tableEta,
  /*
              ThermoPower.Gas.SourceMassFlow SourceP1(
                redeclare package Medium = Medium, 
                T = from_degC(730.43), 
                p0 = bc_heater.st_cold_out.p,
                //h = bc_heater.st_cold_out.h, 
                use_in_T = false, 
                w0 = bc_heater.st_cold_out.mdot,    
                gas(p(nominal = bc_heater.st_cold_out.p), 
                T(nominal=bc_heater.st_cold_out.T))) 
              annotation(
                Placement(transformation(origin = {0, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 270)));
    
              ThermoPower.Gas.SinkPressure SinkP1(
              redeclare package Medium = Medium, 
              p0 =  bc_HTR.st_hot_in.p, 
              T =  bc_HTR.st_hot_in.T,
              gas(
                p(nominal = Turbine1.pstart_out), 
                T(nominal = Turbine1.Tstart_out)))
              // h = bc_HTR.st_hot_in.h) 
              annotation(
                Placement(visible = true, transformation(extent = {{50, 6}, {70, 26}}, rotation = 0)));
              */
  //ThermoPower.Gas.Utility.ClosedSystemInit sys_init(redeclare package Medium = Medium, pstart = bc_heater.st_cold_out.p);
  //w_b(nominal = bc_heater.st_cold_out.mdot),
  Modelica.Mechanics.Rotational.Sources.Speed speed1 annotation(
    Placement(visible = true, transformation(origin = {38, 0}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.Constant const1(k = 60000) annotation(
    Placement(visible = true, transformation(origin = {84, 0}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
  inner ThermoPower.System system(allowFlowReversal = false, initOpt = ThermoPower.Choices.Init.Options.noInit) annotation(
    Placement(transformation(extent = {{80, 80}, {100, 100}})));
  //ThermoPower.Gas.SensT T_gasIn(redeclare package Medium = Medium);
  //ThermoPower.Gas.SensT T_gasOut(redeclare package Medium = Medium);
  //, w0 = bc_heater.st_cold_out.mdot)
  /*
      ThermoPower.Gas.Utility.ClosedSystemInit closedSystemInit(redeclare package Medium = Medium) annotation(
        Placement(visible = true, transformation(origin = {-92, 52}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      */
  ThermoPower.Gas.SourceMassFlow sourceMassFlow(redeclare package Medium = Medium, use_in_T = false, use_in_w0 = false) annotation(
    Placement(visible = true, transformation(origin = {-78, 40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  ThermoPower.Gas.SinkMassFlow sinkMassFlow annotation(
    Placement(visible = true, transformation(origin = {-6, 42}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
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
//connect(sys_init.flange, Turbine1.outlet) annotation(    Line);
  connect(speed1.w_ref, const1.y) annotation(
    Line(points = {{50, 0}, {73, 0}}, color = {0, 0, 127}));
  connect(Turbine1.shaft_b, speed1.flange) annotation(
    Line(points = {{-30, -4}, {9, -4}, {9, 0}, {28, 0}}));
/*connect(closedSystemInit.flange, throughMassFlow.inlet) annotation(
    Line(points = {{-92, 44}, {-34, 44}, {-34, 46}, {-34, 46}}, color = {159, 159, 223}));
  */
  connect(sourceMassFlow.flange, Turbine1.inlet) annotation(
    Line(points = {{-68, 40}, {-58, 40}, {-58, 12}, {-58, 12}}, color = {159, 159, 223}));
  connect(Turbine1.outlet, sinkMassFlow.flange) annotation(
    Line(points = {{-26, 12}, {-16, 12}, {-16, 42}, {-16, 42}}, color = {159, 159, 223}));
  annotation(
    Diagram(graphics),
    experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-3, Interval = 1),
    __OpenModelica_commandLineOptions = "--matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection -d=initialization,NLSanalyticJacobian,aliasConflicts,bltdump",
    __OpenModelica_simulationFlags(lv = "LOG_DEBUG,LOG_NLS,LOG_NLS_V,LOG_STATS,LOG_INIT,LOG_STDOUT, -w", outputFormat = "mat", s = "dassl", nls = "homotopy"));
end TestTP_ClosedSysInit;
