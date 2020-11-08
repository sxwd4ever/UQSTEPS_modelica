within Steps.Test;

model TestTP_Compressor
  "Test for HE in ThermoPower"  
  import Modelica.SIunits.Conversions.{from_degC, from_deg};
  import Modelica.SIunits.{Temperature, Pressure, SpecificEnthalpy};
  import Util = Utilities.Util;
  import Steps.Utilities.CoolProp.PropsSI;  
  import Steps.Components.{PCHEGeoParam};
  import Steps.Model.PBConfiguration;
  import Steps.Model.{PBConfiguration, SimParam, EntityConfig, EntityGeoParam, EntityThermoParam, ThermoState, HEBoundaryCondition} ;
  import ThermoPower.Choices.Init.Options;
  import ThermoPower.System;
  import ThermoPower.Gas;
  
  // package medium_hot = Steps.Media.CO2;
  // package medium_cooler = ThermoPower.Water.StandardWater;//Modelica.Media.IdealGases.SingleGases.CO2;  
  
  package Medium = Steps.Media.CO2; //Modelica.Media.IdealGases.SingleGases.CO2;
  
  parameter Model.PBConfiguration cfg_on_design( 
  p_pump_in = 8e6,
  p_pump_out = 20e6,  
  mdot_main = 51.51,
  mdot_pump = 31.31,
  mdot_heater = 20,
  T_HTR_cold_in = from_degC(141.3), 
  T_HTR_cold_out = from_degC(495.302),
  T_HTR_hot_out = from_degC(141.041),
  T_LTR_hot_out = from_degC(63.6726)); 
  
  // select the configuration of parameters
  parameter Model.PBConfiguration cfg = cfg_on_design;
  
    // set the values of parameters accordingly
  parameter HEBoundaryCondition bc = cfg.bc_LTR;  
  
  ThermoPower.Gas.SourcePressure SourceP1(
    redeclare package Medium = Medium,
    p0=bc.st_hot_out.p,
    T=bc.st_hot_out.T) annotation (Placement(transformation(extent={{-80,6},{-60,26}},
          rotation=0)));
  ThermoPower.Gas.SinkPressure SinkP1(
    redeclare package Medium = Medium,
    p0=bc.st_cold_in.p,
    T=bc.st_cold_in.T) annotation (Placement(transformation(extent={{40,6},{60,26}},
          rotation=0)));
          
  ThermoPower.Gas.Compressor Compressor(
    redeclare package Medium = Medium,
    pstart_in=bc.st_hot_out.p,
    pstart_out=bc.st_cold_in.p,
    Tstart_in=bc.st_hot_out.T,
    Tstart_out=bc.st_cold_in.T,
    tablePhic=tablePhic,
    tableEta=tableEta,
    tablePR=tablePR,
    Table=ThermoPower.Choices.TurboMachinery.TableTypes.matrix,
    Ndesign=523.3,
    Tdes_in=244.4) annotation (Placement(transformation(extent={{-20,-20},{
            20,20}}, rotation=0)));
            
  Modelica.Mechanics.Rotational.Sources.ConstantSpeed ConstantSpeed1(
      w_fixed=523.3, useSupport=false) annotation (Placement(transformation(
          extent={{-50,-10},{-30,10}}, rotation=0)));
  inner ThermoPower.System system
    annotation (Placement(transformation(extent={{80,80},{100,100}})));
protected
  parameter Real tableEta[6, 4]=[0, 95, 100, 105; 1, 82.5e-2, 81e-2,
      80.5e-2; 2, 84e-2, 82.9e-2, 82e-2; 3, 83.2e-2, 82.2e-2, 81.5e-2; 4,
      82.5e-2, 81.2e-2, 79e-2; 5, 79.5e-2, 78e-2, 76.5e-2];
  parameter Real tablePhic[6, 4]=[0, 95, 100, 105; 1, 38.3e-3/400, 43e-3/400,
      46.8e-3/400; 2, 39.3e-3/400, 43.8e-3/400, 47.9e-3/400; 3, 40.6e-3/400, 45.2e-3/400, 48.4e-3/400;
      4, 41.6e-3/400, 46.1e-3/400, 48.9e-3/400; 5, 42.3e-3/400, 46.6e-3/400, 49.3e-3/400];

  parameter Real tablePR[6, 4]=[0, 95, 100, 105; 1, 22.6/12, 27/12, 32/12; 2, 22/12,
      26.6/12, 30.8/12; 3, 20.8/12, 25.5/12, 29/12; 4, 19/12, 24.3/12, 27.1/12; 5, 17/12, 21.5/12, 24.2/12];
equation
  connect(SourceP1.flange, Compressor.inlet) annotation (Line(
      points={{-60,16},{-16,16}},
      color={159,159,223},
      thickness=0.5));
  connect(Compressor.outlet, SinkP1.flange) annotation (Line(
      points={{16,16},{40,16}},
      color={159,159,223},
      thickness=0.5));
  connect(ConstantSpeed1.flange, Compressor.shaft_a) annotation (Line(
      points={{-30,0},{-30,0},{-26,-0.2},{-12,0}},
      color={0,0,0},
      thickness=0.5));
  annotation (
    Diagram(graphics),
    experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-2, Interval = 1),
    __OpenModelica_commandLineOptions = "--matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection -d=initialization,NLSanalyticJacobian,aliasConflicts,bltdump",    
    __OpenModelica_simulationFlags(lv = "LOG_DEBUG,LOG_NLS,LOG_NLS_V,LOG_STATS,LOG_INIT,LOG_STDOUT, -w", outputFormat = "mat", s = "dassl", nls = "homotopy"));
end TestTP_Compressor;
