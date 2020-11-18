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
  
  // package Medium = Steps.Media.CO2; 
  // package Medium = ExternalMedia.Examples.CO2CoolProp;
  package Medium = Steps.Media.SCO2; 
  // package Medium = Modelica.Media.IdealGases.SingleGases.CO2;
  
  parameter Model.PBConfiguration cfg_default;
  
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
  parameter Model.PBConfiguration cfg = cfg_default;
  
    // set the values of parameters accordingly
  parameter HEBoundaryCondition bc_LTR = cfg.bc_LTR;  
  parameter HEBoundaryCondition bc_cooler = cfg.bc_cooler;
  parameter Real N_des = 523.3;
 
  //configuration for re-compressor 
  parameter Modelica.SIunits.Temperature T_comp_des = from_degC(90) "Design temperature for recompressor";   
  
  ThermoPower.Gas.SourcePressure SourceP1(
    redeclare package Medium = Medium,
    p0=bc_LTR.st_hot_out.p,
    T=bc_LTR.st_hot_out.T,
    gas(p(nominal = SourceP1.p0), 
    T(nominal=SourceP1.T))) annotation (Placement(transformation(extent={{-80,6},{-60,26}},
          rotation=0)));
  
  ThermoPower.Gas.SinkPressure SinkP1(
    redeclare package Medium = Medium,
    p0=bc_LTR.st_cold_out.p,
    T=bc_LTR.st_cold_out.T,
    gas(
    p(nominal = SinkP1.p0), 
    T(nominal = SinkP1.T))) annotation (Placement(transformation(extent={{40,6},{60,26}},
          rotation=0)));
        
  ThermoPower.Gas.Compressor Compressor(
    redeclare package Medium = Medium,
    pstart_in=bc_LTR.st_hot_out.p,
    pstart_out=bc_LTR.st_cold_out.p,
    Tstart_in=bc_LTR.st_hot_out.T,
    Tstart_out=bc_LTR.st_cold_out.T,
    tablePhic=tablePhic_rc,
    tableEta=tableEta_rc,
    tablePR=tablePR_rc,
    Table=ThermoPower.Choices.TurboMachinery.TableTypes.matrix,
    Ndesign=N_des,
    Tdes_in=T_comp_des,
    explicitIsentropicEnthalpy = false,
    gas_in(
      p(nominal = Compressor.pstart_in), 
      T(nominal = Compressor.Tstart_in)),
    gas_iso(
      p(nominal = Compressor.pstart_out), 
      T(nominal = Compressor.Tstart_out))) annotation (Placement(transformation(extent={{-20,-20},{
              20,20}}, rotation=0)));
  /*
  //configuration for main compressor          
  parameter Modelica.SIunits.Temperature T_comp_des = from_degC(45) "Design temperature for main compressor";
  
  ThermoPower.Gas.SourcePressure SourceP1(
    redeclare package Medium = Medium,
    p0=bc_cooler.st_hot_out.p,
    T=bc_cooler.st_hot_out.T,
    gas(p(nominal = SourceP1.p0), 
    T(nominal=SourceP1.T))) annotation (Placement(transformation(extent={{-80,6},{-60,26}},
          rotation=0)));
  
  ThermoPower.Gas.SinkPressure SinkP1(
    redeclare package Medium = Medium,
    p0=bc_LTR.st_cold_in.p,
    T=bc_LTR.st_cold_in.T,
    gas(
    p(nominal = SinkP1.p0), 
    T(nominal = SinkP1.T))) annotation (Placement(transformation(extent={{40,6},{60,26}},
          rotation=0)));
  
  ThermoPower.Gas.Compressor Compressor(
    redeclare package Medium = Medium,
    pstart_in=bc_cooler.st_hot_out.p,
    pstart_out=bc_LTR.st_cold_in.p,
    Tstart_in=bc_cooler.st_hot_out.T,
    Tstart_out=bc_LTR.st_cold_in.T,
    tablePhic=tablePhic_mc,
    tableEta=tableEta_mc,
    tablePR=tablePR_mc,
    Table=ThermoPower.Choices.TurboMachinery.TableTypes.matrix,
    Ndesign=N_des,
    Tdes_in=T_comp_des,
    explicitIsentropicEnthalpy = false,
    gas_in(
      p(nominal = Compressor.pstart_in), 
      T(nominal = Compressor.Tstart_in)),
    gas_iso(
      p(nominal = Compressor.pstart_out), 
      T(nominal = Compressor.Tstart_out))) annotation (Placement(transformation(extent={{-20,-20},{
              20,20}}, rotation=0)));  
  */          
  Modelica.Mechanics.Rotational.Sources.ConstantSpeed ConstantSpeed1(
      w_fixed=N_des, useSupport=false) annotation (Placement(transformation(
          extent={{-50,-10},{-30,10}}, rotation=0)));
  inner ThermoPower.System system
    annotation (Placement(transformation(extent={{80,80},{100,100}})));
protected
  // performance map for main compressor
  parameter Real tableEta_mc[5, 4]=[0, 95, 100, 105; 1, 0.85310219, 0.837591241, 0.832420925; 2, 0.868613139, 0.857238443, 0.851034063; 3, 0.860340633, 0.85, 0.842761557; 4, 0.85310219, 0.839659367, 0.816909976];
  
  parameter Real tablePhic_mc[5, 4]=[0, 95, 100, 105; 1, 0.000134346, 0.000150832, 0.000164161; 2, 0.000137854, 0.000153638, 0.00016802; 3, 0.000142414, 0.000158549, 0.000169774; 4, 0.000145921, 0.000161706, 0.000171528];

  parameter Real tablePR_mc[5, 4]=[0, 95, 100, 105; 1, 1.967529638, 2.350588505, 2.785882673; 2, 1.915294338, 2.315764972, 2.681412073; 3, 1.810823737, 2.220000255, 2.524706172; 4, 1.654117837, 2.115529655, 2.359294389];

  // performance map for re compressor      
  parameter Real tableEta_rc[5, 4]=[0, 95, 100, 105; 1, 0.85310219, 0.837591241, 0.832420925; 2, 0.868613139, 0.857238443, 0.851034063; 3, 0.860340633, 0.85, 0.842761557; 4, 0.85310219, 0.839659367, 0.816909976];
  parameter Real tablePhic_rc[5, 4]=[0, 95, 100, 105; 1, 7.17663E-05, 8.05731E-05, 8.76935E-05; 2, 7.36401E-05, 8.20721E-05, 8.97547E-05; 3, 7.6076E-05, 8.46954E-05, 9.06916E-05; 4, 7.79498E-05, 8.63819E-05, 9.16285E-05];

  parameter Real tablePR_rc[5, 4]=[0, 95, 100, 105; 1, 1.967529638, 2.350588505, 2.785882673; 2, 1.915294338, 2.315764972, 2.681412073; 3, 1.810823737, 2.220000255, 2.524706172; 4, 1.654117837, 2.115529655, 2.359294389];   
  
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
