within Steps.Test.TPComponents;

model TestTP_Components_Comp_Path
  "Test for combination of components (Cooler + Spliter + Main/Re compressor) in ThermoPower"  
    
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
  
  package medium_hot = Steps.Media.CO2;
  package medium_cooler = ThermoPower.Water.StandardWater;//Modelica.Media.IdealGases.SingleGases.CO2;  
  
  // parameter for C++ implementation of PCHE - based on Modelica impl's result
  parameter Model.PBConfiguration cfg_tune(
    mdot_cooler = 40,
    r_i_c = 15e-3,
    r_t_c = cfg.r_i_c + 5e-3,
    r_o_c = cfg.r_t_c + 200e-3,
    L_c = 20);  
 
  // default configuration of parameters
  parameter Model.PBConfiguration cfg;
  
  // set the values of parameters accordingly
  parameter HEBoundaryCondition bc_cooler = cfg.bc_cooler;
  parameter HEBoundaryCondition bc_LTR = cfg.bc_LTR;
  parameter HEBoundaryCondition bc_HTR = cfg.bc_HTR;
  
  parameter ThermoState st_bypass = cfg.st_bypass;  
  
  ThermoPower.Gas.SourceMassFlow source(
  redeclare package Medium = medium_hot,
    w0 = bc_LTR.st_hot_out.mdot,
    p0 = bc_LTR.st_hot_out.p,
    T = bc_LTR.st_hot_out.T) 
    annotation(
    Placement(transformation(origin = {0, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 270)));
    
  ThermoPower.Gas.SinkPressure sink(
  redeclare package Medium = medium_hot,
    p0 = bc_LTR.st_cold_in.p,
    T =bc_LTR.st_cold_in.T) 
    annotation(
    Placement(transformation(origin = {0, -80}, extent = {{-10, -10}, {10, 10}}, rotation = 270)));
    
  ThermoPower.Water.SourceMassFlow source_cooler_cold(
  redeclare package Medium = medium_cooler,
    T = bc_cooler.st_cold_in.T, 
    p0 = bc_cooler.st_cold_in.p, 
    use_T = true,
    use_in_T = false, 
    w0 = bc_cooler.st_cold_in.mdot) 
  annotation(
    Placement(transformation(origin = {0, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 270)));      
  
  ThermoPower.Water.SinkPressure sink_cooler_cold(
  redeclare package Medium = medium_cooler,
    p0 = bc_cooler.st_cold_out.p, 
    T = bc_cooler.st_cold_out.T,
    use_T = true) 
  annotation(
    Placement(transformation(origin = {0, -80}, extent = {{-10, -10}, {10, 10}}, rotation = 270)));  
  
  ThermoPower.Water.SinkPressure sink_bypass(
  redeclare package Medium = medium_cooler,
    p0 = st_bypass.p, 
    T = st_bypass.T,
    use_T = true) 
  annotation(
    Placement(transformation(origin = {0, -80}, extent = {{-10, -10}, {10, 10}}, rotation = 270)));    
  
  TPComponents.HE cooler(
  //Components.HEG2G cooler(
    redeclare package FluidMedium = medium_cooler, 
    redeclare package FlueGasMedium = medium_hot, 
    fluidFlow(fixedMassFlowSimplified = true, hstartin = bc_cooler.st_cold_in.h, hstartout=bc_cooler.st_cold_out.h), // set the fluid flow as fixed mdot for simplarity
    gasFlow(QuasiStatic = true,Tstartin = bc_cooler.st_hot_in.T, Tstartout = bc_cooler.st_hot_out.T),
    bc = bc_cooler, 
    geo_hot = cfg.cfg_cooler_hot.geo,
    geo_cold = cfg.cfg_cooler_cold.geo,
    geo_tube = cfg.cfg_cooler_tube.geo,  
    thermo_hot = cfg.cfg_cooler_hot.thermo,
    thermo_cold = cfg.cfg_cooler_cold.thermo,
    thermo_tube = cfg.cfg_cooler_tube.thermo, 
    SSInit=true,
    FluidPhaseStart=ThermoPower.Choices.FluidPhase.FluidPhases.Liquid,    
    redeclare replaceable model HeatTransfer_F = ThermoPower.Thermal.HeatTransferFV.IdealHeatTransfer, //ConstantHeatTransferCoefficientTwoGrids(gamma = thermo_hot.gamma_he),     
    redeclare replaceable model HeatTransfer_G = ThermoPower.Thermal.HeatTransferFV.IdealHeatTransfer, //ConstantHeatTransferCoefficient(gamma =  thermo_cold.gamma_he),     
    redeclare model HeatExchangerTopology = ThermoPower.Thermal.HeatExchangerTopologies.CounterCurrentFlow,
    metalTube(WallRes=false, Tstartbar=bc_cooler.st_hot_in.T - 50)) annotation(
    Placement(transformation(extent = {{-20, -20}, {20, 20}}, rotation = 0)));

  inner ThermoPower.System system(allowFlowReversal = false, initOpt=ThermoPower.Choices.Init.Options.noInit) annotation(
    Placement(transformation(extent = {{80, 80}, {100, 100}})));
  
  //ThermoPower.Water.SensT T_waterOut(redeclare package Medium = medium_heater) annotation(
  //  Placement(transformation(origin = {4, -50}, extent = {{-10, -10}, {10, 10}}, rotation = 270)));
  ThermoPower.Gas.SensT T_gasOut(redeclare package Medium = medium_hot) annotation(
    Placement(transformation(extent = {{30, -6}, {50, 14}}, rotation = 0)));
  
  ThermoPower.Gas.FlowSplit splitter(redeclare package Medium = medium_hot
  );
  
  ThermoPower.Gas.Compressor compressor(
    redeclare package Medium = medium_hot,
    pstart_in=bc_cooler.st_hot_out.p,
    pstart_out=bc_LTR.st_cold_in.p,
    Tstart_in=bc_cooler.st_hot_out.T,
    Tstart_out=bc_LTR.st_cold_in.T,
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

  ThermoPower.Gas.Compressor compressor_bypass(
    redeclare package Medium = medium_hot,
    pstart_in=bc_HTR.st_hot_out.p,
    pstart_out=st_bypass.p,
    Tstart_in=bc_HTR.st_hot_out.T,
    Tstart_out=st_bypass.T,
    tablePhic=tablePhic,
    tableEta=tableEta,
    tablePR=tablePR,
    Table=ThermoPower.Choices.TurboMachinery.TableTypes.matrix,
    Ndesign=523.3,
    Tdes_in=244.4) annotation (Placement(transformation(extent={{-20,-20},{
            20,20}}, rotation=0)));
            
  Modelica.Mechanics.Rotational.Sources.ConstantSpeed ConstantSpeed2(
      w_fixed=523.3, useSupport=false) annotation (Placement(transformation(
          extent={{-50,-10},{-30,10}}, rotation=0)));

protected
  parameter Real tableEta[6, 4]=[0, 95, 100, 105; 1, 82.5e-2, 81e-2,
      80.5e-2; 2, 84e-2, 82.9e-2, 82e-2; 3, 83.2e-2, 82.2e-2, 81.5e-2; 4,
      82.5e-2, 81.2e-2, 79e-2; 5, 79.5e-2, 78e-2, 76.5e-2];
  parameter Real tablePhic[6, 4]=[0, 95, 100, 105; 1, 38.3e-3/400, 43e-3/400,
      46.8e-3/400; 2, 39.3e-3/400, 43.8e-3/400, 47.9e-3/400; 3, 40.6e-3/400, 45.2e-3/400, 48.4e-3/400;
      4, 41.6e-3/400, 46.1e-3/400, 48.9e-3/400; 5, 42.3e-3/400, 46.6e-3/400, 49.3e-3/400];

  parameter Real tablePR[6, 4]=[0, 95, 100, 105; 1, 22.6/12, 27/12, 32/12; 2, 22/12,
      26.6/12, 30.8/12; 3, 20.8/12, 25.5/12, 29/12; 4, 19/12, 24.3/12, 27.1/12; 5, 17/12, 21.5/12, 24.2/12];


initial equation
//hstart_F_Out = cooler.waterOut.h_outflow;
equation
  
  // Cooler + Main compressor + recompressor + splitter
  connect(source.flange, splitter.inlet);
  connect(splitter.outlet1, cooler.gasIn) annotation(
    Line(points = {{-1.83697e-015, 50}, {-1.83697e-015, 20}, {0, 20}}, color = {0, 0, 255}, thickness = 0.5, smooth = Smooth.None));
  
  connect(cooler.gasOut, T_gasOut.inlet);
  
  connect(T_gasOut.outlet, compressor.inlet);
  connect(compressor.shaft_a, ConstantSpeed1.flange) annotation (Line(
      points={{-30,0},{-30,0},{-26,-0.2},{-12,0}},
      color={0,0,0},
      thickness=0.5));
  connect(compressor.outlet, sink.flange) annotation(
    Line(points = {{1.83697e-015, -70}, {1.83697e-015, -56}, {-8.88178e-016, -56}}, thickness = 0.5, color = {0, 0, 255}));
  
  // bypass path 
  connect(splitter.outlet2, compressor_bypass.inlet);  
  connect(compressor_bypass.outlet, sink_bypass.flange); 
  connect(compressor_bypass.shaft_a, ConstantSpeed2.flange) annotation (Line(
      points={{-30,0},{-30,0},{-26,-0.2},{-12,0}},
      color={0,0,0},
      thickness=0.5));  
    
  // water cooling path
  connect(source_cooler_cold.flange, cooler.waterIn) annotation(
    Line(points = {{-50, 0}, {-20, 0}}, color = {159, 159, 223}, thickness = 0.5, smooth = Smooth.None));
  connect(cooler.waterOut, sink_cooler_cold.flange) annotation(
    Line(points = {{46, 0}, {46, 0}, {60, 0}}, color = {159, 159, 223}, thickness = 0.5));

annotation(
    Diagram(graphics),
    experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-3, Interval = 1),
    __OpenModelica_commandLineOptions = "--matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection -d=initialization,NLSanalyticJacobian,aliasConflicts,bltdump",    
    __OpenModelica_simulationFlags(lv = "LOG_DEBUG,LOG_NLS,LOG_NLS_V,LOG_STATS,LOG_INIT,LOG_STDOUT, -w", outputFormat = "mat", s = "dassl", nls = "homotopy"));
end TestTP_Components_Comp_Path;
