within Steps.Test.TPComponents;

model TestTP_HE
  "Test against HE as a Heater or Cooler in ThermoPower"  
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
  
  // Alternative mediums  
  // package medium_HE = SolarTherm.Media.Sodium.Sodium_pT; 
  // package medium_HE = ThermoPower.Water.StandardWater; 
  // package medium_HE = Modelica.Media.IdealGases.SingleGases.CO2;
  // package medium_HE = Steps.Media.ThermiaOilD;
  // package medium_main = Steps.Media.CO2;
  // package medium_main = Modelica.Media.Ideal; 
  // package medium_main = Gases.SingleGases.CO2; 
  // package medium_main = Steps.Media.SCO2;
  // package medium_main = ExternalMedia.Examples.CO2CoolProp;   
     
  // Two paths in this test
  // main path -> main stream of the working fluid <-> gas side of HE
  // HE path -> heater's hot side or cooler's cold side, <-> fluid side of HE
 
    
  parameter Real Cf_C1 = 1.626, Cf_C2 = 1, Cf_C3 = 1;
  
  // configs for heater test
  parameter Model.RCBCycleConfig cfg(   
    // results calculated at 2021-05-21 20:07, RCBC without recompressor, Open LOOP
    p_comp_in  = 109.59e5,
    p_comp_out = 20e6,    
    p_heater   = 20e6,    
    T_HTR_hot_in      = from_degC(556.322),
    T_HTR_cold_out    = from_degC(521.234),
    T_HTR_hot_out     = from_degC(330.103),
    T_HTR_cold_in     = from_degC(303.425),
    T_LTR_cold_in     = from_degC(119.011),
    T_LTR_hot_out     = from_degC(164.419),
    T_heater_hot_in   = from_degC(800),
    T_heater_hot_out  = from_degC(523.547),
    T_heater_cold_out = from_degC(608.148),
    T_cooler_cold_out = from_degC(112.138),
    T_cooler_hot_out  = from_degC(59.279),
    
    mdot_main   = 128.774,    
    mdot_comp   = 88.0661,    
    mdot_heater = 40,
    mdot_cooler = 40.7188  
  );    
  /*
  // For Heater test    
  package medium_HE   = SolarTherm.Media.Sodium.Sodium_pT;
  package medium_main = Steps.Media.SCO2;
    
  parameter Model.HeatExchangerConfig cfg_HE  = cfg.cfg_heater;
  
  parameter Model.ThermoState st_source_main  = cfg_HE.cfg_cold.st_in;
  parameter Model.ThermoState st_sink_main    = cfg_HE.cfg_cold.st_out;
  parameter Model.ThermoState st_source_HE    = cfg_HE.cfg_hot.st_in;
  parameter Model.ThermoState st_sink_HE      = cfg_HE.cfg_hot.st_out;
  
  parameter Integer N_seg_HE                  = cfg.cfg_heater.cfg_hot.geo_path.N_seg;  
  */
  
  
  // For Cooler test
  package medium_HE = ThermoPower.Water.StandardWater;
  package medium_main = Steps.Media.SCO2;         

  parameter Model.HeatExchangerConfig cfg_HE = cfg.cfg_cooler;
  
  parameter Model.ThermoState st_source_main  = cfg_HE.cfg_hot.st_in;
  parameter Model.ThermoState st_sink_main    = cfg_HE.cfg_hot.st_out;
  parameter Model.ThermoState st_source_HE = cfg_HE.cfg_cold.st_in;
  parameter Model.ThermoState st_sink_HE   = cfg_HE.cfg_cold.st_out;
  
  parameter Integer N_seg_HE                 = cfg.cfg_cooler.cfg_hot.geo_path.N_seg;
  
    
  ThermoPower.Water.SourceMassFlow source_HE(
  redeclare package Medium = medium_main, //medium_HE,
    w0 = st_source_HE.mdot,
    p0 = st_source_HE.p,
    h  = st_source_HE.h,
    T  = st_source_HE.T
    // use_T = true,
    // use_in_T = false 
  ) 
  annotation(
    Placement(transformation(origin = {0, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 270)));
    
  ThermoPower.Water.SinkPressure sink_HE(
  redeclare package Medium = medium_main, //medium_HE,
    p0    = st_sink_HE.p,
    T     = st_sink_HE.T,
    h     = st_sink_HE.h,
    use_T = false
  ) 
  annotation(
    Placement(transformation(origin = {0, -80}, extent = {{-10, -10}, {10, 10}}, rotation = 270)));
    
  ThermoPower.Gas.SourceMassFlow source_main(
  redeclare package Medium = medium_main, //medium_main,
    T        = st_source_main.T,
    p0       = st_source_main.p,
    use_in_T = false,
    w0       = st_source_main.mdot,
    gas(
      p(nominal = st_source_main.p),
      T(nominal = st_source_main.T)
    )
  )
  annotation(
    Placement(transformation(origin = {0, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 270)));      
  
  ThermoPower.Gas.SinkPressure sink_main(
  redeclare package Medium = medium_main, //medium_main,
    p0 = st_sink_main.p,
    T  = st_sink_main.T
  )
  annotation(
    Placement(transformation(origin = {0, -80}, extent = {{-10, -10}, {10, 10}}, rotation = 270)));  
  
  
  // main side = gas side
  // HE side = fluid side 
  Steps.TPComponents.HE hE(
  //Components.HEG2G hE(
    redeclare package FluidMedium              = medium_HE,
    redeclare package FlueGasMedium            = medium_main,
    
    // Ngo 
    redeclare replaceable model HeatTransfer_F = Steps.TPComponents.NgoHeatTransferFV(),
    redeclare replaceable model HeatTransfer_G = Steps.TPComponents.NgoHeatTransferFV(),
    gasFlow(      
      heatTransfer(
        Cf_C1 = Cf_C1,
        Cf_C2 = Cf_C2,
        Cf_C3 = Cf_C3,
        gamma_min = 1000,
        gamma_max = 5000        
      ),
      Tstartin  = cfg_HE.cfg_gas.st_in.T,
      Tstartout = cfg_HE.cfg_gas.st_out.T,
      Nt        = cfg_HE.cfg_gas.N_ch,
      Dhyd      = cfg_HE.cfg_gas.geo_area.d_hyd
    ),
    fluidFlow(      
      heatTransfer(
        Cf_C1 = Cf_C1,
        Cf_C2 = Cf_C2,
        Cf_C3 = Cf_C3,
        gamma_min = 2000,
        gamma_max = 1.5e6
      ), 
      fixedMassFlowSimplified = true,
      hstartin                = cfg_HE.cfg_fluid.st_in.h,
      hstartout               = cfg_HE.cfg_fluid.st_out.h
    ),     
    /*    
    redeclare replaceable model HeatTransfer_F = ThermoPower.Thermal.HeatTransferFV.IdealHeatTransfer,             //ConstantHeatTransferCoefficientTwoGrids(gamma = thermo_hot.gamma_he),     
    redeclare replaceable model HeatTransfer_G = ThermoPower.Thermal.HeatTransferFV.IdealHeatTransfer,             //ConstantHeatTransferCoefficient(gamma =  thermo_cold.gamma_he),     
    redeclare model HeatExchangerTopology      = ThermoPower.Thermal.HeatExchangerTopologies.CounterCurrentFlow,
    fluidFlow(
      fixedMassFlowSimplified = true,              // set the fluid flow as fixed mdot for simplarity
      hstartin                = st_source_HE.h,
      hstartout               = st_sink_HE.h
    ),
    gasFlow(
      QuasiStatic = true,
      Tstartin    = st_source_main.T,
      Tstartout   = st_sink_main.T
    ),
    */
    cfg             = cfg_HE,
    N_G             = N_seg_HE,                             // N_G and N_F has to be assigned explicitly since they set dimension of Gas and Fluid cells as gasflow.gas[N_G] and fluidflow.gas[N_F]
    N_F             = N_seg_HE,
    SSInit          = false,
    FluidPhaseStart = ThermoPower.Choices.FluidPhase.FluidPhases.Liquid,
    metalTube(
      WallRes   = false,
      Tstartbar = st_source_main.T - 50
    )
  )
  annotation(
    Placement(transformation(extent = {{-20, -20}, {20, 20}}, rotation = 0)));

  inner ThermoPower.System system(
    allowFlowReversal = false,
    initOpt           = ThermoPower.Choices.Init.Options.noInit
    ) 
  annotation(
    Placement(transformation(extent = {{80, 80}, {100, 100}})));
  
  ThermoPower.Water.SensT T_waterIn(redeclare package Medium = medium_HE) annotation(
    Placement(transformation(origin = {4, -50}, extent = {{-10, -10}, {10, 10}}, rotation = 270)));
  ThermoPower.Water.SensT T_waterOut(redeclare package Medium = medium_HE) annotation(
    Placement(transformation(origin = {4, -50}, extent = {{-10, -10}, {10, 10}}, rotation = 270)));
  
  ThermoPower.Gas.SensT T_gasIn(redeclare package Medium = medium_main) annotation(
    Placement(transformation(extent = {{30, -6}, {50, 14}}, rotation = 0)));
  ThermoPower.Gas.SensT T_gasOut(redeclare package Medium = medium_main) annotation(
    Placement(transformation(extent = {{30, -6}, {50, 14}}, rotation = 0)));
  
  // variable for validation
  Modelica.SIunits.Power Q_out    = (hE.gasIn.h_outflow - hE.gasOut.h_outflow) * hE.gasIn.m_flow;
  Modelica.SIunits.Power Q_in     = (hE.waterOut.h_outflow - hE.waterIn.h_outflow) * hE.waterIn.m_flow;
  Boolean                isQMatch = abs(Q_out -Q_in) < 1e-3;
initial equation
//hstart_F_Out = hE.waterOut.h_outflow;
equation

/* 
  // Heater test   
  connect(source_main.flange, T_waterIn.inlet);
  connect(T_waterIn.outlet, hE.waterIn) annotation(
    Line(points = {{-1.83697e-015, 50}, {-1.83697e-015, 20}, {0, 20}}, color = {0, 0, 255}, thickness = 0.5, smooth = Smooth.None));
  connect(hE.waterOut, T_waterOut.inlet);
  connect(T_waterOut.outlet, sink_main.flange) annotation(
    Line(points = {{1.83697e-015, -70}, {1.83697e-015, -56}, {-8.88178e-016, -56}}, thickness = 0.5, color = {0, 0, 255}));
   
  connect(source_HE.flange, T_gasIn.inlet);
  connect(T_gasIn.outlet, hE.gasIn) annotation(
    Line(points = {{-50, 0}, {-20, 0}}, color = {159, 159, 223}, thickness = 0.5, smooth = Smooth.None));
  connect(hE.gasOut, T_gasOut.inlet);
  connect(T_gasOut.outlet, sink_HE.flange) annotation(
    Line(points = {{46, 0}, {46, 0}, {60, 0}}, color = {159, 159, 223}, thickness = 0.5));
*/

  // Coolear test  
  connect(source_HE.flange, T_waterIn.inlet);
  connect(T_waterIn.outlet, hE.waterIn) annotation(
    Line(points = {{-1.83697e-015, 50}, {-1.83697e-015, 20}, {0, 20}}, color = {0, 0, 255}, thickness = 0.5, smooth = Smooth.None));
  connect(hE.waterOut, T_waterOut.inlet);
  connect(T_waterOut.outlet, sink_HE.flange) annotation(
    Line(points = {{1.83697e-015, -70}, {1.83697e-015, -56}, {-8.88178e-016, -56}}, thickness = 0.5, color = {0, 0, 255}));
   
  connect(source_main.flange, T_gasIn.inlet);
  connect(T_gasIn.outlet, hE.gasIn) annotation(
    Line(points = {{-50, 0}, {-20, 0}}, color = {159, 159, 223}, thickness = 0.5, smooth = Smooth.None));
  connect(hE.gasOut, T_gasOut.inlet);
  connect(T_gasOut.outlet, sink_main.flange) annotation(
    Line(points = {{46, 0}, {46, 0}, {60, 0}}, color = {159, 159, 223}, thickness = 0.5));


annotation(
    Diagram(graphics),
    experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-2, Interval = 1),
    __OpenModelica_commandLineOptions = "--matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection -d=initialization,NLSanalyticJacobian,aliasConflicts",    
    __OpenModelica_simulationFlags(lv = "LOG_DEBUG,LOG_NLS,LOG_NLS_V,LOG_STATS,LOG_INIT,LOG_STDOUT, -w", outputFormat = "mat", s = "dassl", nls = "homotopy"));
end TestTP_HE;
