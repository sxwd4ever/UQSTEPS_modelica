within Steps.Test.TPComponents;

model TestTP_HE
  "Test against HE as a Heater in ThermoPower"  
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
  
  // package medium_main = Steps.Media.CO2;
  // package medium_heater = SolarTherm.Media.Sodium.Sodium_pT; // ThermoPower.Water.StandardWater; //Modelica.Media.IdealGases.SingleGases.CO2;
  package medium_heater = Steps.Media.MoltenSalt.MoltenSalt_pT;
  // package medium_heater = Steps.Media.ThermiaOilD;
  // package medium_cold = Steps.Media.CO2; // Modelica.Media.Ideal//Gases.SingleGases.CO2;
  package medium_main = Steps.Media.SCO2;//ExternalMedia.Examples.CO2CoolProp;
  package medium_cold = Steps.Media.SCO2;//ExternalMedia.Examples.CO2CoolProp;    
  
  parameter Model.RCBCycleConfig cfg(
    redeclare package medium_heater = medium_heater,
    redeclare package medium_main = medium_main,
    // mdot_heater      = 40,
    // T_heater_hot_in  = from_degC(800),
    // T_heater_hot_out = from_degC(600),
    r_i_heater       = 20e-3,
    r_t_heater       = cfg.r_i_heater + 10e-3,
    r_o_heater       = 1/2,                      // agree with the final parameter Dhyd = 1 in HE, should be checked to see if it is capable of containing all fluid-metal tubes
    N_ch_heater      = 100,
    L_heater         = 1
  );
 
  parameter Model.RCBCycleConfig cfg_test(
    redeclare package medium_heater = medium_heater,
    redeclare package medium_main = medium_main,  
    mdot_main       = 93.75,
    mdot_heater     = 55,
    T_heater_hot_in = from_degC(550),
    L_heater        = 1
  );        

  parameter Model.HeatExchangerConfig cfg_heater = cfg.cfg_heater;
  parameter Model.ThermoState st_source_hot      = cfg_heater.cfg_hot.st_in;
  parameter Model.ThermoState st_sink_hot        = cfg_heater.cfg_hot.st_out;
  parameter Model.ThermoState st_source_cold     = cfg_heater.cfg_cold.st_in;
  parameter Model.ThermoState st_sink_cold       = cfg_heater.cfg_cold.st_out;
  
  ThermoPower.Water.SourceMassFlow source_heater_hot(
  redeclare package Medium = medium_heater,
    w0 = st_source_hot.mdot,
    p0 = st_source_hot.p,
    h  = st_source_hot.h,
    T  = st_source_hot.T
    // use_T = true,
    // use_in_T = false 
  ) 
  annotation(
    Placement(transformation(origin = {0, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 270)));
    
  ThermoPower.Water.SinkPressure sink_heater_hot(
  redeclare package Medium = medium_heater,
    p0    = st_sink_hot.p,
    T     = st_sink_hot.T,
    h     = st_sink_hot.h,
    use_T = true
  ) 
  annotation(
    Placement(transformation(origin = {0, -80}, extent = {{-10, -10}, {10, 10}}, rotation = 270)));
    
  ThermoPower.Gas.SourceMassFlow source_cold(
  redeclare package Medium = medium_cold,
    T        = st_source_cold.T,
    p0       = st_source_cold.p,
    use_in_T = false,
    w0       = st_source_cold.mdot,
    gas(
      p(nominal = st_source_cold.p),
      T(nominal = st_source_cold.T)
    )
  )
  annotation(
    Placement(transformation(origin = {0, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 270)));      
  
  ThermoPower.Gas.SinkPressure sink_cold(
  redeclare package Medium = medium_cold,
    p0 = st_sink_cold.p,
    T  = st_sink_cold.T
  )
  annotation(
    Placement(transformation(origin = {0, -80}, extent = {{-10, -10}, {10, 10}}, rotation = 270)));  
  
  Steps.TPComponents.HE hE(
  //Components.HEG2G hE(
    redeclare package FluidMedium              = medium_heater,
    redeclare package FlueGasMedium            = medium_cold,
    redeclare replaceable model HeatTransfer_F = ThermoPower.Thermal.HeatTransferFV.IdealHeatTransfer,             //ConstantHeatTransferCoefficientTwoGrids(gamma = thermo_hot.gamma_he),     
    redeclare replaceable model HeatTransfer_G = ThermoPower.Thermal.HeatTransferFV.IdealHeatTransfer,             //ConstantHeatTransferCoefficient(gamma =  thermo_cold.gamma_he),     
    redeclare model HeatExchangerTopology      = ThermoPower.Thermal.HeatExchangerTopologies.CounterCurrentFlow,
    fluidFlow(
      fixedMassFlowSimplified = true,              // set the fluid flow as fixed mdot for simplarity
      hstartin                = st_source_hot.h,
      hstartout               = st_sink_hot.h
    ),
    gasFlow(
      QuasiStatic = true,
      Tstartin    = st_source_cold.T,
      Tstartout   = st_sink_cold.T
    ),
    cfg             = cfg_heater,
    N_G             = cfg.gp_heater_hot.N_seg,                             // N_G and N_F has to be assigned explicitly since they set dimension of Gas and Fluid cells as gasflow.gas[N_G] and fluidflow.gas[N_F]
    N_F             = cfg.gp_heater_hot.N_seg,
    SSInit          = true,
    FluidPhaseStart = ThermoPower.Choices.FluidPhase.FluidPhases.Liquid,
    metalTube(
      WallRes   = false,
      Tstartbar = st_source_hot.T - 50
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
  
  ThermoPower.Water.SensT T_waterIn(redeclare package Medium = medium_heater) annotation(
    Placement(transformation(origin = {4, -50}, extent = {{-10, -10}, {10, 10}}, rotation = 270)));
  ThermoPower.Water.SensT T_waterOut(redeclare package Medium = medium_heater) annotation(
    Placement(transformation(origin = {4, -50}, extent = {{-10, -10}, {10, 10}}, rotation = 270)));
  
  ThermoPower.Gas.SensT T_gasIn(redeclare package Medium = medium_cold) annotation(
    Placement(transformation(extent = {{30, -6}, {50, 14}}, rotation = 0)));
  ThermoPower.Gas.SensT T_gasOut(redeclare package Medium = medium_cold) annotation(
    Placement(transformation(extent = {{30, -6}, {50, 14}}, rotation = 0)));
  
  // variable for validation
  Modelica.SIunits.Power Q_out    = (hE.gasIn.h_outflow - hE.gasOut.h_outflow) * hE.gasIn.m_flow;
  Modelica.SIunits.Power Q_in     = (hE.waterOut.h_outflow - hE.waterIn.h_outflow) * hE.waterIn.m_flow;
  Boolean                isQMatch = abs(Q_out -Q_in) < 1e-3;
initial equation
//hstart_F_Out = hE.waterOut.h_outflow;
equation
  connect(source_heater_hot.flange, T_waterIn.inlet);
  connect(T_waterIn.outlet, hE.waterIn) annotation(
    Line(points = {{-1.83697e-015, 50}, {-1.83697e-015, 20}, {0, 20}}, color = {0, 0, 255}, thickness = 0.5, smooth = Smooth.None));
  connect(hE.waterOut, T_waterOut.inlet);
  connect(T_waterOut.outlet, sink_heater_hot.flange) annotation(
    Line(points = {{1.83697e-015, -70}, {1.83697e-015, -56}, {-8.88178e-016, -56}}, thickness = 0.5, color = {0, 0, 255}));
   
  connect(source_cold.flange, T_gasIn.inlet);
  connect(T_gasIn.outlet, hE.gasIn) annotation(
    Line(points = {{-50, 0}, {-20, 0}}, color = {159, 159, 223}, thickness = 0.5, smooth = Smooth.None));
  connect(hE.gasOut, T_gasOut.inlet);
  connect(T_gasOut.outlet, sink_cold.flange) annotation(
    Line(points = {{46, 0}, {46, 0}, {60, 0}}, color = {159, 159, 223}, thickness = 0.5));

annotation(
    Diagram(graphics),
    experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-2, Interval = 1),
    __OpenModelica_commandLineOptions = "--matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection -d=initialization,NLSanalyticJacobian,aliasConflicts,bltdump",    
    __OpenModelica_simulationFlags(lv = "LOG_DEBUG,LOG_NLS,LOG_NLS_V,LOG_STATS,LOG_INIT,LOG_STDOUT, -w", outputFormat = "mat", s = "dassl", nls = "homotopy"));
end TestTP_HE;
