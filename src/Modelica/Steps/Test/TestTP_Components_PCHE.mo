within Steps.Test;

model TestTP_Components_PCHE
  "Test for combination of components: heater + turbine + HTR + LTR in ThermoPower, PCHE as recuperator"  
    
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
  
  //package medium_hot = Modelica.Media.IdealGases.SingleGases.CO2; //Steps.Media.CO2;
  //package medium_cold = Modelica.Media.IdealGases.SingleGases.CO2; //Steps.Media.CO2;
  // package medium_hot = Steps.Media.CO2;
  // package medium_cold = Steps.Media.CO2;
  package medium_hot = Steps.Media.SCO2;
  package medium_cold = Steps.Media.SCO2;  
  package medium_heater = SolarTherm.Media.Sodium.Sodium_pT;
  //package medium_heater = ThermoPower.Water.StandardWater;// Modelica.Media.IdealGases.SingleGases.CO2;

  parameter Model.PBConfiguration cfg_default(
    mdot_heater = 80
  );
  
  parameter Model.PBConfiguration cfg_tune(
  r_i_HTR = 5e-3,
  r_t_HTR = cfg_tune.r_i_HTR + 2e-3,
  r_o_HTR = cfg_tune.r_t_HTR + 2e-3,
  N_ch_HTR = 200,
  L_HTR = 1 "Don't modify this, since L in HE model is fixed as 1m. Modify Nt instead",
  r_i_LTR = 2e-3,
  r_t_LTR = cfg_tune.r_i_LTR + 1e-3,
  r_o_LTR = cfg_tune.r_t_LTR + 1e-3,
  N_ch_LTR = 200,
  L_LTR = 1);    
 
  // select the configuration of parameters
  parameter Model.PBConfiguration cfg = cfg_default;
  
  // set the values of parameters accordingly
  parameter HEBoundaryCondition bc_HTR = cfg.bc_HTR;  
  parameter HEBoundaryCondition bc_LTR = cfg.bc_LTR;
  parameter HEBoundaryCondition bc_heater = cfg.bc_heater;
  
  parameter ThermoState st_bypass = cfg.st_bypass;
  
  parameter EntityGeoParam geo_HTR_hot = cfg.cfg_PCHE_HTR_hot.geo;
  parameter EntityGeoParam geo_HTR_cold = cfg.cfg_PCHE_HTR_cold.geo;
  parameter EntityGeoParam geo_HTR_tube = cfg.cfg_PCHE_HTR_tube.geo;
  
  parameter EntityGeoParam geo_LTR_hot = cfg.cfg_PCHE_LTR_hot.geo;
  parameter EntityGeoParam geo_LTR_cold = cfg.cfg_PCHE_LTR_cold.geo;
  parameter EntityGeoParam geo_LTR_tube = cfg.cfg_PCHE_LTR_tube.geo;  
  
  parameter EntityThermoParam thermo_HTR_hot = cfg.cfg_PCHE_HTR_hot.thermo;
  parameter EntityThermoParam thermo_HTR_cold = cfg.cfg_PCHE_HTR_cold.thermo;
  parameter EntityThermoParam thermo_HTR_tube = cfg.cfg_PCHE_HTR_tube.thermo;  
  
  parameter EntityThermoParam thermo_LTR_hot = cfg.cfg_PCHE_LTR_hot.thermo;
  parameter EntityThermoParam thermo_LTR_cold = cfg.cfg_PCHE_LTR_cold.thermo;
  parameter EntityThermoParam thermo_LTR_tube = cfg.cfg_PCHE_LTR_tube.thermo;  

  // use HTR's geo parameters as default 
  parameter EntityGeoParam geo_heater_hot = cfg.cfg_HTR_hot.geo;
  parameter EntityGeoParam geo_heater_cold = cfg.cfg_HTR_cold.geo;
  parameter EntityGeoParam geo_heater_tube = cfg.cfg_HTR_tube.geo;

  parameter EntityThermoParam thermo_heater_tube = cfg.cfg_heater_tube.thermo;

  parameter EntityThermoParam thermo_mixer = cfg.cfg_mixer.thermo;
  //Components
  inner ThermoPower.System system(allowFlowReversal = false, initOpt=ThermoPower.Choices.Init.Options.noInit) annotation(
    Placement(transformation(extent = {{80, 80}, {100, 100}})));
  
  // global init opition (system.initOpt) leads to order reduction error
  // use this flag to control the initialization of all components instead. 
  parameter Boolean SSInit = false "Steady-state initialization";

  ThermoPower.Water.SourceMassFlow source_heater_hot(
    redeclare package Medium = medium_heater,
    w0 = bc_heater.st_hot_in.mdot,
    p0 = bc_heater.st_hot_in.p,
    h = bc_heater.st_hot_in.h,
    T = bc_heater.st_hot_in.T) 
    annotation(
    Placement(transformation(origin = {0, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 270)));
    
  ThermoPower.Water.SinkPressure sink_heater_hot(
    redeclare package Medium = medium_heater, 
    p0 = bc_heater.st_hot_out.mdot, 
    T = bc_heater.st_hot_out.T,
    use_T = true) 
    annotation(
    Placement(transformation(origin = {0, -80}, extent = {{-10, -10}, {10, 10}}, rotation = 270)));

  ThermoPower.PowerPlants.HRSG.Components.HE heater(
    redeclare package FluidMedium = medium_heater, 
    redeclare package FlueGasMedium = medium_cold, 
    fluidFlow(fixedMassFlowSimplified = true, hstartin = bc_heater.st_hot_in.h, hstartout=bc_heater.st_hot_out.h), // set the fluid flow as fixed mdot for simplarity
    gasFlow(QuasiStatic = true, Tstartin = bc_heater.st_cold_in.T, Tstartout = bc_heater.st_cold_out.T),
    N_G=geo_heater_cold.N_seg,
    N_F=geo_heater_hot.N_seg,
    Nw_G=geo_heater_tube.N_seg,
    gasNomFlowRate=bc_heater.st_cold_in.mdot,
    gasNomPressure=bc_heater.st_cold_in.p,
    fluidNomFlowRate=bc_heater.st_hot_in.mdot,
    fluidNomPressure=bc_heater.st_hot_in.p,
    exchSurface_G=geo_heater_cold.A_ex,
    exchSurface_F=geo_heater_hot.A_ex,
    extSurfaceTub=geo_heater_tube.A_ex,
    gasVol=geo_heater_cold.V,
    fluidVol=geo_heater_hot.V,
    metalVol=geo_heater_tube.V,
    SSInit=SSInit,
    rhomcm=thermo_heater_tube.rho_mcm,
    lambda=thermo_heater_tube.lambda,
    Tstartbar_G=bc_heater.st_cold_in.T,
    Tstartbar_M=bc_heater.st_hot_in.T - 50,
    pstart_F = bc_heater.st_hot_in.p, 
    FluidPhaseStart=ThermoPower.Choices.FluidPhase.FluidPhases.Liquid,    
    redeclare replaceable model HeatTransfer_F =  ThermoPower.Thermal.HeatTransferFV.IdealHeatTransfer, // ConstantHeatTransferCoefficientTwoGrids(gamma = thermo_HTR_hot.gamma_he),    
    redeclare replaceable model HeatTransfer_G = ThermoPower.Thermal.HeatTransferFV.IdealHeatTransfer, // ConstantHeatTransferCoefficient(gamma = thermo_HTR_cold.gamma_he),     
    redeclare model HeatExchangerTopology = ThermoPower.Thermal.HeatExchangerTopologies.CounterCurrentFlow,
    metalTube(WallRes=false, Tstartbar=heater.Tstartbar_M)) annotation(
    Placement(transformation(extent = {{-20, -20}, {20, 20}}, rotation = 0)));

  ThermoPower.Gas.SourceMassFlow source_cold(
    redeclare package Medium = medium_cold, 
    T = bc_LTR.st_cold_in.T, 
    p0 = bc_LTR.st_cold_in.p, 
    use_in_T = false, 
    w0 = bc_LTR.st_cold_in.mdot) 
  annotation(
    Placement(transformation(origin = {0, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 270)));

  ThermoPower.Gas.SourceMassFlow source_mixer_in(
    redeclare package Medium = medium_cold,
    T = st_bypass.T,
    p0 = st_bypass.p,
    use_in_T = false,
    w0 = st_bypass.mdot    
  );
/*
  ThermoPower.Gas.SinkPressure sink_cold(
    redeclare package Medium = medium_cold, 
    // p0 = bc_HTR.st_cold_out.p, 
    p0 = bc_heater.st_cold_out.p,
    // T = bc_HTR.st_cold_out.T)
    T = bc_heater.st_cold_out.T) 
  annotation(
    Placement(transformation(origin = {0, -80}, extent = {{-10, -10}, {10, 10}}, rotation = 270)));
*/
  ThermoPower.Gas.SinkPressure sink_hot(
    redeclare package Medium = medium_hot,
    T = bc_LTR.st_hot_out.T, 
    p0 = bc_LTR.st_hot_out.p) 
  annotation(
    Placement(transformation(extent = {{60, -10}, {80, 10}}, rotation = 0)));
/*
  ThermoPower.Gas.SourceMassFlow source_hot(
    redeclare package Medium = medium_hot, 
    // T = from_degC(730.32), 
    //T = bc_heater.st_cold_out.T, 
    T = bc_HTR.st_hot_in.T, 
    //p0 = bc_heater.st_cold_out.p, 
    p0 = bc_HTR.st_hot_in.p, 
    //w0 = bc_heater.st_cold_out.mdot,
    w0 = bc_HTR.st_hot_in.mdot,
    use_in_T = false) 
  annotation(
    Placement(transformation(extent = {{-70, -10}, {-50, 10}}, rotation = 0))); 
*/
  // use FlowJoin to mix flow
  Gas.FlowJoin mixer(redeclare package Medium = medium_cold);  
 /*
  ThermoPower.Gas.SensT T_waterOut(
    redeclare package Medium = medium_cold) 
  annotation(
    Placement(transformation(origin = {4, -50}, extent = {{-10, -10}, {10, 10}}, rotation = 270)));

  ThermoPower.Gas.SensT T_gasOut(
    redeclare package Medium = medium_hot) 
  annotation(
    Placement(transformation(extent = {{30, -6}, {50, 14}}, rotation = 0)));
*/
  Components.TP_PCHE HTR(
  redeclare package FluidMedium = medium_cold, 
  redeclare package FlueGasMedium = medium_hot, 
  gasFlow(QuasiStatic = true, Tstartin = bc_HTR.st_hot_in.T, Tstartout = bc_HTR.st_hot_out.T),
  fluidFlow(Tstartin = bc_HTR.st_cold_in.T, Tstartout = bc_HTR.st_cold_out.T),     
  redeclare replaceable model HeatTransfer_F = Steps.Components.KimPCHEHeatTransferFV(
  pitch = pitch,
  phi = phi), 
  redeclare replaceable model HeatTransfer_G = Steps.Components.KimPCHEHeatTransferFV(
  pitch = pitch,
  phi = phi), 
  redeclare model HeatExchangerTopology = ThermoPower.Thermal.HeatExchangerTopologies.CounterCurrentFlow,  
  N_F = geo_HTR_cold.N_seg, 
  N_G = geo_HTR_hot.N_seg,    
  Nt = geo_HTR_hot.N_ch,   
  Tstartbar_G = bc_HTR.st_hot_in.T, 
  Tstartbar_F = bc_HTR.st_cold_in.T, 
  exchSurface_F = geo_HTR_cold.A_ex, 
  exchSurface_G = geo_HTR_hot.A_ex, 
  extSurfaceTub = geo_HTR_tube.A_ex, 
  fluidNomFlowRate = bc_HTR.st_cold_in.mdot, 
  fluidNomPressure = bc_HTR.st_cold_in.p, 
  fluidVol = geo_HTR_cold.V, 
  gasNomFlowRate = bc_HTR.st_hot_in.mdot, 
  gasNomPressure = bc_HTR.st_hot_in.p, 
  gasVol = geo_HTR_hot.V, 
  lambda = thermo_HTR_tube.lambda, 
  metalVol = geo_HTR_tube.V, 
  pstart_F = bc_HTR.st_cold_in.p, 
  pstart_G = bc_HTR.st_hot_in.p,
  rhomcm = thermo_HTR_tube.rho_mcm,
  SSInit=SSInit,
  gasQuasiStatic = true,
  fluidQuasiStatic = true,  
  metalTube(WallRes=false)) annotation(
    Placement(transformation(extent = {{-20, -20}, {20, 20}}, rotation = 0)));

  //Real kim_cor_coe[4] = {kim.a, kim.b, kim.c, kim.d};
  parameter SI.Length pitch = 12.3e-3 "pitch length";
  parameter Real phi = 35 "pitch angle °";
  
  Components.TP_PCHE LTR(
  redeclare package FluidMedium = medium_cold, 
  redeclare package FlueGasMedium = medium_hot, 
  redeclare replaceable model HeatTransfer_F = Steps.Components.KimPCHEHeatTransferFV(
  pitch = pitch,
  phi = phi),
  // ThermoPower.Thermal.HeatTransferFV.ConstantHeatTransferCoefficient(gamma = thermo_LTR_cold.gamma_he),
  // redeclare replaceable model HeatTransfer_G = ThermoPower.Thermal.HeatTransferFV.IdealHeatTransfer, 
  redeclare replaceable model HeatTransfer_G = Steps.Components.KimPCHEHeatTransferFV(
  pitch = pitch,
  phi = phi),  
  gasFlow(QuasiStatic = true, Tstartin = bc_LTR.st_hot_in.T, Tstartout = bc_LTR.st_hot_out.T),
  fluidFlow(Tstartin = bc_LTR.st_cold_in.T, Tstartout = bc_LTR.st_cold_out.T),   
  // redeclare replaceable model HeatTransfer_F = ThermoPower.Thermal.HeatTransferFV.IdealHeatTransfer,   
  // ThermoPower.Thermal.HeatTransferFV.ConstantHeatTransferCoefficient(gamma = thermo_LTR_cold.gamma_he),
  // redeclare replaceable model HeatTransfer_G = ThermoPower.Thermal.HeatTransferFV.IdealHeatTransfer,  
  // ThermoPower.Thermal.HeatTransferFV.ConstantHeatTransferCoefficientTwoGrids(gamma = thermo_LTR_hot.gamma_he),
  redeclare model HeatExchangerTopology = ThermoPower.Thermal.HeatExchangerTopologies.CounterCurrentFlow,  
  N_F = geo_LTR_cold.N_seg, 
  N_G = geo_LTR_hot.N_seg,   
  Nt = geo_LTR_hot.N_ch,  
  Tstartbar_G = bc_LTR.st_hot_in.T, 
  Tstartbar_F = bc_LTR.st_cold_in.T, 
  exchSurface_F = geo_LTR_cold.A_ex, 
  exchSurface_G = geo_LTR_hot.A_ex, 
  extSurfaceTub = geo_LTR_tube.A_ex, 
  fluidNomFlowRate = bc_LTR.st_cold_in.mdot, 
  fluidNomPressure = bc_LTR.st_cold_in.p, 
  fluidVol = geo_LTR_cold.V, 
  gasNomFlowRate = bc_LTR.st_hot_in.mdot, 
  gasNomPressure = bc_LTR.st_hot_in.p, 
  gasVol = geo_LTR_hot.V, 
  lambda = thermo_LTR_tube.lambda, 
  metalVol = geo_LTR_tube.V, 
  pstart_F = bc_LTR.st_cold_in.p, 
  pstart_G = bc_LTR.st_hot_in.p,
  rhomcm = thermo_LTR_tube.rho_mcm,
  SSInit=SSInit,
  gasQuasiStatic = true,
  fluidQuasiStatic = true,  
  metalTube(WallRes=false)) annotation(
    Placement(transformation(extent = {{-20, -20}, {20, 20}}, rotation = 0)));

  ThermoPower.Gas.Turbine Turbine1(
  redeclare package Medium = medium_hot, 
  fileName = Modelica.Utilities.Files.loadResource("modelica://Steps/Resources/Data/turbine_map.txt"),   
  tablePhic = fill(0.0, 14, 12), //tablePhic, 
  tableEta = fill(0.0, 14, 12), //tableEta, 
  pstart_in = bc_heater.st_cold_out.p, 
  pstart_out = bc_HTR.st_hot_in.p, 
  Tstart_in = bc_heater.st_cold_out.T, 
  Tstart_out = bc_HTR.st_hot_in.T, 
  Ndesign = 60000.0, 
  Tdes_in = bc_heater.st_cold_out.T,  
  Table = ThermoPower.Choices.TurboMachinery.TableTypes.file,
  explicitIsentropicEnthalpy = false,
  gas_in(
    p(nominal = Turbine1.pstart_in), 
    T(nominal = Turbine1.Tstart_in)),
  gas_iso(
    p(nominal = Turbine1.pstart_out), 
    T(nominal = Turbine1.Tstart_out))) 
  annotation(
    Placement(transformation(extent = {{-40, -20}, {0, 20}}, rotation = 0)));

  Modelica.Mechanics.Rotational.Sources.Speed speed1 annotation(
    Placement(visible = true, transformation(origin = {84, 0}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.Constant const1(k = 60000) annotation(
    Placement(visible = true, transformation(origin = {130, 0}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));  
 
  ThermoPower.Gas.SensT sens_turbine(redeclare package Medium = medium_hot) annotation(
    Placement(visible = true, transformation(origin = {20, 44}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));

protected
  parameter Real tablePhic[5, 4] = [1, 37, 80, 100; 1.5, 7.10E-05, 7.10E-05, 7.10E-05; 2, 8.40E-05, 8.40E-05, 8.40E-05; 2.5, 8.70E-05, 8.70E-05, 8.70E-05; 3, 1.04E-04, 1.04E-04, 1.04E-04];
  parameter Real tableEta[5, 4] = [1, 37, 80, 100; 1.5, 0.57, 0.89, 0.81; 2, 0.46, 0.82, 0.88; 2.5, 0.41, 0.76, 0.85; 3, 0.38, 0.72, 0.82];
  
equation

  //HTR + mixer + LTR + Heater + Turbine - OPEN LOOP
  // main stream, water/cold side  
  connect(source_mixer_in.flange, mixer.inlet1);  
  connect(source_cold.flange, LTR.waterIn);  
  connect(LTR.waterOut, mixer.inlet2);  
  connect(mixer.outlet, HTR.waterIn);
  connect(HTR.waterOut, heater.gasIn);    
  connect(heater.gasOut, Turbine1.inlet);

  // connect(heater.gasOut, sink_cold.flange); 
  // connect(source_hot.flange, Turbine1.inlet);
  
  connect(Turbine1.outlet, sens_turbine.inlet) annotation(
   Line(points = {{-50, 0}, {-20, 0}}, color = {159, 159, 223}, thickness = 0.5, smooth = Smooth.None));

  connect(Turbine1.shaft_b, speed1.flange) annotation(
    Line(points = {{30, 0}, {74, 0}, {74, 0}, {74, 0}}));
  connect(speed1.w_ref, const1.y) annotation(
    Line(points = {{96, 0}, {120, 0}, {120, 0}, {118, 0}}, color = {0, 0, 127}));
  
  // connect(sens_turbine.outlet, sink_cold.flange);  
  // connect(source_hot.flange, HTR.gasIn);  
  
  connect(sens_turbine.outlet, HTR.gasIn);   

  connect(HTR.gasOut, LTR.gasIn);
  connect(LTR.gasOut, sink_hot.flange) annotation(
    Line(points = {{46, 0}, {46, 0}, {60, 0}}, color = {159, 159, 223}, thickness = 0.5));
   
  // hot stream for heater
  connect(source_heater_hot.flange, heater.waterIn);
  connect(heater.waterOut, sink_heater_hot.flange);

/*
  //HTR + mixer + LTR + Heater
  // main stream, water/cold side  
  connect(source_mixer_in.flange, mixer.inlet1);  
  connect(source_cold.flange, LTR.waterIn);  
  connect(LTR.waterOut, mixer.inlet2);  
  connect(mixer.outlet, HTR.waterIn);
  connect(HTR.waterOut, heater.gasIn);  
  connect(heater.gasOut, sink_cold.flange) annotation(
    Line(points = {{1.83697e-015, -70}, {1.83697e-015, -56}, {-8.88178e-016, -56}}, thickness = 0.5, color = {0, 0, 255}));
  
  // main stream, gas/hot side
  connect(source_hot.flange, HTR.gasIn) annotation(
   Line(points = {{-50, 0}, {-20, 0}}, color = {159, 159, 223}, thickness = 0.5, smooth = Smooth.None));
  connect(HTR.gasOut, LTR.gasIn);  
  connect(LTR.gasOut, sink_hot.flange) annotation(
    Line(points = {{46, 0}, {46, 0}, {60, 0}}, color = {159, 159, 223}, thickness = 0.5));
   
  // hot stream for heater
  connect(source_heater_hot.flange, heater.waterIn);
  connect(heater.waterOut, sink_heater_hot.flange);
*/

/*
  // mixer alone
  connect(source_mixer_in.flange, mixer.in1);
  connect(source_cold.flange, mixer.in2);
  connect(mixer.out, sink_cold.flange);
*/  

/*
  // HTR alone
  connect(source_cold.flange, HTR.waterIn) annotation(
    Line(points = {{-1.83697e-015, 50}, {-1.83697e-015, 20}, {0, 20}}, color = {0, 0, 255}, thickness = 0.5, smooth = Smooth.None));
  
  connect(T_waterOut.inlet, HTR.waterOut) annotation(
    Line(points = {{8.88178e-016, -44}, {8.88178e-016, -20}, {0, -20}}, thickness = 0.5, color = {0, 0, 255}));      
   
  connect(sink_cold.flange, T_waterOut.outlet) annotation(
    Line(points = {{1.83697e-015, -70}, {1.83697e-015, -56}, {-8.88178e-016, -56}}, thickness = 0.5, color = {0, 0, 255}));
  
  connect(source_hot.flange, HTR.gasIn) annotation(
        Line(points = {{-50, 0}, {-20, 0}}, color = {159, 159, 223}, thickness = 0.5, smooth = Smooth.None));
    
  connect(T_gasOut.inlet, HTR.gasOut) annotation(
    Line(points = {{34, 0}, {34, 0}, {20, 0}}, color = {159, 159, 223}, thickness = 0.5));
    
  connect(T_gasOut.outlet, sink_hot.flange) annotation(
    Line(points = {{46, 0}, {46, 0}, {60, 0}}, color = {159, 159, 223}, thickness = 0.5));    

  // LTR alone
  // water/cold side  
  connect(source_cold.flange, LTR.waterIn);
  
  connect(LTR.waterOut, sink_cold.flange) annotation(
    Line(points = {{1.83697e-015, -70}, {1.83697e-015, -56}, {-8.88178e-016, -56}}, thickness = 0.5, color = {0, 0, 255}));
  
  // gas/hot side
  connect(source_hot.flange, LTR.gasIn);
  
  connect(LTR.gasOut, sink_hot.flange) annotation(
    Line(points = {{46, 0}, {46, 0}, {60, 0}}, color = {159, 159, 223}, thickness = 0.5));
  */  
 
/*
  //HTR + mixer

  connect(source_mixer_in.flange, mixer.in1);  
  connect(source_cold.flange, mixer.in2);  
  connect(mixer.out, HTR.waterIn);

  connect(source_cold.flange, HTR.waterIn);
  connect(HTR.waterOut, T_waterOut.inlet) annotation(
    Line(points = {{8.88178e-016, -44}, {8.88178e-016, -20}, {0, -20}}, thickness = 0.5, color = {0, 0, 255}));
  connect(sink_cold.flange, T_waterOut.outlet) annotation(
    Line(points = {{1.83697e-015, -70}, {1.83697e-015, -56}, {-8.88178e-016, -56}}, thickness = 0.5, color = {0, 0, 255}));
    
  connect(source_hot.flange, HTR.gasIn) annotation(
   Line(points = {{-50, 0}, {-20, 0}}, color = {159, 159, 223}, thickness = 0.5, smooth = Smooth.None));   
  connect(HTR.gasOut, T_gasOut.inlet ) annotation(
    Line(points = {{34, 0}, {34, 0}, {20, 0}}, color = {159, 159, 223}, thickness = 0.5));
  connect(T_gasOut.outlet, sink_hot.flange) annotation(
    Line(points = {{46, 0}, {46, 0}, {60, 0}}, color = {159, 159, 223}, thickness = 0.5));
*/
/*
  // mixer + LTR
  // water/cold side  
  connect(source_mixer_in.flange, mixer.inlet1);
  
  connect(source_cold.flange, LTR.waterIn);
  
  connect(LTR.waterOut, mixer.inlet2);
  
  connect(mixer.outlet, sink_cold.flange) annotation(
    Line(points = {{1.83697e-015, -70}, {1.83697e-015, -56}, {-8.88178e-016, -56}}, thickness = 0.5, color = {0, 0, 255}));
  
  // gas/hot side
  connect(source_hot.flange, LTR.gasIn);
  
  connect(LTR.gasOut, sink_hot.flange) annotation(
    Line(points = {{46, 0}, {46, 0}, {60, 0}}, color = {159, 159, 223}, thickness = 0.5));

*/
/*
  //HTR + mixer + LTR 
  // water/cold side  
  connect(source_mixer_in.flange, mixer.inlet1);
  
  connect(source_cold.flange, LTR.waterIn);
  
  connect(LTR.waterOut, mixer.inlet2);
  
  connect(mixer.outlet, HTR.waterIn);

  connect(HTR.waterOut, sink_cold.flange) annotation(
    Line(points = {{1.83697e-015, -70}, {1.83697e-015, -56}, {-8.88178e-016, -56}}, thickness = 0.5, color = {0, 0, 255}));
  
  // gas/hot side
  connect(source_hot.flange, HTR.gasIn) annotation(
   Line(points = {{-50, 0}, {-20, 0}}, color = {159, 159, 223}, thickness = 0.5, smooth = Smooth.None));

  connect(HTR.gasOut, LTR.gasIn);
  
  connect(LTR.gasOut, sink_hot.flange) annotation(
    Line(points = {{46, 0}, {46, 0}, {60, 0}}, color = {159, 159, 223}, thickness = 0.5));
*/


annotation(
    Diagram(graphics),
    experiment(StartTime = 0, StopTime = 300, Tolerance = 1e-3, Interval = 10),
    __OpenModelica_commandLineOptions = "--matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection -d=initialization,NLSanalyticJacobian,aliasConflicts,bltdump",    
    __OpenModelica_simulationFlags(lv = "LOG_DEBUG,LOG_NLS,LOG_NLS_V,LOG_STATS,LOG_INIT,LOG_STDOUT, -w", outputFormat = "mat", s = "dassl", nls = "homotopy"));
end TestTP_Components_PCHE;
