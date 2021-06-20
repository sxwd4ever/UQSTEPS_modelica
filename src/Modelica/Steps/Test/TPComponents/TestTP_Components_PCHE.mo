within Steps.Test.TPComponents;

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
  
  //package medium_main = Modelica.Media.IdealGases.SingleGases.CO2; //Steps.Media.CO2;
  //package medium_main = Modelica.Media.IdealGases.SingleGases.CO2; //Steps.Media.CO2;
  // package medium_main = Steps.Media.CO2;
  // package medium_main = Steps.Media.CO2;
  package medium_main = Steps.Media.SCO2(
    // inputChoice = ExternalMedia.Common.InputChoice.pT
  );
  package medium_heater = SolarTherm.Media.Sodium.Sodium_pT;
  //package medium_heater = ThermoPower.Water.StandardWater;// Modelica.Media.IdealGases.SingleGases.CO2;

  parameter Real table_k_metalwall[:,:] = [293.15, 12.1; 373.15, 16.3; 773.15, 21.5]; 
  parameter Real Cf_C1 = 1.626, Cf_C2 = 1, Cf_C3 = 1;
  // parameter Real Cf_C1_cold = 1, Cf_C2_cold = 1, Cf_C3_cold = 1;
  parameter Real use_rho_bar = -1.0;  
  parameter Real rho_bar_hot = 1.0;
  parameter Real rho_bar_cold = 1.0;  
  
  // select the configuration of parameters
  parameter Model.RCBCycleConfig cfg(
    redeclare package medium_heater = medium_heater,
    redeclare package medium_main   = medium_main,
    // mdot_heater      = 40,
    // T_heater_hot_in  = from_degC(800),
    // T_heater_hot_out = from_degC(600),
    r_i_heater  = 20e-3,
    r_t_heater  = cfg.r_i_heater + 10e-3,
    r_o_heater  = 1/2,                      // agree with the final parameter Dhyd = 1 in HE, should be checked to see if it is capable of containing all fluid-metal tubes
    N_ch_heater = 100,
    L_heater    = 1,
    N_ch_HTR    = 30000,
    L_HTR       = 2.5,
    r_i_HTR     = 1.5e-3,
    r_o_HTR     = 1.5e-3,
    N_ch_LTR    = 30000,
    L_LTR       = 2.5,
    r_i_LTR     = 1.5e-3,
    r_o_LTR     = 1.5e-3,
    Ns_turb     = 30000,
    // mdot_main = 100,
    mdot_main   = 128.774,   // this value is a simulation result with sourcePressure.p = 200 bar
    mdot_heater = 40,
    table_k_LTR_wall = table_k_metalwall,
    table_k_HTR_wall = table_k_metalwall
  );
  
  parameter Model.HeatExchangerConfig cfg_heater = cfg.cfg_heater;
  parameter Model.HeatExchangerConfig cfg_LTR    = cfg.cfg_LTR;
  parameter Model.HeatExchangerConfig cfg_HTR    = cfg.cfg_HTR;
  parameter Model.TurbomachineryConfig cfg_turb  = cfg.cfg_turb;
  parameter Model.SplitterConfig cfg_merger      = cfg.cfg_merger;
  
  parameter Model.ThermoState st_bypass      = cfg.st_recomp_out;
  parameter Model.ThermoState st_source_hot  = cfg_HTR.cfg_hot.st_in;
  parameter Model.ThermoState st_sink_hot    = cfg_LTR.cfg_hot.st_out;
  parameter Model.ThermoState st_source_cold = cfg_LTR.cfg_cold.st_in;
  parameter Model.ThermoState st_sink_cold   = cfg_heater.cfg_cold.st_out;

  parameter Integer N_seg_heater = cfg.cfg_heater.cfg_hot.geo_path.N_seg; 
  parameter Integer N_seg_LTR    = cfg.cfg_LTR.cfg_hot.geo_path.N_seg; 
  parameter Integer N_seg_HTR    = cfg.cfg_HTR.cfg_hot.geo_path.N_seg; 

  //Components
  // global init opition (system.initOpt) leads to order reduction error
  // use this flag to control the initialization of all components instead. 
  parameter Boolean SSInit = false "Steady-state initialization";  
  inner ThermoPower.System system(allowFlowReversal = false, initOpt=ThermoPower.Choices.Init.Options.noInit) annotation(
    Placement(transformation(extent = {{80, 80}, {100, 100}})));

  ThermoPower.Water.SourceMassFlow source_heater_hot(
    redeclare package Medium = medium_heater,
    w0 = cfg_heater.cfg_hot.st_in.mdot,
    p0 = cfg_heater.cfg_hot.st_in.p,
    h  = cfg_heater.cfg_hot.st_in.h,
    T  = cfg_heater.cfg_hot.st_in.T) 
    annotation(
    Placement(transformation(origin = {0, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 270)));
    
  ThermoPower.Water.SinkPressure sink_heater_hot(
    redeclare package Medium = medium_heater, 
    p0    = cfg_heater.cfg_hot.st_out.p,
    T     = cfg_heater.cfg_hot.st_out.T,
    use_T = false) 
    annotation(
    Placement(transformation(origin = {0, -80}, extent = {{-10, -10}, {10, 10}}, rotation = 270)));

  Steps.TPComponents.HE heater(
    redeclare package FluidMedium = medium_heater, 
    redeclare package FlueGasMedium = medium_main, 
    
    // // DittusBoelter
    // redeclare replaceable model HeatTransfer_F = ThermoPower.Thermal.HeatTransferFV.DittusBoelter,    
    // redeclare replaceable model HeatTransfer_G = ThermoPower.Thermal.HeatTransferFV.DittusBoelter,
    // fluidFlow(
    //   heatTransfer(heating=true), 
    //   fixedMassFlowSimplified = true,
    //   hstartin                = cfg_heater.cfg_hot.st_in.h,
    //   hstartout               = cfg_heater.cfg_hot.st_out.h),   // set the fluid flow as fixed mdot for simplarity
    // gasFlow(
    //   heatTransfer(heating=false), 
    //   Tstartin  = cfg_heater.cfg_cold.st_in.T,
    //   Tstartout = cfg_heater.cfg_cold.st_out.T),
    
    
    redeclare replaceable model HeatTransfer_F = ThermoPower.Thermal.HeatTransferFV.IdealHeatTransfer,    
    redeclare replaceable model HeatTransfer_G = ThermoPower.Thermal.HeatTransferFV.IdealHeatTransfer, 
    fluidFlow(
      fixedMassFlowSimplified = true,
      hstartin                = cfg_heater.cfg_hot.st_in.h,
      hstartout               = cfg_heater.cfg_hot.st_out.h),   // set the fluid flow as fixed mdot for simplarity
    gasFlow(
      Tstartin    = cfg_heater.cfg_cold.st_in.T,
      Tstartout   = cfg_heater.cfg_cold.st_out.T),
    
    // redeclare replaceable model HeatTransfer_F = ThermoPower.Thermal.HeatTransferFV.ConstantHeatTransferCoefficientTwoGrids(gamma = thermo_hot.gamma_he),          
    // redeclare replaceable model HeatTransfer_G = ThermoPower.Thermal.HeatTransferFV.ConstantHeatTransferCoefficient(gamma =  thermo_cold.gamma_he),           
    redeclare model HeatExchangerTopology = ThermoPower.Thermal.HeatExchangerTopologies.CounterCurrentFlow,
    cfg             = cfg_heater,
    N_G             = N_seg_heater,
    N_F             = N_seg_heater,
    SSInit          = false,
    gasQuasiStatic  = true,
    FluidPhaseStart = ThermoPower.Choices.FluidPhase.FluidPhases.Liquid,
    metalTube(WallRes=false))
    annotation(
      Placement(transformation(extent = {{-20, -20}, {20, 20}}, rotation = 0)));
  
  Steps.TPComponents.WaterStateReader r_heater_hin(redeclare package Medium = medium_heater) annotation(
    Placement(visible = true, transformation(origin = {66, -56}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
  Steps.TPComponents.WaterStateReader r_heater_hout(redeclare package Medium = medium_heater) annotation(
    Placement(visible = true, transformation(origin = {96, 44}, extent = {{-6, -6}, {6, 6}}, rotation = 180)));
  Steps.TPComponents.GasStateReader r_heater_cin(redeclare package Medium = medium_main) annotation(
    Placement(visible = true, transformation(origin = {-22, 50}, extent = {{-6, -6}, {6, 6}}, rotation = -90)));
  Steps.TPComponents.GasStateReader r_heater_cout(redeclare package Medium = medium_main) annotation(
    Placement(visible = true, transformation(origin = {-32, 20}, extent = {{-6, -6}, {6, 6}}, rotation = 180)));

  
  // for HTR/LRT 's cold side, a gas source/sink
  // for PCHE's cold side
  /*
  ThermoPower.Gas.SourceMassFlow source_cold(
    redeclare package Medium = medium_main, 
    T        = st_source_cold.T,
    p0       = st_source_cold.p,
    use_in_T = false,
    w0       = st_source_cold.mdot,
    gas(
      p(nominal = st_source_cold.p),
      T(nominal = st_source_cold.T)
    )) 
  annotation(
    Placement(transformation(origin = {0, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 270)));
  */
  
  ThermoPower.Gas.SourcePressure source_cold(
    redeclare package Medium = medium_main, 
    T        = st_source_cold.T,
    p0       = st_source_cold.p,
    use_in_T = false,
    gas(
      p(start = st_source_cold.p, nominal = st_source_cold.p), 
      T(start = st_source_cold.T, nominal = st_source_cold.T))) 
  annotation(
    Placement(transformation(origin = {0, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 270)));

/*
  ThermoPower.Gas.SinkPressure sink_cold(
    redeclare package Medium = medium_main,     
    // for HTR + Mixer + LTR  
    // p0 = bc_HTR.st_cold_out.p,
    // T = bc_HTR.st_cold_out.T)
    // for HTR + Mixer + LTR + Heater
    p0 = st_sink_cold.p,
    T  = st_sink_cold.T
  )
  annotation(
    Placement(transformation(origin = {0, -80}, extent = {{-10, -10}, {10, 10}}, rotation = 270)));
  

  ThermoPower.Gas.SourceMassFlow source_hot(
    redeclare package Medium = medium_main, 
    // T = from_degC(730.32), 
    //T = bc_heater.st_cold_out.T, 
    //p0 = bc_heater.st_cold_out.p, 
    //w0 = bc_heater.st_cold_out.mdot,
    T        = st_source_hot.T,
    p0       = st_source_hot.p,
    w0       = st_source_hot.mdot,
    use_in_T = false,
    gas(
      p(nominal = st_source_hot.p),
      T(nominal = st_source_hot.T)
    ))
  annotation(
    Placement(transformation(extent = {{-70, -10}, {-50, 10}}, rotation = 0))); 
*/

  ThermoPower.Gas.SinkPressure sink_hot(
    redeclare package Medium = medium_main,
    T  = st_sink_hot.T,
    p0 = st_sink_hot.p) 
  annotation(
    Placement(transformation(extent = {{60, -10}, {80, 10}}, rotation = 0)));

  ThermoPower.Gas.SourceMassFlow source_mixer_in(
    redeclare package Medium = medium_main,
    T        = st_bypass.T,
    p0       = st_bypass.p,
    use_in_T = false,
    w0       = st_bypass.mdot
  );

  // use FlowJoin to mix flow
  Gas.FlowJoin mixer(redeclare package Medium = medium_main);  


  Steps.TPComponents.PCHE HTR(
    redeclare package FluidMedium   = medium_main,
    redeclare package FlueGasMedium = medium_main,
    // with Marchionni Correlation    
    redeclare replaceable model HeatTransfer_F = Steps.TPComponents.GnielinskiHeatTransferFV(),
    redeclare replaceable model HeatTransfer_G = Steps.TPComponents.GnielinskiHeatTransferFV(),
    gasFlow(
      heatTransfer(
        pitch                 = cfg_HTR.cfg_hot.l_pitch,
        phi                   = cfg_HTR.cfg_hot.a_phi,
        Cf_C1                 = Cf_C1,
        Cf_C2                 = Cf_C2,
        Cf_C3                 = Cf_C3,
        use_rho_bar           = use_rho_bar,
        rho_bar               = rho_bar_hot,
        useAverageTemperature = false)),
    fluidFlow(
      heatTransfer(
        pitch                 = cfg_HTR.cfg_cold.l_pitch,
        phi                   = cfg_HTR.cfg_cold.a_phi,
        Cf_C1                 = Cf_C1,
        Cf_C2                 = Cf_C2,
        Cf_C3                 = Cf_C3,
        use_rho_bar           = use_rho_bar,
        rho_bar               = rho_bar_cold,
        useAverageTemperature = false)
    ),  
    redeclare model HeatExchangerTopology = ThermoPower.Thermal.HeatExchangerTopologies.CounterCurrentFlow,  
    cfg              = cfg_HTR,
    N_G              = N_seg_HTR,
    N_F              = N_seg_HTR,
    SSInit           = SSInit,
    gasQuasiStatic   = true,
    fluidQuasiStatic = true,
    table_k_metalwall =  table_k_metalwall,
    //metalQuasiStatic = false,   
    metalWall(WallRes=false)
  ) 
  annotation(
    Placement(transformation(extent = {{-20, -20}, {20, 20}}, rotation = 0)));
  
  Steps.TPComponents.GasStateReader r_HTR_hin(redeclare package Medium = medium_main) annotation(
    Placement(visible = true, transformation(origin = {66, -56}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
  Steps.TPComponents.GasStateReader r_HTR_hout(redeclare package Medium = medium_main) annotation(
    Placement(visible = true, transformation(origin = {96, 44}, extent = {{-6, -6}, {6, 6}}, rotation = 180)));
  Steps.TPComponents.GasStateReader r_HTR_cin(redeclare package Medium = medium_main) annotation(
    Placement(visible = true, transformation(origin = {-22, 50}, extent = {{-6, -6}, {6, 6}}, rotation = -90)));
  Steps.TPComponents.GasStateReader r_HTR_cout(redeclare package Medium = medium_main) annotation(
    Placement(visible = true, transformation(origin = {-32, 20}, extent = {{-6, -6}, {6, 6}}, rotation = 180)));

  Steps.TPComponents.PCHE LTR(
    redeclare package FluidMedium   = medium_main,
    redeclare package FlueGasMedium = medium_main,
    // with Marchionni Correlation    
    redeclare replaceable model HeatTransfer_F = Steps.TPComponents.GnielinskiHeatTransferFV(),
    redeclare replaceable model HeatTransfer_G = Steps.TPComponents.GnielinskiHeatTransferFV(),
    gasFlow(
      heatTransfer(
        pitch                 = cfg_LTR.cfg_hot.l_pitch,
        phi                   = cfg_LTR.cfg_hot.a_phi,
        Cf_C1                 = Cf_C1,
        Cf_C2                 = Cf_C2,
        Cf_C3                 = Cf_C3,
        use_rho_bar           = use_rho_bar,
        rho_bar               = rho_bar_hot,
        useAverageTemperature = false)),
    fluidFlow(
      heatTransfer(
        pitch                 = cfg_LTR.cfg_cold.l_pitch,
        phi                   = cfg_LTR.cfg_cold.a_phi,
        Cf_C1                 = Cf_C1,
        Cf_C2                 = Cf_C2,
        Cf_C3                 = Cf_C3,
        use_rho_bar           = use_rho_bar,
        rho_bar               = rho_bar_cold,
        useAverageTemperature = false)
    ),    
    cfg              = cfg_LTR,
    N_G              = N_seg_LTR,
    N_F              = N_seg_LTR,
    SSInit           = SSInit,
    gasQuasiStatic   = true,
    fluidQuasiStatic = true,
    table_k_metalwall =   cfg.table_k_LTR_wall,
    //metalQuasiStatic = false,   
    metalWall(WallRes=false))
    
  annotation(
    Placement(transformation(extent = {{-20, -20}, {20, 20}}, rotation = 0)));

  Steps.TPComponents.GasStateReader r_LTR_hin(redeclare package Medium = medium_main) annotation(
    Placement(visible = true, transformation(origin = {66, -56}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
  Steps.TPComponents.GasStateReader r_LTR_hout(redeclare package Medium = medium_main) annotation(
    Placement(visible = true, transformation(origin = {96, 44}, extent = {{-6, -6}, {6, 6}}, rotation = 180)));
  Steps.TPComponents.GasStateReader r_LTR_cin(redeclare package Medium = medium_main) annotation(
    Placement(visible = true, transformation(origin = {-22, 50}, extent = {{-6, -6}, {6, 6}}, rotation = -90)));
  Steps.TPComponents.GasStateReader r_LTR_cout(redeclare package Medium = medium_main) annotation(
    Placement(visible = true, transformation(origin = {-32, 20}, extent = {{-6, -6}, {6, 6}}, rotation = 180)));

  ThermoPower.Gas.Turbine Turbine1(
    redeclare package Medium = medium_main, 
    fileName   = Modelica.Utilities.Files.loadResource("modelica://Steps/Resources/Data/turbine_map.txt"),
    tablePhic  = fill(0.0, 14, 12),                                                                          //tablePhic, 
    tableEta   = fill(0.0, 14, 12),                                                                          //tableEta, 
    pstart_in  = cfg_turb.st_in.p,
    pstart_out = cfg_turb.st_out.p,
    Tstart_in  = cfg_turb.st_in.T,
    Tstart_out = cfg_turb.st_out.T,
    Ndesign    = cfg_turb.N,
    Tdes_in    = cfg_turb.st_in.T,
    Table      = ThermoPower.Choices.TurboMachinery.TableTypes.file,
    //explicitIsentropicEnthalpy = false,
    gas_in(
      p(nominal = Turbine1.pstart_in), 
      T(nominal = Turbine1.Tstart_in)),
    gas_iso(
      p(nominal = Turbine1.pstart_out), 
      T(nominal = Turbine1.Tstart_out),
      h(start   = cfg_turb.st_out.h)))    
    annotation(
      Placement(transformation(extent = {{-40, -20}, {0, 20}}, rotation = 0)));

  Modelica.Mechanics.Rotational.Sources.ConstantSpeed const_speed_comp(
      w_fixed=cfg_turb.N, useSupport=false) annotation(
    Placement(visible = true, transformation(origin = {81, -7}, extent = {{-5, -5}, {5, 5}}, rotation = 0)));  
/* 
  ThermoPower.Gas.SensT sens_turbine(redeclare package Medium = medium_main) annotation(
    Placement(visible = true, transformation(origin = {20, 44}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
*/
protected
  parameter Real tablePhic[5, 4] = [1, 37, 80, 100; 1.5, 7.10E-05, 7.10E-05, 7.10E-05; 2, 8.40E-05, 8.40E-05, 8.40E-05; 2.5, 8.70E-05, 8.70E-05, 8.70E-05; 3, 1.04E-04, 1.04E-04, 1.04E-04];
  parameter Real tableEta[5, 4] = [1, 37, 80, 100; 1.5, 0.57, 0.89, 0.81; 2, 0.46, 0.82, 0.88; 2.5, 0.41, 0.76, 0.85; 3, 0.38, 0.72, 0.82];
  
equation

  // HTR + mixer + LTR + Heater + Turbine - OPEN LOOP
  // main stream, water/cold side  
  connect(source_mixer_in.flange, mixer.inlet1);  
  connect(source_cold.flange, r_LTR_cin.inlet);
  connect(r_LTR_cin.outlet, LTR.waterIn);  
  connect(LTR.waterOut, r_LTR_cout.inlet);
  connect(r_LTR_cout.outlet, mixer.inlet2);  
  connect(mixer.outlet, r_HTR_cin.inlet);
  connect(r_HTR_cin.outlet, HTR.waterIn);
  connect(HTR.waterOut, r_HTR_cout.inlet);
  connect(r_HTR_cout.outlet, r_heater_cin.inlet);
  connect(r_heater_cin.outlet, heater.gasIn);    
  connect(heater.gasOut, r_heater_cout.inlet);
  connect(r_heater_cout.outlet, Turbine1.inlet);
  
  connect(Turbine1.outlet, sens_turbine.inlet) annotation(
   Line(points = {{-50, 0}, {-20, 0}}, color = {159, 159, 223}, thickness = 0.5, smooth = Smooth.None));   
  connect(sens_turbine.outlet, r_HTR_hin.inlet);
  connect(r_HTR_hin.outlet, HTR.gasIn);  
  connect(HTR.gasOut, r_HTR_hout.inlet);
  connect(r_HTR_hout.outlet, r_LTR_hin.inlet);
  connect(r_LTR_hin.outlet, LTR.gasIn);
  connect(LTR.gasOut, r_LTR_hout.inlet);
  connect(r_LTR_hout.outlet, sink_hot.flange) annotation(
    Line(points = {{46, 0}, {46, 0}, {60, 0}}, color = {159, 159, 223}, thickness = 0.5));
   
  connect(Turbine1.shaft_b, const_speed_comp.flange) annotation(
    Line(points = {{30, 0}, {74, 0}, {74, 0}, {74, 0}}));
  
  // connect(sens_turbine.outlet, sink_cold.flange);  
  // connect(source_hot.flange, HTR.gasIn);    

  // hot stream for heater
  connect(source_heater_hot.flange, r_heater_hin.inlet);
  connect(r_heater_hin.outlet, heater.waterIn);
  connect(heater.waterOut, r_heater_hout.inlet);
  connect(r_heater_hout.outlet, sink_heater_hot.flange);

/*
  //HTR + mixer + LTR + Heater
  // main stream, water/cold side  
  connect(source_mixer_in.flange, mixer.inlet1);  
  connect(source_cold.flange, r_LTR_cin.inlet);
  connect(r_LTR_cin.outlet, LTR.waterIn);  
  connect(LTR.waterOut, r_LTR_cout.inlet);
  connect(r_LTR_cout.outlet, mixer.inlet2);  
  connect(mixer.outlet, r_HTR_cin.inlet);
  connect(r_HTR_cin.outlet, HTR.waterIn);
  connect(HTR.waterOut, r_HTR_cout.inlet);
  connect(r_HTR_cout.outlet, r_heater_cin.inlet);
  connect(r_heater_cin.outlet, heater.gasIn);  
  connect(heater.gasOut, r_heater_cout.inlet);
  connect(r_heater_cout.outlet, sink_cold.flange) annotation(
    Line(points = {{1.83697e-015, -70}, {1.83697e-015, -56}, {-8.88178e-016, -56}}, thickness = 0.5, color = {0, 0, 255}));
  
  // main stream, gas/hot side
  connect(source_hot.flange, r_HTR_hin.inlet);
  connect(r_HTR_hin.outlet, HTR.gasIn) annotation(
   Line(points = {{-50, 0}, {-20, 0}}, color = {159, 159, 223}, thickness = 0.5, smooth = Smooth.None));
  connect(HTR.gasOut, r_HTR_hout.inlet);
  connect(r_HTR_hout.outlet, r_LTR_hin.inlet);
  connect(r_LTR_hin.outlet,LTR.gasIn);  
  connect(LTR.gasOut, r_LTR_hout.inlet);
  connect(r_LTR_hout.outlet, sink_hot.flange) annotation(
    Line(points = {{46, 0}, {46, 0}, {60, 0}}, color = {159, 159, 223}, thickness = 0.5));
   
  // hot stream for heater
  connect(source_heater_hot.flange, r_heater_hin.inlet);
  connect(r_heater_hin.outlet, heater.waterIn);
  connect(heater.waterOut, r_heater_hout.inlet);
  connect(r_heater_hout.outlet, sink_heater_hot.flange);
*/

/*
  // HTR + Heater + turbine
  
  connect(source_cold.flange, r_HTR_cin.inlet);
  connect(r_HTR_cin.outlet, HTR.waterIn);
  connect(HTR.waterOut, r_HTR_cout.inlet);
  connect(r_HTR_cout.outlet, r_heater_cin.inlet);
  connect(r_heater_cin.outlet, heater.gasIn);    
  connect(heater.gasOut, r_heater_cout.inlet);
  connect(r_heater_cout.outlet, Turbine1.inlet);
  
  connect(Turbine1.outlet, sens_turbine.inlet) annotation(
   Line(points = {{-50, 0}, {-20, 0}}, color = {159, 159, 223}, thickness = 0.5, smooth = Smooth.None));   
  connect(sens_turbine.outlet, r_HTR_hin.inlet);
  connect(r_HTR_hin.outlet, HTR.gasIn);  
  connect(HTR.gasOut, r_HTR_hout.inlet);
  connect(r_HTR_hout.outlet, sink_hot.flange) annotation(
    Line(points = {{46, 0}, {46, 0}, {60, 0}}, color = {159, 159, 223}, thickness = 0.5));
   
  connect(Turbine1.shaft_b, const_speed_comp.flange) annotation(
    Line(points = {{30, 0}, {74, 0}, {74, 0}, {74, 0}}));
  
  // connect(sens_turbine.outlet, sink_cold.flange);  
  // connect(source_hot.flange, HTR.gasIn);    

  // hot stream for heater
  connect(source_heater_hot.flange, r_heater_hin.inlet);
  connect(r_heater_hin.outlet, heater.waterIn);
  connect(heater.waterOut, r_heater_hout.inlet);
  connect(r_heater_hout.outlet, sink_heater_hot.flange);
*/

/*
  //HTR + mixer + LTR 
  // water/cold side  
  connect(source_mixer_in.flange, mixer.inlet1);  
  connect(source_cold.flange, r_LTR_cin.inlet);
  connect(r_LTR_cin.outlet, LTR.waterIn);  
  connect(LTR.waterOut, r_LTR_cout.inlet);
  connect(r_LTR_cout.outlet, mixer.inlet2);  
  connect(mixer.outlet, r_HTR_cin.inlet);
  connect(r_HTR_cin.outlet, HTR.waterIn);
  connect(HTR.waterOut, r_HTR_cout.inlet);
  connect(r_HTR_cout.outlet, sink_cold.flange) 
  annotation(
    Line(points = {{1.83697e-015, -70}, {1.83697e-015, -56}, {-8.88178e-016, -56}}, thickness = 0.5, color = {0, 0, 255}));
  
  // gas/hot side
  connect(source_hot.flange, r_HTR_hin.inlet);
  connect(r_HTR_hin.outlet, HTR.gasIn) 
  annotation(
   Line(points = {{-50, 0}, {-20, 0}}, color = {159, 159, 223}, thickness = 0.5, smooth = Smooth.None));
  connect(HTR.gasOut, r_HTR_hout.inlet);
  connect(r_HTR_hout.outlet, r_LTR_hin.inlet);
  connect(r_LTR_hin.outlet, LTR.gasIn);  
  connect(LTR.gasOut, r_LTR_hout.inlet);
  connect(r_LTR_hout.outlet, sink_hot.flange) annotation(
    Line(points = {{46, 0}, {46, 0}, {60, 0}}, color = {159, 159, 223}, thickness = 0.5));
*/

/*
  // Heater + turbine
  
  connect(source_cold.flange, r_heater_cin.inlet);
  connect(r_heater_cin.outlet, heater.gasIn);    
  connect(heater.gasOut, r_heater_cout.inlet);
  connect(r_heater_cout.outlet, Turbine1.inlet);
  
  connect(Turbine1.outlet, sens_turbine.inlet) annotation(
   Line(points = {{-50, 0}, {-20, 0}}, color = {159, 159, 223}, thickness = 0.5, smooth = Smooth.None));   
  connect(sens_turbine.outlet, sink_hot.flange) annotation(
    Line(points = {{46, 0}, {46, 0}, {60, 0}}, color = {159, 159, 223}, thickness = 0.5));
   
  connect(Turbine1.shaft_b, const_speed_comp.flange) annotation(
    Line(points = {{30, 0}, {74, 0}, {74, 0}, {74, 0}}));
  
  // connect(sens_turbine.outlet, sink_cold.flange);  
  // connect(source_hot.flange, HTR.gasIn);    

  // hot stream for heater
  connect(source_heater_hot.flange, r_heater_hin.inlet);
  connect(r_heater_hin.outlet, heater.waterIn);
  connect(heater.waterOut, r_heater_hout.inlet);
  connect(r_heater_hout.outlet, sink_heater_hot.flange);
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
  // turbine alone
  connect(source_cold.flange, Turbine1.inlet);  
  connect(Turbine1.outlet, sens_turbine.inlet) annotation(
   Line(points = {{-50, 0}, {-20, 0}}, color = {159, 159, 223}, thickness = 0.5, smooth = Smooth.None));   
  connect(sens_turbine.outlet, sink_hot.flange) annotation(
    Line(points = {{46, 0}, {46, 0}, {60, 0}}, color = {159, 159, 223}, thickness = 0.5));   
  connect(Turbine1.shaft_b, const_speed_comp.flange) annotation(
    Line(points = {{30, 0}, {74, 0}, {74, 0}, {74, 0}}));
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
*/
/*
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
  // Heater alone
  connect(source_heater_hot.flange, heater.waterIn);
  connect(heater.waterOut, T_waterOut.inlet);
  connect(T_waterOut.outlet, sink_heater_hot.flange);
  
  connect(source_cold.flange, heater.gasIn);
  connect(heater.gasOut, T_gasOut.inlet);
  connect(T_gasOut.outlet, sink_cold.flange);
 */  

annotation(
    Diagram(graphics),
    experiment(StartTime = 0, StopTime = 300, Tolerance = 1e-3, Interval = 10),
    __OpenModelica_commandLineOptions = "--matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection -d=initialization,NLSanalyticJacobian,aliasConflicts",    
    __OpenModelica_simulationFlags(lv = "LOG_DEBUG,LOG_NLS,LOG_NLS_V,LOG_STATS,LOG_INIT,LOG_STDOUT, -w", outputFormat = "mat", s = "dassl", nls = "homotopy"));
end TestTP_Components_PCHE;
