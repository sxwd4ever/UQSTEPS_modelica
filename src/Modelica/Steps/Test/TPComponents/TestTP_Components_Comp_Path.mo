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
  
  // package medium_main = Modelica.Media.IdealGases.SingleGases.CO2;
  package medium_main = Steps.Media.CO2;
  package medium_cooler = ThermoPower.Water.StandardWater;//Modelica.Media.IdealGases.SingleGases.CO2;  
  
  parameter Real Cf_C1                  = 1, Cf_C2 = 1, Cf_C3 = 1;
  
  parameter Model.RCBCycleConfig cfg(
    redeclare package medium_cooler = medium_cooler,
    redeclare package medium_main   = medium_main,
    
    r_i_heater  = 1e-3,
    r_t_heater  = 2e-3, 
    r_o_heater  = 3e-3,                      // agree with the final parameter Dhyd = 1 in HE, should be checked to see if it is capable of containing all fluid-metal tubes
    N_ch_heater = 10000,
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
    Ns_comp     = 30000,
    Ns_recomp   = 30000,
    N_ch_cooler = 50000,
    r_i_cooler  = 0.5e-3,
    r_t_cooler  = 0.7e-3,
    r_o_cooler  = 1e-3,  
        
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
    
    //mdot_main   = 128.774,
    mdot_main   = 125.774,
    // mdot_comp   = 88.0661,
    mdot_comp   = 85,
    mdot_heater = 40,
    mdot_cooler = 40.7188    
  );
  
  // set the values of parameters accordingly
  parameter Model.TurbomachineryConfig cfg_comp   = cfg.cfg_comp;
  parameter Model.TurbomachineryConfig cfg_recomp = cfg.cfg_recomp;
  parameter Model.HeatExchangerConfig cfg_cooler  = cfg.cfg_cooler;
  parameter Model.SplitterConfig cfg_splitter     = cfg.cfg_splitter;
  
  parameter Model.ThermoState st_source      = cfg_splitter.st_in;
  parameter Model.ThermoState st_sink        = cfg_comp.st_out;
  parameter Model.ThermoState st_sink_bypass = cfg_recomp.st_out;
  
  parameter Integer N_seg_cooler = cfg.cfg_cooler.cfg_hot.geo_path.N_seg; 
  
  inner ThermoPower.System system(allowFlowReversal = false, initOpt=ThermoPower.Choices.Init.Options.noInit) annotation(
    Placement(transformation(extent = {{80, 80}, {100, 100}})));
  
  
  ThermoPower.Gas.SourceMassFlow source(
  redeclare package Medium = medium_main,
    w0 = st_source.mdot,
    p0 = st_source.p,
    T  = st_source.T)
    annotation(
    Placement(transformation(origin = {0, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 270)));
    
  ThermoPower.Gas.SinkPressure sink(
  redeclare package Medium = medium_main,
    p0 = st_sink.p,
    T  = st_sink.T) 
    annotation(
    Placement(transformation(origin = {0, -80}, extent = {{-10, -10}, {10, 10}}, rotation = 270)));

  ThermoPower.Gas.SinkPressure sink_bypass(
  redeclare package Medium = medium_main,
    p0 = st_sink_bypass.p,
    T  = st_sink_bypass.T)
    annotation(
    Placement(transformation(origin = {0, -80}, extent = {{-10, -10}, {10, 10}}, rotation = 270)));    

  ThermoPower.Water.SourceMassFlow source_cooler(
    redeclare package Medium = medium_cooler,
    w0 = cfg_cooler.cfg_cold.st_in.mdot,
    p0 = cfg_cooler.cfg_cold.st_in.p,
    h  = cfg_cooler.cfg_cold.st_in.h,
    T  = cfg_cooler.cfg_cold.st_in.T) 
    annotation(
    Placement(transformation(origin = {0, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 270)));
    
  ThermoPower.Water.SinkPressure sink_cooler(
    redeclare package Medium = medium_cooler, 
    p0    = cfg_cooler.cfg_cold.st_out.p,
    T     = cfg_cooler.cfg_cold.st_out.T,
    use_T = false) 
    annotation(
    Placement(transformation(origin = {0, -80}, extent = {{-10, -10}, {10, 10}}, rotation = 270)));
  
  Steps.TPComponents.HE cooler(
    redeclare package FluidMedium = medium_cooler, 
    redeclare package FlueGasMedium = medium_main, 
       
    // Ngo 
    redeclare replaceable model HeatTransfer_F = Steps.TPComponents.NgoHeatTransferFV(),
    redeclare replaceable model HeatTransfer_G = Steps.TPComponents.NgoHeatTransferFV(),
    gasFlow(      
      heatTransfer(
        Cf_C1 = Cf_C1,
        Cf_C2 = Cf_C2,
        Cf_C3 = Cf_C3,
        gamma_min = 1500,
        gamma_max = 5000        
      ),
      Tstartin  = cfg_cooler.cfg_gas.st_in.T,
      Tstartout = cfg_cooler.cfg_gas.st_out.T,
      Nt        = cfg_cooler.cfg_gas.N_ch,
      Dhyd      = cfg_cooler.cfg_gas.geo_area.d_hyd
    ),
    fluidFlow(      
      heatTransfer(
        Cf_C1 = Cf_C1,
        Cf_C2 = Cf_C2,
        Cf_C3 = Cf_C3,
        gamma_min = 2000,
        gamma_max = 18000
      ), 
      fixedMassFlowSimplified = true,
      hstartin                = cfg_cooler.cfg_fluid.st_in.h,
      hstartout               = cfg_cooler.cfg_fluid.st_out.h
    ),    
  
    // // DittusBoelter
    // redeclare replaceable model HeatTransfer_F = ThermoPower.Thermal.HeatTransferFV.DittusBoelter,    
    // redeclare replaceable model HeatTransfer_G = ThermoPower.Thermal.HeatTransferFV.DittusBoelter,
    // fluidFlow(
    //   heatTransfer(heating=false, Re_min = 200), 
    //   fixedMassFlowSimplified = true,
    //   hstartin                = cfg_cooler.cfg_fluid.st_in.h,
    //   hstartout               = cfg_cooler.cfg_fluid.st_out.h),   // set the fluid flow as fixed mdot for simplarity
    // gasFlow(
    //   heatTransfer(heating=true, Re_min = 2300), 
    //   Tstartin  = cfg_cooler.cfg_gas.st_in.T,
    //   Tstartout = cfg_cooler.cfg_gas.st_out.T,
    //   Nt        = cfg_cooler.cfg_gas.N_ch,
    //   Dhyd      = cfg_cooler.cfg_gas.geo_area.d_hyd),    

    // redeclare replaceable model HeatTransfer_F = ThermoPower.Thermal.HeatTransferFV.IdealHeatTransfer,    
    // redeclare replaceable model HeatTransfer_G = ThermoPower.Thermal.HeatTransferFV.IdealHeatTransfer, 
    // fluidFlow(
    //   fixedMassFlowSimplified = true,
    //   hstartin                = cfg_cooler.cfg_cold.st_in.h,
    //   hstartout               = cfg_cooler.cfg_cold.st_out.h),   // set the fluid flow as fixed mdot for simplarity
    // gasFlow(
    //   Tstartin    = cfg_cooler.cfg_hot.st_in.T,
    //   Tstartout   = cfg_cooler.cfg_hot.st_out.T),
    
    // redeclare replaceable model HeatTransfer_F = ThermoPower.Thermal.HeatTransferFV.ConstantHeatTransferCoefficient(),      
    // redeclare replaceable model HeatTransfer_G = ThermoPower.Thermal.HeatTransferFV.ConstantHeatTransferCoefficient(), 
        
    // fluidFlow(
    //   heatTransfer(gamma = cfg_cooler.cfg_cold.gamma_HE),
    //   fixedMassFlowSimplified = true,
    //   hstartin                = cfg_cooler.cfg_cold.st_in.h,
    //   hstartout               = cfg_cooler.cfg_cold.st_out.h),   // set the fluid flow as fixed mdot for simplarity
    // gasFlow(
    //   heatTransfer(gamma = cfg_cooler.cfg_hot.gamma_HE),
    //   Tstartin    = cfg_cooler.cfg_hot.st_in.T,
    //   Tstartout   = cfg_cooler.cfg_hot.st_out.T),   
           
    redeclare model HeatExchangerTopology = ThermoPower.Thermal.HeatExchangerTopologies.CounterCurrentFlow,
    cfg             = cfg_cooler,
    N_G             = N_seg_cooler,
    N_F             = N_seg_cooler,
    // Nt = cfg_cooler.cfg_cold.N_ch,
    SSInit          = false,
    gasQuasiStatic  = true,
    FluidPhaseStart = ThermoPower.Choices.FluidPhase.FluidPhases.Liquid,
    metalTube(WallRes=false))
  annotation(
    Placement(visible = true, transformation(origin = {56, -46}, extent = {{-8, -8}, {8, 8}}, rotation = 0)));
 /*
  Steps.TPComponents.HE cooler(
    redeclare package FluidMedium = medium_cooler, 
    redeclare package FlueGasMedium = medium_main, 
    
    // // DittusBoelter
    // redeclare replaceable model HeatTransfer_F = ThermoPower.Thermal.HeatTransferFV.DittusBoelter,    
    // redeclare replaceable model HeatTransfer_G = ThermoPower.Thermal.HeatTransferFV.DittusBoelter,
    // fluidFlow(
    //   heatTransfer(heating=true), 
    //   fixedMassFlowSimplified = true,
    //   hstartin                = cfg_cooler.cfg_hot.st_in.h,
    //   hstartout               = cfg_cooler.cfg_hot.st_out.h),   // set the fluid flow as fixed mdot for simplarity
    // gasFlow(
    //   heatTransfer(heating=false), 
    //   Tstartin  = cfg_cooler.cfg_cold.st_in.T,
    //   Tstartout = cfg_cooler.cfg_cold.st_out.T),
    
    redeclare replaceable model HeatTransfer_F = ThermoPower.Thermal.HeatTransferFV.IdealHeatTransfer,    
    redeclare replaceable model HeatTransfer_G = ThermoPower.Thermal.HeatTransferFV.IdealHeatTransfer, 
    fluidFlow(
      fixedMassFlowSimplified = true,
      hstartin                = cfg_cooler.cfg_cold.st_in.h,
      hstartout               = cfg_cooler.cfg_cold.st_out.h),   // set the fluid flow as fixed mdot for simplarity
    gasFlow(
      Tstartin    = cfg_cooler.cfg_hot.st_in.T,
      Tstartout   = cfg_cooler.cfg_hot.st_out.T),
    
    // redeclare replaceable model HeatTransfer_F = ThermoPower.Thermal.HeatTransferFV.ConstantHeatTransferCoefficientTwoGrids(gamma = thermo_hot.gamma_he),          
    // redeclare replaceable model HeatTransfer_G = ThermoPower.Thermal.HeatTransferFV.ConstantHeatTransferCoefficient(gamma =  thermo_cold.gamma_he),           
    redeclare model HeatExchangerTopology = ThermoPower.Thermal.HeatExchangerTopologies.CounterCurrentFlow,
    cfg             = cfg_cooler,
    N_G             = N_seg_cooler,
    N_F             = N_seg_cooler,
    SSInit          = false,
    gasQuasiStatic  = true,
    FluidPhaseStart = ThermoPower.Choices.FluidPhase.FluidPhases.Liquid,
    metalTube(WallRes=false))
    annotation(
      Placement(transformation(extent = {{-20, -20}, {20, 20}}, rotation = 0)));
  */
  
  
    Steps.TPComponents.GasStateReader r_cooler_hin(redeclare package Medium = medium_main) annotation(
      Placement(visible = true, transformation(origin = {66, -56}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
    Steps.TPComponents.GasStateReader r_cooler_hout(redeclare package Medium = medium_main) annotation(
      Placement(visible = true, transformation(origin = {96, 44}, extent = {{-6, -6}, {6, 6}}, rotation = 180)));
    Steps.TPComponents.WaterStateReader r_cooler_cin(redeclare package Medium = medium_cooler) annotation(
      Placement(visible = true, transformation(origin = {-22, 50}, extent = {{-6, -6}, {6, 6}}, rotation = -90)));
    Steps.TPComponents.WaterStateReader r_cooler_cout(redeclare package Medium = medium_cooler) annotation(
      Placement(visible = true, transformation(origin = {-32, 20}, extent = {{-6, -6}, {6, 6}}, rotation = 180)));
 
  ThermoPower.Gas.Compressor compressor(
  // Steps.TPComponents.CompressorFixedP compressor(
    redeclare package Medium = medium_main,    
    pstart_in  = cfg_comp.st_in.p,
    pstart_out = cfg_comp.st_out.p,
    Tstart_in  = cfg_comp.st_in.T,
    Tstart_out = cfg_comp.st_out.T,
    tablePhic  = tablePhic_comp_mc,
    tableEta   = tableEta_comp_mc,
    tablePR    = tablePR_comp_mc,
    Table      = ThermoPower.Choices.TurboMachinery.TableTypes.matrix,
    Ndesign    = cfg_comp.N,
    Tdes_in    = cfg_comp.st_in.T,
    gas_in(
      p(start = cfg_comp.st_in.p, nominal = cfg_comp.st_in.p), 
      T(start = cfg_comp.st_in.T, nominal = cfg_comp.st_in.T)),    
    gas_iso(
      p(nominal = cfg_comp.st_in.p), 
      T(nominal = cfg_comp.st_in.T),
      h(start   = cfg_comp.st_out.h)))
  annotation (Placement(transformation(extent={{-20,-20},{
            20,20}}, rotation=0)));
  
  /*       
  Modelica.Mechanics.Rotational.Sources.ConstantSpeed ConstantSpeed_mc(
    w_fixed    = cfg_comp.N,
    useSupport = false)
  annotation (Placement(transformation(
        extent={{-50,-10},{-30,10}}, rotation=0)));
  */
  
  Steps.TPComponents.GasStateReader r_comp_in(redeclare package Medium = medium_main); 

  Steps.TPComponents.FixedPController FPC_mc(
    redeclare package Medium   = medium_main,
    tablePhic                  = tablePhic_comp_mc,
    tablePR                    = tablePR_comp_mc,
    Table                      = ThermoPower.Choices.TurboMachinery.TableTypes.matrix,
    Ndesign                    = cfg_comp.N,
    Tdes_in                    = cfg_comp.st_in.T,
    T_inlet(start = cfg_comp.st_in.T, nominal = cfg_comp.st_in.T),
    p_inlet(start = cfg_comp.st_in.p, nominal = cfg_comp.st_in.p),
    mdot_inlet(start = cfg_comp.st_in.mdot, nominal = cfg_comp.st_in.mdot)
  ) "Fixed pressure controller for main compressor";
  
  Modelica.Mechanics.Rotational.Sources.Speed speed_mc(
    exact = false,
    w_ref(start = cfg_comp.N, nominal= cfg_comp.N),
    w(start = cfg_comp.N, nominal = cfg_comp.N)
  );  


  ThermoPower.Gas.Compressor recompressor(
  // Steps.TPComponents.CompressorFixedP recompressor(
    redeclare package Medium = medium_main,    
    pstart_in  = cfg_recomp.st_in.p,
    pstart_out = cfg_recomp.st_out.p,
    Tstart_in  = cfg_recomp.st_in.T,
    Tstart_out = cfg_recomp.st_out.T,
    tablePhic  = tablePhic_comp_rc,
    tableEta   = tableEta_comp_rc,
    tablePR    = tablePR_comp_rc,
    Table      = ThermoPower.Choices.TurboMachinery.TableTypes.matrix,
    Ndesign    = cfg_recomp.N,
    Tdes_in    = cfg_recomp.st_in.T,
    gas_in(
      p(start = cfg_recomp.st_in.p, nominal = cfg_recomp.st_in.p), 
      T(start = cfg_recomp.st_in.T, nominal = cfg_recomp.st_in.T)),    
    gas_iso(
      p(nominal = cfg_recomp.st_in.p), 
      T(nominal = cfg_recomp.st_in.T),
      h(start   = cfg_recomp.st_out.h)))
  annotation (Placement(transformation(extent={{-20,-20},{
            20,20}}, rotation=0)));

  /*          
  Modelica.Mechanics.Rotational.Sources.ConstantSpeed ConstantSpeed_rc(
    w_fixed    = cfg_recomp.N,
    useSupport = false)
  annotation (Placement(transformation(
      extent={{-50,-10},{-30,10}}, rotation=0)));
  */

  Steps.TPComponents.GasStateReader r_recomp_in(redeclare package Medium = medium_main); 
  
  Steps.TPComponents.FixedPController FPC_rc(
    redeclare package Medium   = medium_main,
    tablePhic                  = tablePhic_comp_rc,
    tablePR                    = tablePR_comp_rc,
    Table                      = ThermoPower.Choices.TurboMachinery.TableTypes.matrix,
    Ndesign                    = cfg_recomp.N,
    Tdes_in                    = cfg_recomp.st_in.T,
    T_inlet(start = cfg_recomp.st_in.T, nominal = cfg_recomp.st_in.T),
    p_inlet(start = cfg_recomp.st_in.p, nominal = cfg_recomp.st_in.p),
    mdot_inlet(start = cfg_recomp.st_in.mdot, nominal = cfg_recomp.st_in.mdot)
  ) "Fixed pressure controller for recompressor";
  
  Modelica.Mechanics.Rotational.Sources.Speed speed_rc(
    exact = false,
    w_ref(start = cfg_recomp.N, nominal= cfg_recomp.N),
    w(start = cfg_recomp.N, nominal = cfg_recomp.N)
  );    
  
  ThermoPower.Gas.FlowSplit splitter(redeclare package Medium = medium_main);
  
protected

  parameter Real tableEta_comp_mc[5, 4]  = [0, 95, 100, 105; 1, 0.85310219, 0.837591241, 0.832420925; 2, 0.868613139, 0.857238443, 0.851034063; 3, 0.860340633, 0.85, 0.842761557; 4, 0.85310219, 0.839659367, 0.816909976];
  parameter Real tablePhic_comp_mc[5, 4] = [0, 95, 100, 105; 1, 0.000134346, 0.000150832, 0.000164161; 2, 0.000137854, 0.000153638, 0.00016802; 3, 0.000142414, 0.000158549, 0.000169774; 4, 0.000145921, 0.000161706, 0.000171528];
  parameter Real tablePR_comp_mc[5, 4]   = [0, 95, 100, 105; 1, 1.967529638, 2.350588505, 2.785882673; 2, 1.915294338, 2.315764972, 2.681412073; 3, 1.810823737, 2.220000255, 2.524706172; 4, 1.654117837, 2.115529655, 2.359294389];

  // performance map for re compressor      
  parameter Real tableEta_comp_rc[5, 4]  = [0, 95, 100, 105; 1, 0.85310219, 0.837591241, 0.832420925; 2, 0.868613139, 0.857238443, 0.851034063; 3, 0.860340633, 0.85, 0.842761557; 4, 0.85310219, 0.839659367, 0.816909976];
  parameter Real tablePhic_comp_rc[5, 4] = [0, 95, 100, 105; 1, 7.17663E-05, 8.05731E-05, 8.76935E-05; 2, 7.36401E-05, 8.20721E-05, 8.97547E-05; 3, 7.6076E-05, 8.46954E-05, 9.06916E-05; 4, 7.79498E-05, 8.63819E-05, 9.16285E-05];
  parameter Real tablePR_comp_rc[5, 4]   = [0, 95, 100, 105; 1, 1.967529638, 2.350588505, 2.785882673; 2, 1.915294338, 2.315764972, 2.681412073; 3, 1.810823737, 2.220000255, 2.524706172; 4, 1.654117837, 2.115529655, 2.359294389];

initial equation
//hstart_F_Out = cooler.waterOut.h_outflow;
equation

/*
  // main compressor alone
  
  connect(source.flange, compressor.inlet);
  connect(compressor.outlet, sink.flange);
  
  connect(compressor.shaft_a, ConstantSpeed_mc.flange);
*/
/*
  // re-compressor alone
  
  connect(source.flange, recompressor.inlet);
  connect(recompressor.outlet, sink.flange);
  
  connect(recompressor.shaft_a, ConstantSpeed_rc.flange);
*/
/*
  // Cooler alone
  connect(source.flange, r_cooler_hin.inlet);
  connect(r_cooler_hin.outlet, cooler.gasIn);
  connect(cooler.gasOut, r_cooler_hout.inlet);
  connect(r_cooler_hout.outlet, sink.flange);
  
  connect(source_cooler.flange, r_cooler_cin.inlet);
  connect(r_cooler_cin.outlet, cooler.waterIn);
  connect(cooler.waterOut, r_cooler_cout.inlet);
  connect(r_cooler_cout.outlet, sink_cooler.flange);  
*/

/*
  // cooler + compressor
  // hot side 
  connect(source.flange, r_cooler_hin.inlet);
  connect(r_cooler_hin.outlet, cooler.gasIn);
  connect(cooler.gasOut, r_cooler_hout.inlet);
  connect(r_cooler_hout.outlet, compressor.inlet);
  connect(compressor.outlet, sink.flange);  
  
  connect(compressor.shaft_a, ConstantSpeed_mc.flange);
  
  // cold side
  connect(source_cooler.flange, r_cooler_cin.inlet);
  connect(r_cooler_cin.outlet, cooler.waterIn);
  connect(cooler.waterOut, r_cooler_cout.inlet);
  connect(r_cooler_cout.outlet, sink_cooler.flange);
*/
 
  // splitter + cooler + compressor + recompressor
  // hot side / main path
  connect(source.flange, splitter.inlet);  
  connect(splitter.outlet1, r_cooler_hin.inlet);
  
  connect(r_cooler_hin.outlet, cooler.gasIn);
  connect(cooler.gasOut, r_cooler_hout.inlet);
  connect(r_cooler_hout.outlet, r_comp_in.inlet);
  connect(r_comp_in.outlet, compressor.inlet);

    connect(r_comp_in.p, FPC_mc.p_inlet);
    connect(r_comp_in.T, FPC_mc.T_inlet);
    connect(r_comp_in.w, FPC_mc.mdot_inlet);    
    connect(FPC_mc.omega, speed_mc.w_ref);

    connect(speed_mc.flange, compressor.shaft_a);
  
  connect(compressor.outlet, sink.flange);  
  
    // connect(compressor.shaft_a, ConstantSpeed_mc.flange);
    // connect(compressor.speed.flange, compressor.shaft_a);    
    
  // bypass path
  
  connect(splitter.outlet2, r_recomp_in.inlet);
  connect(r_recomp_in.outlet, recompressor.inlet);
  
    connect(r_recomp_in.p, FPC_rc.p_inlet);
    connect(r_recomp_in.T, FPC_rc.T_inlet);
    connect(r_recomp_in.w, FPC_rc.mdot_inlet);    
    connect(FPC_rc.omega, speed_rc.w_ref);
    connect(speed_rc.flange, recompressor.shaft_a);
    
  connect(recompressor.outlet, sink_bypass.flange);
  
    //connect(recompressor.speed.flange, recompressor.shaft_a);
    
    //connect(recompressor.shaft_a, ConstantSpeed_rc.flange);
  
  // cold side
  connect(source_cooler.flange, r_cooler_cin.inlet);
  connect(r_cooler_cin.outlet, cooler.waterIn);
  connect(cooler.waterOut, r_cooler_cout.inlet);
  connect(r_cooler_cout.outlet, sink_cooler.flange);

annotation(
    Diagram(graphics),
    experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-3, Interval = 1),
    __OpenModelica_commandLineOptions = "--matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection -d=initialization,NLSanalyticJacobian,aliasConflicts,bltdump",    
    __OpenModelica_simulationFlags(lv = "LOG_DEBUG,LOG_NLS,LOG_NLS_V,LOG_STATS,LOG_INIT,LOG_STDOUT, -w", outputFormat = "mat", s = "dassl", nls = "homotopy"));
end TestTP_Components_Comp_Path;
