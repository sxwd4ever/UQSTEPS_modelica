within Steps.Test.TPComponents.Transient;

model TestDyn_Comps_Comp_Path
  "Test for combination of components (Cooler + Spliter + Main/Re compressor) in ThermoPower in transient simulation"  
    
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
  
  package medium_main = Steps.Media.SCO2;
  package medium_cooler = ThermoPower.Water.StandardWater;//Modelica.Media.IdealGases.SingleGases.CO2;  
  
  // geometry parameters
  constant  Real pi           = Modelica.Constants.pi;
  // parameter Integer N_ch      = integer(2400*1e1) "channel number";
  parameter Integer N_ch      = integer(2e6) "channel number";
  parameter Integer N_seg     = 10 "number of segments in one tube";
  parameter SI.Length D_ch    = 1.72e-3  "channel diameter, semi circular tube";
  parameter SI.Length r_ch    = D_ch / 2 "channel radiaus";
  parameter SI.Length L_fp    = 270 * 1e-3 * 2 "equivalent valid channel flow path length";
  parameter SI.Length L_pitch = 12e-3 "pitch length";
  parameter Real a_phi        = 0 "pitch angle °";
  parameter SI.Length H_ch    = 4.17e-3 "Height of the solid domain, containing one cold tube and one hot tube";
  parameter SI.Length W_ch    = 2.3e-3 "Width of the solid domain";
  parameter SI.Length T_wall  = 0.51e-3 "Wall thinckness";
  parameter SI.Length L_wall  = 420e-3 * 2 "Length of wall, not necessarily equals to length of flow path";
  parameter SI.Area A         = pi * r_ch ^2 / 2 "Area of cross section of semi circular tube";
  
  // Stainless 316, 316L, 317, 317L
  parameter Modelica.SIunits.Density rho_wall             = 8030 "density of wall, kg/m3";
  parameter Modelica.SIunits.SpecificHeatCapacity cp_wall = 485 "cp of wall, J/kg-K";
  // thermal conductivity (T in K) https://www.theworldmaterial.com/aisi-316-ss316-stainless-steel-properties-composition/
  // parameter Real table_k_metalwall[:,:] = [20, 12.1; 100, 16.3; 500, 21.5];
  parameter Real table_k_metalwall[:,:] = [293.15, 12.1; 373.15, 16.3; 773.15, 21.5];
  parameter Real Cf_C1                  = 1.626, Cf_C2 = 1, Cf_C3 = 1;
  // parameter Real Cf_C1_cold = 1, Cf_C2_cold = 1, Cf_C3_cold = 1;
  parameter Real use_rho_bar  = -1.0;
  parameter Real rho_bar_hot  = 1.0;
  parameter Real rho_bar_cold = 1.0;

  // input parameters of the power block
  parameter Modelica.SIunits.MassFlowRate mdot_main = 125 "kg/s, mass flow in the main path of PB, which follows the power demand";
  parameter Modelica.SIunits.MassFlowRate mdot_heater_hot = 90 "kg/s, mass flow rate of heater's hot fluid";
  parameter Real gamma = 0.4 "split ratio, mdot_bypass/mdot_main";
  
  parameter Modelica.SIunits.Temperature T_heater_hot = from_degC(800) "K, Temperature of heater's hot fluid";  
  parameter Modelica.SIunits.Temperature T_cooler_cold = from_degC(45) "K, Temperature of cooler's cold fluid";
  parameter Integer T_step = 3;
  
  parameter Model.RCBCycleConfig cfg(
    redeclare package medium_cooler = medium_cooler,
    redeclare package medium_main   = medium_main,
    r_i_heater  = 1e-3,
    r_t_heater  = 2e-3, //cfg.r_i_heater + 10e-3,
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
    N_ch_cooler = 50000,
    r_i_cooler  = 0.5e-3,
    r_t_cooler  = 0.7e-3,
    r_o_cooler  = 1e-3,    
    table_k_LTR_wall = table_k_metalwall,
    table_k_HTR_wall = table_k_metalwall,    
    // latest boundary conditions, following values are simulation results with sourcePressure.p = 200 bar and above geometry params
    p_comp_in  = 109.59e5,
    p_comp_out = 20e6,    
    p_heater   = 20e6,    
    T_HTR_hot_in      = from_degC(556.322),
    T_HTR_cold_out    = from_degC(515.234),
    T_HTR_hot_out     = from_degC(296.146),
    T_HTR_cold_in     = from_degC(266.181),
    T_LTR_cold_in     = from_degC(118.13),
    T_LTR_hot_out     = from_degC(154.909),
    T_heater_hot_in   = from_degC(800),
    T_heater_hot_out  = from_degC(707.918),
    T_heater_cold_out = from_degC(637.551),
    T_cooler_cold_out = from_degC(61.7479462700001),
    T_cooler_hot_out  = from_degC(68.494),
    
    mdot_main   = 128.774,
    mdot_comp   = 88.0661,
    mdot_heater = 40,
    mdot_cooler = 40.7188
  );
  
  // set the values of parameters accordingly
  parameter Model.TurbomachineryConfig cfg_comp   = cfg.cfg_comp;
  parameter Model.TurbomachineryConfig cfg_recomp = cfg.cfg_recomp;
  parameter Model.HeatExchangerConfig cfg_cooler  = cfg.cfg_cooler;
  parameter Model.SplitterConfig cfg_splitter     = cfg.cfg_splitter;
  
  parameter Model.ThermoState st_source      = cfg_splitter.st_in;
  parameter Model.ThermoState st_sink        = cfg_comp.st_out;//cfg_cooler.cfg_hot.st_out;
  parameter Model.ThermoState st_sink_bypass = cfg_recomp.st_out;
  
  parameter Integer N_seg_cooler = cfg.cfg_cooler.cfg_hot.geo_path.N_seg; 
  
  parameter Boolean SSInit = true "Steady-state initialization";  
  inner ThermoPower.System system(
    allowFlowReversal = false,
    initOpt           = ThermoPower.Choices.Init.Options.steadyState)
    annotation(
    Placement(transformation(extent = {{80, 80}, {100, 100}})));  
  
 
  ThermoPower.Gas.SourceMassFlow source(
    redeclare package Medium = medium_main, 
    T        = st_source.T,
    p0       = st_source.p,
    w0       = st_source.mdot,
    use_in_T = false,
    use_in_w0 = true,
    gas(
      p(start = st_source.p, nominal = st_source.p), 
      T(start = st_source.T, nominal = st_source.T))) 
    annotation(
    Placement(transformation(origin = {0, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 270)));
 
/*  
  ThermoPower.Gas.SourcePressure source(
    redeclare package Medium = medium_main, 
    T        = st_source.T,
    p0       = st_source.p,
    use_in_T = true,
    gas(
      p(start = st_source.p, nominal = st_source.p), 
      T(start = st_source.T, nominal = st_source.T))); 
*/
     
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
       
    // DittusBoelter
    redeclare replaceable model HeatTransfer_F = ThermoPower.Thermal.HeatTransferFV.DittusBoelter,    
    redeclare replaceable model HeatTransfer_G = ThermoPower.Thermal.HeatTransferFV.DittusBoelter,
    fluidFlow(
      heatTransfer(heating=false, Re_min = 200), 
      fixedMassFlowSimplified = true,
      hstartin                = cfg_cooler.cfg_fluid.st_in.h,
      hstartout               = cfg_cooler.cfg_fluid.st_out.h),   // set the fluid flow as fixed mdot for simplarity
    gasFlow(
      heatTransfer(heating=true, Re_min = 2300), 
      Tstartin  = cfg_cooler.cfg_gas.st_in.T,
      Tstartout = cfg_cooler.cfg_gas.st_out.T,
      Nt        = cfg_cooler.cfg_gas.N_ch,
      Dhyd      = cfg_cooler.cfg_gas.geo_area.d_hyd),    
    
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
    SSInit          = true,
    gasQuasiStatic  = false,
    FluidPhaseStart = ThermoPower.Choices.FluidPhase.FluidPhases.Liquid,
    metalTube(WallRes=false))
    annotation(
      Placement(transformation(extent = {{-20, -20}, {20, 20}}, rotation = 0)));
  
    Steps.TPComponents.GasStateReader r_cooler_hin(redeclare package Medium = medium_main) annotation(
      Placement(visible = true, transformation(origin = {66, -56}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
    Steps.TPComponents.GasStateReader r_cooler_hout(redeclare package Medium = medium_main) annotation(
      Placement(visible = true, transformation(origin = {96, 44}, extent = {{-6, -6}, {6, 6}}, rotation = 180)));
    Steps.TPComponents.WaterStateReader r_cooler_cin(redeclare package Medium = medium_cooler) annotation(
      Placement(visible = true, transformation(origin = {-22, 50}, extent = {{-6, -6}, {6, 6}}, rotation = -90)));
    Steps.TPComponents.WaterStateReader r_cooler_cout(redeclare package Medium = medium_cooler) annotation(
      Placement(visible = true, transformation(origin = {-32, 20}, extent = {{-6, -6}, {6, 6}}, rotation = 180)));
 
  ThermoPower.Gas.Compressor compressor(
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
    gas_iso(
      p(nominal = cfg_comp.st_in.p), 
      T(nominal = cfg_comp.st_in.T),
      h(start   = cfg_comp.st_out.h)))
  annotation (Placement(transformation(extent={{-20,-20},{
            20,20}}, rotation=0)));
            
  Modelica.Mechanics.Rotational.Sources.ConstantSpeed ConstantSpeed_mc(
    w_fixed    = cfg_comp.N,
    useSupport = false)
  annotation (Placement(transformation(
        extent={{-50,-10},{-30,10}}, rotation=0)));

  ThermoPower.Gas.Compressor recompressor(
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
    gas_iso(
      p(nominal = cfg_recomp.st_in.p), 
      T(nominal = cfg_recomp.st_in.T),
      h(start   = cfg_recomp.st_out.h)))
  annotation (Placement(transformation(extent={{-20,-20},{
            20,20}}, rotation=0)));
            
  Modelica.Mechanics.Rotational.Sources.ConstantSpeed ConstantSpeed_rc(
    w_fixed    = cfg_recomp.N,
    useSupport = false)
  annotation (Placement(transformation(
      extent={{-50,-10},{-30,10}}, rotation=0)));
  
  ThermoPower.Gas.FlowSplit splitter(redeclare package Medium = medium_main);

  // ramp component for transient test
  Modelica.Blocks.Sources.Ramp ramp1(
  final height    = st_source.mdot * 0.1,
  final duration  = 9 * 60, // 9 mins
  final startTime = 3,
  final offset    = st_source.mdot);
  
  /*
  // hot/gas side
  Modelica.Blocks.Sources.IntegerConstant const_T_offset_h(k = integer(st_source.T)) annotation(
    Placement(visible = true, transformation(origin = {-128, 72}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Interaction.Show.IntegerValue disp_T_h annotation(
    Placement(visible = true, transformation(origin = {-144, 158}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.MathInteger.TriggeredAdd triadd_T_h annotation(
    Placement(visible = true, transformation(origin = {-128, 126}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
  Modelica.Blocks.Sources.BooleanPulse en_triadd_T_h(period = 10, startTime = 3, width = 10) annotation(
    Placement(visible = true, transformation(origin = {-156, 102}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.MathInteger.Sum sum_T_h(nu = 2) annotation(
    Placement(visible = true, transformation(origin = {-94, 126}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Interaction.Show.IntegerValue disp_T_step_h annotation(
    Placement(visible = true, transformation(origin = {-56, 158}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Math.IntegerToReal I2R_T_h annotation(
    Placement(visible = true, transformation(origin = {-48, 124}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.IntegerConstant const_T_step_h(k = T_step) annotation(
    Placement(visible = true, transformation(origin = {-194, 126}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  */
  
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
  connect(r_cooler_hout.outlet, compressor.inlet);
  connect(compressor.outlet, sink.flange);  
    connect(compressor.shaft_a, ConstantSpeed_mc.flange);    
    
  // bypass path
  connect(splitter.outlet2, recompressor.inlet);
  connect(recompressor.outlet, sink_bypass.flange);
    connect(recompressor.shaft_a, ConstantSpeed_rc.flange);
  
  // cold side
  connect(source_cooler.flange, r_cooler_cin.inlet);
  connect(r_cooler_cin.outlet, cooler.waterIn);
  connect(cooler.waterOut, r_cooler_cout.inlet);
  connect(r_cooler_cout.outlet, sink_cooler.flange);
  
  /*
  // hot / gas side 
  connect(en_triadd_T_h.y, triadd_T_h.trigger) annotation(
    Line(points = {{-145, 102}, {-132, 102}, {-132, 119}}, color = {255, 0, 255}));
  connect(triadd_T_h.y, sum_T_h.u[1]) annotation(
    Line(points = {{-121, 126}, {-104, 126}}, color = {255, 127, 0}));
  connect(triadd_T_h.y, disp_T_h.numberPort) annotation(
    Line(points = {{-121, 126}, {-160.5, 126}, {-160.5, 158}, {-155.5, 158}}, color = {255, 127, 0}));
  connect(sum_T_h.y, I2R_T_h.u) annotation(
    Line(points = {{-82.5, 126}, {-71.25, 126}, {-71.25, 124}, {-60, 124}}, color = {255, 127, 0}));
  connect(sum_T_h.y, disp_T_step_h.numberPort) annotation(
    Line(points = {{-82.5, 126}, {-71, 126}, {-71, 158}, {-67.5, 158}}, color = {255, 127, 0}));
  connect(const_T_step_h.y, triadd_T_h.u) annotation(
    Line(points = {{-183, 126}, {-136, 126}}, color = {255, 127, 0}));
  connect(const_T_offset_h.y, sum_T_h.u[2]) annotation(
    Line(points = {{-117, 72}, {-104, 72}, {-104, 126}}, color = {255, 127, 0}));
  connect(I2R_T_h.y, source.in_T) annotation(
    Line(points = {{-37, 124}, {6, 124}, {6, 64}}, color = {0, 0, 127}));  
  */
  
  connect(ramp1.y, source.in_w0);
annotation(
    Diagram(graphics),
    experiment(StartTime = 0, StopTime = 600, Tolerance = 1e-3, Interval = 6),
    __OpenModelica_commandLineOptions = "--matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection -d=initialization,NLSanalyticJacobian,aliasConflicts,bltdump",    
    __OpenModelica_simulationFlags(lv = "LOG_DEBUG,LOG_NLS,LOG_NLS_V,LOG_STATS,LOG_INIT,LOG_STDOUT, -w", outputFormat = "mat", s = "dassl", nls = "homotopy"));
end TestDyn_Comps_Comp_Path;
