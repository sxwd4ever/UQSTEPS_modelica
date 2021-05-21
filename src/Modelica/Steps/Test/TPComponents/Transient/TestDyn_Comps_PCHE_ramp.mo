within Steps.Test.TPComponents.Transient;

model TestDyn_Comps_PCHE_ramp
  "Test for ThermoPower based components (heater + turbine + HTR + LTR, mainly PCHE) with ramp change in mass flow rate"  
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
  
  // package medium_main = Modelica.Media.IdealGases.SingleGases.CO2; //Steps.Media.CO2;
  package medium_main = Steps.Media.SCO2(
    // inputChoice = ExternalMedia.Common.InputChoice.pT,
    substanceNames = {"CO2|debug=40"}
  );  
  // package medium_main = ExternalMedia.Examples.CO2CoolProp;
  // package medium_heater = Steps.Media.ThermiaOilD; // out of working range of this 10Mw high T loop
  package medium_heater = SolarTherm.Media.Sodium.Sodium_pT;

  
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
  parameter Integer T_step = 10;
 
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
    mdot_heater = 40
    // T_HTR_cold_in = from_degC(216.593),
    // T_HTR_hot_out = from_degC(263.249)
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
  parameter Boolean SSInit = true "Steady-state initialization";  
  inner ThermoPower.System system(
    allowFlowReversal = false,
    initOpt           = ThermoPower.Choices.Init.Options.steadyState)
    annotation(
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
    
    
    // redeclare replaceable model HeatTransfer_F = ThermoPower.Thermal.HeatTransferFV.IdealHeatTransfer,    
    // redeclare replaceable model HeatTransfer_G = ThermoPower.Thermal.HeatTransferFVH.IdealHeatTransfer, 
    fluidFlow(
      heatTransfer(gamma = cfg_heater.cfg_hot.gamma_HE),
      fixedMassFlowSimplified = true,
      hstartin                = cfg_heater.cfg_hot.st_in.h,
      hstartout               = cfg_heater.cfg_hot.st_out.h),   // set the fluid flow as fixed mdot for simplarity
    gasFlow(
      heatTransfer(gamma = cfg_heater.cfg_cold.gamma_HE),
      Tstartin    = cfg_heater.cfg_cold.st_in.T,
      Tstartout   = cfg_heater.cfg_cold.st_out.T),
    
    redeclare replaceable model HeatTransfer_F = ThermoPower.Thermal.HeatTransferFV.ConstantHeatTransferCoefficientTwoGrids(),     
    // redeclare replaceable model HeatTransfer_G = ThermoPower.Thermal.HeatTransferFV.IdealHeatTransfer,      
    redeclare replaceable model HeatTransfer_G = ThermoPower.Thermal.HeatTransferFV.ConstantHeatTransferCoefficient(),       

    redeclare model HeatExchangerTopology = ThermoPower.Thermal.HeatExchangerTopologies.CounterCurrentFlow,
    cfg             = cfg_heater,
    N_G             = N_seg_heater,
    N_F             = N_seg_heater,
    SSInit          = true,
    gasQuasiStatic  = false,
    FluidPhaseStart = ThermoPower.Choices.FluidPhase.FluidPhases.Liquid
  )
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

/*
  ThermoPower.Gas.SourceMassFlow source_cold(
    redeclare package Medium = medium_main, 
    T        = st_source_cold.T,
    p0       = st_source_cold.p,
    use_in_T = true,
    w0       = st_source_cold.mdot)
  annotation(
    Placement(transformation(origin = {0, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 270)));
  
*/

  ThermoPower.Gas.SourcePressure source_cold(
    redeclare package Medium = medium_main, 
    T        = st_source_cold.T,
    p0       = st_source_cold.p,
    use_in_T = true,
    gas(
      p(start = st_source_cold.p, nominal = st_source_cold.p), 
      T(start = st_source_cold.T, nominal = st_source_cold.T))) 
  annotation(
    Placement(transformation(origin = {0, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 270)));

  ThermoPower.Gas.SourceMassFlow source_mixer_in(
    redeclare package Medium = medium_main,
    T        = st_bypass.T,
    p0       = st_bypass.p,
    use_in_T = false,
    w0       = st_bypass.mdot
  );  

/*  
  ThermoPower.Gas.SinkPressure sink_cold(
    redeclare package Medium = medium_main, 
    p0 = st_sink_cold.p,
    T  = st_sink_cold.T)
  annotation(
    Placement(transformation(origin = {0, -80}, extent = {{-10, -10}, {10, 10}}, rotation = 270)));
*/
  ThermoPower.Gas.SinkPressure sink_hot(
    redeclare package Medium = medium_main,
    T  = st_sink_hot.T,
    p0 = st_sink_hot.p)
  annotation(
    Placement(transformation(extent = {{60, -10}, {80, 10}}, rotation = 0)));
/*
  ThermoPower.Gas.SourceMassFlow source_hot(
    redeclare package Medium = medium_main, 
    T        = st_source_hot.T,
    p0       = st_source_hot.p,
    w0       = st_source_hot.mdot,
    use_in_T = false)
  annotation(
    Placement(transformation(extent = {{-70, -10}, {-50, 10}}, rotation = 0))); 
*/
  // use FlowJoin to mix flow
  Gas.FlowJoin mixer(redeclare package Medium = medium_main);  

  Steps.TPComponents.PCHE LTR(
    redeclare package FluidMedium = medium_main, 
    redeclare package FlueGasMedium = medium_main, 
     
    // use Marchionni PCHE HeatTransfer
    // slow but can have a result - set a_phi = 0 to use Gnielinski's correlation 
    redeclare replaceable model HeatTransfer_F = Steps.TPComponents.MarchionniPCHEHeatTransferFV(),
    redeclare replaceable model HeatTransfer_G = Steps.TPComponents.MarchionniPCHEHeatTransferFV(),
    gasFlow(      
      heatTransfer(
        pitch = cfg_LTR.cfg_hot.l_pitch,
        phi   = cfg_LTR.cfg_hot.a_phi,
        Cf_C1 = Cf_C1,
        Cf_C2 = Cf_C2,
        Cf_C3 = Cf_C3
      )
    ),
    fluidFlow(      
      heatTransfer(
        pitch = cfg_LTR.cfg_cold.l_pitch,
        phi   = cfg_LTR.cfg_cold.a_phi,
        Cf_C1 = Cf_C1,
        Cf_C2 = Cf_C2,
        Cf_C3 = Cf_C3
      )
    ),    
    
    // fast and works fine for now. Error occurs when mass flow rate is zero, i.e. one flow is shut down. 
    // redeclare replaceable model HeatTransfer_F = ThermoPower.Thermal.HeatTransferFV.IdealHeatTransfer,      
    // redeclare replaceable model HeatTransfer_G = ThermoPower.Thermal.HeatTransferFV.IdealHeatTransfer, 
    cfg              = cfg_LTR,
    N_G              = N_seg_LTR,
    N_F              = N_seg_LTR,
    SSInit           = true,
    gasQuasiStatic   = false,
    fluidQuasiStatic = false,
    table_k_metalwall = cfg.table_k_LTR_wall
    // metalWall(L = L_wall, w_ch = W_ch, h_ch = H_ch, dx = T_wall),
    // metalQuasiStatic = true
    // override the values of Am and L of metaltubeFV
    // to make them agree with semi-circular tube of PCHE
    // ('final' modifier of Am in metalTubeFv was removed as well)
    //metalTube(WallRes=false, L = 1, rhomcm=200, Am = HE.metalVol / 1) 
  )
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
   
  Steps.TPComponents.PCHE HTR(
    redeclare package FluidMedium = medium_main, 
    redeclare package FlueGasMedium = medium_main, 
     
    // use Marchionni PCHE HeatTransfer
    // slow but can have a result - set a_phi = 0 to use Gnielinski's correlation 
    
    redeclare replaceable model HeatTransfer_F = Steps.TPComponents.MarchionniPCHEHeatTransferFV(),
    redeclare replaceable model HeatTransfer_G = Steps.TPComponents.MarchionniPCHEHeatTransferFV(),
    gasFlow(      
      heatTransfer(
        pitch = cfg_HTR.cfg_hot.l_pitch,
        phi   = cfg_HTR.cfg_hot.a_phi,
        Cf_C1 = Cf_C1,
        Cf_C2 = Cf_C2,
        Cf_C3 = Cf_C3
      )
    ),
    fluidFlow(      
      heatTransfer(
        pitch = cfg_HTR.cfg_cold.l_pitch,
        phi   = cfg_HTR.cfg_cold.a_phi,
        Cf_C1 = Cf_C1,
        Cf_C2 = Cf_C2,
        Cf_C3 = Cf_C3
      )
    ),    
    
    /*
    redeclare replaceable model HeatTransfer_F = ThermoPower.Thermal.HeatTransferFV.DittusBoelter,    
    redeclare replaceable model HeatTransfer_G = ThermoPower.Thermal.HeatTransferFV.DittusBoelter,
    gasFlow(heatTransfer(heating=true)),
    fluidFlow(heatTransfer(heating=false)),    
    */
    
    // fast and works fine for now. Error occurs when mass flow rate is zero, i.e. one flow is shut down. 
    // redeclare replaceable model HeatTransfer_F = ThermoPower.Thermal.HeatTransferFV.IdealHeatTransfer,      
    // redeclare replaceable model HeatTransfer_G = ThermoPower.Thermal.HeatTransferFV.IdealHeatTransfer,  
    cfg              = cfg_HTR,
    N_G              = N_seg_HTR,
    N_F              = N_seg_HTR,
    SSInit           = true,
    gasQuasiStatic   = false,
    fluidQuasiStatic = false,
    table_k_metalwall = cfg.table_k_HTR_wall  
    // metalWall(L = L_wall, w_ch = W_ch, h_ch = H_ch, dx = T_wall),
    // table_k_metalwall = table_k_metalwall
    // metalQuasiStatic = true
    // override the values of Am and L of metaltubeFV
    // to make them agree with semi-circular tube of PCHE
    // ('final' modifier of Am in metalTubeFv was removed as well)
    //metalTube(WallRes=false, L = 1, rhomcm=200, Am = HE.metalVol / 1) 
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
 
  ThermoPower.Gas.SensT sens_turbine(redeclare package Medium = medium_main) annotation(
    Placement(visible = true, transformation(origin = {20, 44}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));

/*
  // Input signals for transient simulation
  // hot/gas side
  Modelica.Blocks.Sources.IntegerConstant const_T_step_h(k = 20) annotation(
    Placement(visible = true, transformation(origin = {-230, 10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Interaction.Show.IntegerValue disp_T_h annotation(
    Placement(visible = true, transformation(origin = {-108, 44}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.IntegerConstant const_T_offset_h(k = 630) annotation(
    Placement(visible = true, transformation(origin = {-180, -54}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.MathInteger.TriggeredAdd triadd_T_h annotation(
    Placement(visible = true, transformation(origin = {-182, 8}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
  Modelica.Blocks.MathInteger.Sum sum_T_h(nu = 2) annotation(
    Placement(visible = true, transformation(origin = {-144, 8}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.BooleanPulse en_triadd_T_h(period = 10, width = 10) annotation(
    Placement(visible = true, transformation(origin = {-224, -24}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Interaction.Show.IntegerValue disp_T_step_h annotation(
    Placement(visible = true, transformation(origin = {-140, 42}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Math.IntegerToReal I2R_T_h annotation(
    Placement(visible = true, transformation(origin = {-102, 8}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
*/
  // cold/fluid side
  Modelica.Blocks.Sources.IntegerConstant const_T_offset_c(k = integer(st_source_cold.T)) annotation(
    Placement(visible = true, transformation(origin = {-128, 72}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Interaction.Show.IntegerValue disp_T_c annotation(
    Placement(visible = true, transformation(origin = {-144, 158}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.MathInteger.TriggeredAdd triadd_T_c annotation(
    Placement(visible = true, transformation(origin = {-128, 126}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
  Modelica.Blocks.Sources.BooleanPulse en_triadd_T_c(period = 10, startTime = 3, width = 10) annotation(
    Placement(visible = true, transformation(origin = {-156, 102}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.MathInteger.Sum sum_T_c(nu = 2) annotation(
    Placement(visible = true, transformation(origin = {-94, 126}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Interaction.Show.IntegerValue disp_T_step_c annotation(
    Placement(visible = true, transformation(origin = {-56, 158}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Math.IntegerToReal I2R_T_c annotation(
    Placement(visible = true, transformation(origin = {-48, 124}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.IntegerConstant const_T_step_c(k = T_step) annotation(
    Placement(visible = true, transformation(origin = {-194, 126}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));

protected
  parameter Real tablePhic[5, 4] = [1, 37, 80, 100; 1.5, 7.10E-05, 7.10E-05, 7.10E-05; 2, 8.40E-05, 8.40E-05, 8.40E-05; 2.5, 8.70E-05, 8.70E-05, 8.70E-05; 3, 1.04E-04, 1.04E-04, 1.04E-04];
  parameter Real tableEta[5, 4] = [1, 37, 80, 100; 1.5, 0.57, 0.89, 0.81; 2, 0.46, 0.82, 0.88; 2.5, 0.41, 0.76, 0.85; 3, 0.38, 0.72, 0.82];
  
equation
/*
  // HTR alone
  connect(source_cold.flange, HTR.waterIn);

  connect(HTR.waterOut, T_waterOut.inlet) annotation(
    Line(points = {{8.88178e-016, -44}, {8.88178e-016, -20}, {0, -20}}, thickness = 0.5, color = {0, 0, 255}));
  connect(T_waterOut.outlet, sink_cold.flange) annotation(
    Line(points = {{1.83697e-015, -70}, {1.83697e-015, -56}, {-8.88178e-016, -56}}, thickness = 0.5, color = {0, 0, 255}));

  // connect(HTR.waterOut, sink_cold.flange) 
   //annotation(Line(points = {{1.83697e-015, -70}, {1.83697e-015, -56}, {-8.88178e-016, -56}}, thickness = 0.5, color = {0, 0, 255}));   
  
  connect(source_hot.flange, HTR.gasIn) annotation(
    Line(points = {{-50, 0}, {-20, 0}}, color = {159, 159, 223}, thickness = 0.5, smooth = Smooth.None));   
    
  connect(HTR.gasOut, T_gasOut.inlet ) annotation(
    Line(points = {{34, 0}, {34, 0}, {20, 0}}, color = {159, 159, 223}, thickness = 0.5));
  connect(T_gasOut.outlet, sink_hot.flange) annotation(
    Line(points = {{46, 0}, {46, 0}, {60, 0}}, color = {159, 159, 223}, thickness = 0.5));  

  // LTR alone
  connect(source_cold.flange, LTR.waterIn);

  connect(LTR.waterOut, T_waterOut.inlet) annotation(
    Line(points = {{8.88178e-016, -44}, {8.88178e-016, -20}, {0, -20}}, thickness = 0.5, color = {0, 0, 255}));
  connect(T_waterOut.outlet, sink_cold.flange) annotation(
    Line(points = {{1.83697e-015, -70}, {1.83697e-015, -56}, {-8.88178e-016, -56}}, thickness = 0.5, color = {0, 0, 255}));

  // connect(HTR.waterOut, sink_cold.flange) 
   //annotation(Line(points = {{1.83697e-015, -70}, {1.83697e-015, -56}, {-8.88178e-016, -56}}, thickness = 0.5, color = {0, 0, 255}));   
  
  connect(source_hot.flange, LTR.gasIn) annotation(
    Line(points = {{-50, 0}, {-20, 0}}, color = {159, 159, 223}, thickness = 0.5, smooth = Smooth.None));   
    
  connect(LTR.gasOut, T_gasOut.inlet ) annotation(
    Line(points = {{34, 0}, {34, 0}, {20, 0}}, color = {159, 159, 223}, thickness = 0.5));
  connect(T_gasOut.outlet, sink_hot.flange) annotation(
    Line(points = {{46, 0}, {46, 0}, {60, 0}}, color = {159, 159, 223}, thickness = 0.5));  
*/
   
 
/*
  //HTR + mixer

  connect(source_mixer_in.flange, mixer.inlet1);  
  connect(source_cold.flange, mixer.inlet2);  
  connect(mixer.outlet, HTR.waterIn);

  connect(HTR.waterOut, T_waterOut.inlet) annotation(
    Line(points = {{8.88178e-016, -44}, {8.88178e-016, -20}, {0, -20}}, thickness = 0.5, color = {0, 0, 255}));
  connect(sink_cold.flange, T_waterOut.outlet) annotation(
    Line(points = {{1.83697e-015, -70}, {1.83697e-015, -56}, {-8.88178e-016, -56}}, thickness = 0.5, color = {0, 0, 255}));

  // connect(HTR.waterOut, sink_cold.flange) 
   //annotation(Line(points = {{1.83697e-015, -70}, {1.83697e-015, -56}, {-8.88178e-016, -56}}, thickness = 0.5, color = {0, 0, 255}));   
  
  connect(source_hot.flange, HTR.gasIn) annotation(
    Line(points = {{-50, 0}, {-20, 0}}, color = {159, 159, 223}, thickness = 0.5, smooth = Smooth.None));   
    
  connect(HTR.gasOut, T_gasOut.inlet ) annotation(
    Line(points = {{34, 0}, {34, 0}, {20, 0}}, color = {159, 159, 223}, thickness = 0.5));
  connect(T_gasOut.outlet, sink_hot.flange) annotation(
    Line(points = {{46, 0}, {46, 0}, {60, 0}}, color = {159, 159, 223}, thickness = 0.5));
  
  // connect(HTR.gasOut, sink_hot.flange) annotation(
  //   Line(points = {{46, 0}, {46, 0}, {60, 0}}, color = {159, 159, 223}, thickness = 0.5));
*/  
/*
  // mixer + LTR
  // water/cold side  
  connect(source_mixer_in.flange, mixer.inlet1);
  
  connect(source_cold.flange, LTR.waterIn);
  
  connect(LTR.waterOut, mixer.inlet2);
  
  connect(mixer.outlet, T_waterOut.inlet);
  
  connect(T_waterOut.outlet, sink_cold.flange) annotation(
    Line(points = {{1.83697e-015, -70}, {1.83697e-015, -56}, {-8.88178e-016, -56}}, thickness = 0.5, color = {0, 0, 255}));
  
  // gas/hot side
  connect(source_hot.flange, LTR.gasIn);
  
  connect(LTR.gasOut, T_gasOut.inlet);
  connect(T_gasOut.outlet, sink_hot.flange) annotation(
    Line(points = {{46, 0}, {46, 0}, {60, 0}}, color = {159, 159, 223}, thickness = 0.5));
*/

/*
  //HTR + mixer + LTR 
  // water/cold side    
  connect(source_cold.flange, r_LTR_cin.inlet);  
  connect(r_LTR_cin.outlet, LTR.waterIn);  
  connect(source_mixer_in.flange, mixer.inlet1);  
  connect(LTR.waterOut, r_LTR_cout.inlet);
  connect(r_LTR_cout.outlet, mixer.inlet2);  
  connect(mixer.outlet, r_HTR_cin.inlet);
  connect(r_HTR_cin.outlet, HTR.waterIn);
  connect(HTR.waterOut, r_HTR_cout.inlet);  
  connect(r_HTR_cout.outlet, sink_cold.flange) annotation(
    Line(points = {{1.83697e-015, -70}, {1.83697e-015, -56}, {-8.88178e-016, -56}}, thickness = 0.5, color = {0, 0, 255}));
  
  // gas/hot side
  connect(source_hot.flange, r_HTR_hin.inlet);
  connect(r_HTR_hin.outlet, HTR.gasIn) annotation(
   Line(points = {{-50, 0}, {-20, 0}}, color = {159, 159, 223}, thickness = 0.5, smooth = Smooth.None));
  connect(HTR.gasOut, r_HTR_hout.inlet);
  connect(r_HTR_hout.outlet, r_LTR_hin.inlet);
  connect(r_LTR_hin.outlet, LTR.gasIn);  
  connect(LTR.gasOut, r_LTR_hout.inlet);  
  connect(r_LTR_hout.outlet, sink_hot.flange) annotation(
    Line(points = {{46, 0}, {46, 0}, {60, 0}}, color = {159, 159, 223}, thickness = 0.5));
*/

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
  // temperature input
  // hot / gas side
  connect(const_T_step_h.y, triadd_T_h.u) annotation(
    Line(points = {{-219, 10}, {-202.5, 10}, {-202.5, 8}, {-190, 8}}, color = {255, 127, 0}));
  connect(const_T_offset_h.y, sum_T_h.u[2]) annotation(
    Line(points = {{-169, -54}, {-154, -54}, {-154, 8}}, color = {255, 127, 0}));
  connect(sum_T_h.y, disp_T_h.numberPort) annotation(
    Line(points = {{-132.5, 8}, {-121, 8}, {-121, 44}, {-119.5, 44}}, color = {255, 127, 0}));
  connect(triadd_T_h.y, disp_T_step_h.numberPort) annotation(
    Line(points = {{-175, 8}, {-160.5, 8}, {-160.5, 42}, {-151.5, 42}}, color = {255, 127, 0}));
  connect(triadd_T_h.y, sum_T_h.u[1]) annotation(
    Line(points = {{-175, 8}, {-154, 8}}, color = {255, 127, 0}));
  connect(en_triadd_T_h.y, triadd_T_h.trigger) annotation(
    Line(points = {{-213, -24}, {-186, -24}, {-186, 1}}, color = {255, 0, 255}));
  connect(sum_T_h.y, I2R_T_h.u) annotation(
    Line(points = {{-132.5, 8}, {-114, 8}}, color = {255, 127, 0}));
  connect(I2R_T_h.y, source_hot.in_T) annotation(
    Line(points = {{-91, 8}, {-90, 8}, {-90, 9}, {-80, 9}, {-80, 7.75}, {-78, 7.75}, {-78, 6}}, color = {0, 0, 127}));
*/

  // cold / fluid side 
  connect(en_triadd_T_c.y, triadd_T_c.trigger) annotation(
    Line(points = {{-145, 102}, {-132, 102}, {-132, 119}}, color = {255, 0, 255}));
  connect(triadd_T_c.y, sum_T_c.u[1]) annotation(
    Line(points = {{-121, 126}, {-104, 126}}, color = {255, 127, 0}));
  connect(triadd_T_c.y, disp_T_c.numberPort) annotation(
    Line(points = {{-121, 126}, {-160.5, 126}, {-160.5, 158}, {-155.5, 158}}, color = {255, 127, 0}));
  connect(sum_T_c.y, I2R_T_c.u) annotation(
    Line(points = {{-82.5, 126}, {-71.25, 126}, {-71.25, 124}, {-60, 124}}, color = {255, 127, 0}));
  connect(sum_T_c.y, disp_T_step_c.numberPort) annotation(
    Line(points = {{-82.5, 126}, {-71, 126}, {-71, 158}, {-67.5, 158}}, color = {255, 127, 0}));
  connect(const_T_step_c.y, triadd_T_c.u) annotation(
    Line(points = {{-183, 126}, {-136, 126}}, color = {255, 127, 0}));
  connect(const_T_offset_c.y, sum_T_c.u[2]) annotation(
    Line(points = {{-117, 72}, {-104, 72}, {-104, 126}}, color = {255, 127, 0}));
  connect(I2R_T_c.y, source_cold.in_T) annotation(
    Line(points = {{-37, 124}, {6, 124}, {6, 64}}, color = {0, 0, 127}));


annotation(
    Diagram(graphics),
    experiment(StartTime = 0, StopTime = 100, Tolerance = 1e-3, Interval = 5),    
    // options = "-showErrorMessages -demoMode",
    __OpenModelica_commandLineOptions = "--matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection -d=initialization,NLSanalyticJacobian,aliasConflicts,dumpCSE",    
    __OpenModelica_simulationFlags(lv = "LOG_DEBUG,LOG_NLS,LOG_NLS_V,LOG_STATS,LOG_INIT,LOG_STDOUT, -w", outputFormat = "mat", s = "dassl", nls = "homotopy"));
end TestDyn_Comps_PCHE_ramp;
