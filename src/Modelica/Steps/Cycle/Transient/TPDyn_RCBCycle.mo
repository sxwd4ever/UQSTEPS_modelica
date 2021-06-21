within Steps.Cycle.Transient;

model TPDyn_RCBCycle
  "Recompressoion Brayton Cycle transient simulation based on ThermoPower"
  
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
    // substanceNames = {"CO2|debug=40"}
  );
  
  // package medium_main = ExternalMedia.Examples.CO2CoolProp;
  // package medium_heater = Steps.Media.ThermiaOilD; // out of working range of this 10Mw high T loop
  package medium_heater = SolarTherm.Media.Sodium.Sodium_pT;
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
  // parameter Real Cf_C1                  = 1.626, Cf_C2 = 1, Cf_C3 = 1;
  parameter Real Cf_C1                  = 1, Cf_C2 = 1, Cf_C3 = 1;
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
 
  // select the configuration of parameters
  parameter Model.RCBCycleConfig cfg(
    redeclare package medium_heater = medium_heater,
    redeclare package medium_main   = medium_main,
    redeclare package medium_cooler = medium_cooler,
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
    // table_k_LTR_wall = table_k_metalwall,
    // table_k_HTR_wall = table_k_metalwall,
    
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

    mdot_main   = 128.774, // mdot on 10WMe point
    mdot_comp   = 88.0661, // mdot on 10WMe point 
        
    // mdot_main   = 109.458, // minimum mdot for serial ramp case III in guan's report.    
    // mdot_comp   = 74.856,  // minimum mdot at above point    
    
    mdot_heater = 40,
    mdot_cooler = 40.7188
  );
  
  // set the values of parameters accordingly
  parameter Model.HeatExchangerConfig cfg_heater  = cfg.cfg_heater;
  parameter Model.HeatExchangerConfig cfg_LTR     = cfg.cfg_LTR;
  parameter Model.HeatExchangerConfig cfg_HTR     = cfg.cfg_HTR;
  parameter Model.TurbomachineryConfig cfg_turb   = cfg.cfg_turb;
  parameter Model.TurbomachineryConfig cfg_comp   = cfg.cfg_comp;
  parameter Model.TurbomachineryConfig cfg_recomp = cfg.cfg_recomp;
  parameter Model.HeatExchangerConfig cfg_cooler  = cfg.cfg_cooler;
  parameter Model.SplitterConfig cfg_splitter     = cfg.cfg_splitter;
  parameter Model.SplitterConfig cfg_merger       = cfg.cfg_merger;
  
  parameter Model.ThermoState st_bypass      = cfg_recomp.st_in;
  parameter Model.ThermoState st_source_temp = cfg_cooler.cfg_hot.st_out; //.cfg_hot.st_out;
  parameter Model.ThermoState st_sink_temp   = cfg_cooler.cfg_hot.st_out; //.cfg_hot.st_out;
  parameter Model.ThermoState st_source      = cfg_heater.cfg_cold.st_out;
  parameter Model.ThermoState st_sink        = cfg_heater.cfg_cold.st_out;

  parameter Integer N_seg_heater = cfg.cfg_heater.cfg_hot.geo_path.N_seg; 
  parameter Integer N_seg_LTR    = cfg.cfg_LTR.cfg_hot.geo_path.N_seg; 
  parameter Integer N_seg_HTR    = cfg.cfg_HTR.cfg_hot.geo_path.N_seg;   
  parameter Integer N_seg_cooler = cfg.cfg_cooler.cfg_hot.geo_path.N_seg; 

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
    Placement(visible = true, transformation(origin = {-74, 4}, extent = {{-4, -4}, {4, 4}}, rotation = 90)));
    
  ThermoPower.Water.SinkPressure sink_heater_hot(
    redeclare package Medium = medium_heater, 
    p0    = cfg_heater.cfg_hot.st_out.p,
    T     = cfg_heater.cfg_hot.st_out.T,
    use_T = false) 
  annotation(
    Placement(visible = true, transformation(origin = {-74, 62}, extent = {{-4, -4}, {4, 4}}, rotation = 90)));
    
  Steps.TPComponents.HE heater(
    redeclare package FluidMedium = medium_heater, 
    redeclare package FlueGasMedium = medium_main, 
    
    // Ngo 
    redeclare replaceable model HeatTransfer_F = Steps.TPComponents.NgoHeatTransferFV(),
    redeclare replaceable model HeatTransfer_G = Steps.TPComponents.NgoHeatTransferFV(),
    gasFlow(      
      heatTransfer(
        Cf_C1     = Cf_C1,
        Cf_C2     = Cf_C2,
        Cf_C3     = Cf_C3,
        gamma_min = 1000,
        gamma_max = 5000
      ),
      Tstartin  = cfg_heater.cfg_gas.st_in.T,
      Tstartout = cfg_heater.cfg_gas.st_out.T,
      Nt        = cfg_heater.cfg_gas.N_ch,
      Dhyd      = cfg_heater.cfg_gas.geo_area.d_hyd
    ),
    fluidFlow(      
      heatTransfer(
        Cf_C1     = Cf_C1,
        Cf_C2     = Cf_C2,
        Cf_C3     = Cf_C3,
        gamma_min = 2000,
        gamma_max = 1.5e6
      ), 
      fixedMassFlowSimplified = true,
      hstartin                = cfg_heater.cfg_fluid.st_in.h,
      hstartout               = cfg_heater.cfg_fluid.st_out.h
    ),    
    redeclare model HeatExchangerTopology = ThermoPower.Thermal.HeatExchangerTopologies.CounterCurrentFlow,
    cfg             = cfg_heater,
    N_G             = N_seg_heater,
    N_F             = N_seg_heater,
    SSInit          = true,
    gasQuasiStatic  = false,
    FluidPhaseStart = ThermoPower.Choices.FluidPhase.FluidPhases.Liquid
  )
  annotation(
    Placement(visible = true, transformation(origin = {-74, 32}, extent = {{-8, -8}, {8, 8}}, rotation = 180)));  
  
  Steps.TPComponents.WaterStateReader r_heater_hin(redeclare package Medium = medium_heater) 
  annotation(
    Placement(visible = true, transformation(origin = {-74, 16}, extent = {{-4, -4}, {4, 4}}, rotation = 90)));  
  Steps.TPComponents.WaterStateReader r_heater_hout(redeclare package Medium = medium_heater)
  annotation(
    Placement(visible = true, transformation(origin = {-74, 50}, extent = {{-4, -4}, {4, 4}}, rotation = 90)));
  Steps.TPComponents.GasStateReader r_heater_cin(redeclare package Medium = medium_main) 
  annotation(
    Placement(visible = true, transformation(origin = {-56, 32}, extent = {{-4, -4}, {4, 4}}, rotation = 180)));
  Steps.TPComponents.GasStateReader r_heater_cout(redeclare package Medium = medium_main) 
  annotation(
    Placement(visible = true, transformation(origin = {-90, 32}, extent = {{-4, -4}, {4, 4}}, rotation = 180)));
/*    
  // Case I - sourceMassFlow + sinkPressure
  ThermoPower.Gas.SourceMassFlow source(
    redeclare package Medium = medium_main, 
    T         = st_source.T,
    p0        = st_source.p,
    w0        = st_source.mdot,
    use_in_T  = false,
    use_in_w0 = true,
    gas(
      p(start = st_source.p, nominal = st_source.p), 
      T(start = st_source.T, nominal = st_source.T)))
      //h(start = st_source.h, nominal = st_source.h)))
  annotation(
    Placement(visible = true, transformation(origin = {-107, -7}, extent = {{-5, -5}, {5, 5}}, rotation = 0)));    
*/

  // Case II - sourcePressure + sinkMassFlow
  ThermoPower.Gas.SourcePressure source(
    redeclare package Medium = medium_main, 
    T        = st_source.T,
    p0       = st_source.p,
    use_in_T = false,
    use_in_p0 = true,
    gas(
      p(start = st_source.p, nominal = st_source.p), 
      T(start = st_source.T, nominal = st_source.T))) 
  annotation(
    Placement(visible = true, transformation(origin = {-107, -7}, extent = {{-5, -5}, {5, 5}}, rotation = 0)));


/*
  ThermoPower.Gas.SourceMassFlow source_mixer_in(
    redeclare package Medium = medium_main,
    T        = st_bypass.T,
    p0       = st_bypass.p,
    use_in_T = false,
    use_in_w0 = false,
    w0       = st_bypass.mdot,    
    gas(
      p(start = st_bypass.p, nominal = st_bypass.p), 
      T(start = st_bypass.T, nominal = st_bypass.T),
      h(start = st_bypass.h, nominal = st_bypass.h))
  );  

  ThermoPower.Gas.SinkPressure sink_mixer_in(
    redeclare package Medium = medium_main, 
    p0 = st_bypass.p,
    T  = st_bypass.T)
  annotation(
    Placement(transformation(origin = {0, -80}, extent = {{-10, -10}, {10, 10}}, rotation = 270)));
*/
/* 
  // Case I - sourceMassFlow + sinkPressure
  ThermoPower.Gas.SinkPressure sink(
    redeclare package Medium = medium_main, 
    p0        = st_sink.p,
    T         = st_sink.T,
    use_in_p0 = true,
    gas(
      p(start = st_sink.p, nominal = st_sink.p), 
      T(start = st_sink.T, nominal = st_sink.T)))
  annotation(
    Placement(visible = true, transformation(origin = {-105, 33}, extent = {{-5, 5}, {5, -5}}, rotation = 180)));
*/   

  // Case II - sourcePressure + sinkMassFlow
  ThermoPower.Gas.SinkMassFlow sink(
    redeclare package Medium = medium_main, 
    T         = st_sink.T,
    p0        = st_sink.p,
    use_in_T  = false,
    use_in_w0 = true,
    w0        = st_sink.mdot
  ) 
  annotation(
    Placement(visible = true, transformation(origin = {-105, 33}, extent = {{-5, 5}, {5, -5}}, rotation = 180)));


  // use FlowJoin to mix flow
  Gas.FlowJoin mixer(redeclare package Medium = medium_main) 
  annotation(
    Placement(visible = true, transformation(origin = {-4, 54}, extent = {{-4, -4}, {4, 4}}, rotation = 180))); 

  Steps.TPComponents.PCHE LTR(
    redeclare package FluidMedium = medium_main, 
    redeclare package FlueGasMedium = medium_main, 
     
    // use Ngo PCHE HeatTransfer
    // slow but can have a result - set a_phi = 0 to use Gnielinski's correlation 
    redeclare replaceable model HeatTransfer_F = Steps.TPComponents.NgoHeatTransferFV(),
    redeclare replaceable model HeatTransfer_G = Steps.TPComponents.NgoHeatTransferFV(),
    gasFlow(      
      heatTransfer(
        pitch     = cfg_LTR.cfg_hot.l_pitch,
        phi       = cfg_LTR.cfg_hot.a_phi,
        Cf_C1     = Cf_C1,
        Cf_C2     = Cf_C2,
        Cf_C3     = Cf_C3,
        gamma_min = 2000,
        gamma_max = 5000
      )
    ),
    fluidFlow(      
      heatTransfer(
        pitch     = cfg_LTR.cfg_cold.l_pitch,
        phi       = cfg_LTR.cfg_cold.a_phi,
        Cf_C1     = Cf_C1,
        Cf_C2     = Cf_C2,
        Cf_C3     = Cf_C3,
        gamma_min = 2000,
        gamma_max = 5000
      )
    ),  
    cfg              = cfg_LTR,
    N_G              = N_seg_LTR,
    N_F              = N_seg_LTR,
    SSInit           = true,
    gasQuasiStatic   = false,
    fluidQuasiStatic = false
  )
  annotation(
    Placement(visible = true, transformation(origin = {20, 32}, extent = {{-8, -8}, {8, 8}}, rotation = 0)));
        
  Steps.TPComponents.GasStateReader r_LTR_hin(redeclare package Medium = medium_main)
  annotation(
    Placement(visible = true, transformation(origin = {6, 32}, extent = {{-4, -4}, {4, 4}}, rotation = 0)));
  Steps.TPComponents.GasStateReader r_LTR_hout(redeclare package Medium = medium_main) 
  annotation(
    Placement(visible = true, transformation(origin = {34, 32}, extent = {{-4, -4}, {4, 4}}, rotation = 0)));    
  Steps.TPComponents.GasStateReader r_LTR_cin(redeclare package Medium = medium_main)
  annotation(
    Placement(visible = true, transformation(origin = {20, 46}, extent = {{-4, -4}, {4, 4}}, rotation = -90)));    
  Steps.TPComponents.GasStateReader r_LTR_cout(redeclare package Medium = medium_main) 
  annotation(
    Placement(visible = true, transformation(origin = {20, 16}, extent = {{-4, -4}, {4, 4}}, rotation = -90)));
    
  Steps.TPComponents.PCHE HTR(
    redeclare package FluidMedium = medium_main, 
    redeclare package FlueGasMedium = medium_main, 
     
    // use Gnielinski PCHE HeatTransfer
    // slow but can have a result - set a_phi = 0 to use Gnielinski's correlation     
    redeclare replaceable model HeatTransfer_F = Steps.TPComponents.GnielinskiHeatTransferFV(),
    redeclare replaceable model HeatTransfer_G = Steps.TPComponents.GnielinskiHeatTransferFV(),
    gasFlow(      
      heatTransfer(
        pitch     = cfg_HTR.cfg_gas.l_pitch,
        phi       = cfg_HTR.cfg_gas.a_phi,
        Cf_C1     = Cf_C1,
        Cf_C2     = Cf_C2,
        Cf_C3     = Cf_C3,
        gamma_min = 2000,
        gamma_max = 5000
      )
    ),
    fluidFlow(      
      heatTransfer(
        pitch     = cfg_HTR.cfg_fluid.l_pitch,
        phi       = cfg_HTR.cfg_fluid.a_phi,
        Cf_C1     = Cf_C1,
        Cf_C2     = Cf_C2,
        Cf_C3     = Cf_C3,
        gamma_min = 2000,
        gamma_max = 5000
      )
    ), 
    cfg              = cfg_HTR,
    N_G              = N_seg_HTR,
    N_F              = N_seg_HTR,
    SSInit           = true,
    gasQuasiStatic   = false,
    fluidQuasiStatic = false
     
  )
  annotation(
    Placement(visible = true, transformation(origin = {-22, 32}, extent = {{-8, -8}, {8, 8}}, rotation = 0)));  
    
  Steps.TPComponents.GasStateReader r_HTR_hin(redeclare package Medium = medium_main)
  annotation(
    Placement(visible = true, transformation(origin = {-36, 32}, extent = {{-4, -4}, {4, 4}}, rotation = 0)));    
  Steps.TPComponents.GasStateReader r_HTR_hout(redeclare package Medium = medium_main)
  annotation(
    Placement(visible = true, transformation(origin = {-8, 32}, extent = {{-4, -4}, {4, 4}}, rotation = 0)));
  Steps.TPComponents.GasStateReader r_HTR_cin(redeclare package Medium = medium_main) 
  annotation(
    Placement(visible = true, transformation(origin = {-22, 46}, extent = {{-4, -4}, {4, 4}}, rotation = -90)));
  Steps.TPComponents.GasStateReader r_HTR_cout(
    redeclare package Medium = medium_main,
    gas(
      p(start = cfg_HTR.cfg_cold.st_out.p, nominal = cfg_HTR.cfg_cold.st_out.p),
      h(start = cfg_HTR.cfg_cold.st_out.h, nominal = cfg_HTR.cfg_cold.st_out.h))) 
  annotation(
    Placement(visible = true, transformation(origin = {-22, 16}, extent = {{-4, -4}, {4, 4}}, rotation = -90)));
    
  ThermoPower.Gas.Turbine turbine(
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
      p(nominal = turbine.pstart_in), 
      T(nominal = turbine.Tstart_in),
      h(start   = cfg_turb.st_out.h, nominal=cfg_turb.st_out.h)),
    gas_iso(
      p(nominal = turbine.pstart_out), 
      T(nominal = turbine.Tstart_out),
      h(start   = cfg_turb.st_out.h)))    
  annotation(
    Placement(visible = true, transformation(origin = {-54, -40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    
  Steps.TPComponents.GasStateReader r_turb_in(redeclare package Medium = medium_main) 
  annotation(
    Placement(visible = true, transformation(origin = {-98, -28}, extent = {{-4, -4}, {4, 4}}, rotation = 0)));
  Steps.TPComponents.GasStateReader r_turb_out(redeclare package Medium = medium_main) 
  annotation(
    Placement(visible = true, transformation(origin = {-42, -30}, extent = {{-4, -4}, {4, 4}}, rotation = 90)));
    
  Modelica.Mechanics.Rotational.Sources.ConstantSpeed const_speed_turb(
      w_fixed=cfg_turb.N, useSupport=false) 
  annotation(
    Placement(visible = true, transformation(origin = {-73, -41}, extent = {{-5, -5}, {5, 5}}, rotation = 0)));
 
  // controller FPC (and the speed_turb) only required for Case I - ramp source mdot. Commnent, remove it otherwise. 
  Steps.TPComponents.TurbineFixedPController FPC_turb(
    redeclare package Medium   = medium_main,
    fileName  = Modelica.Utilities.Files.loadResource("modelica://Steps/Resources/Data/turbine_map.txt"),
    tablePhic = fill(0.0, 14, 12), 
    Table     = ThermoPower.Choices.TurboMachinery.TableTypes.file,
    Ndesign   = cfg_turb.N,
    Tdes_in   = cfg_turb.st_in.T,
    use_in_p2 = true,    
    in_T1(start = cfg_turb.st_in.T, nominal = cfg_turb.st_in.T),
    in_p1(start = cfg_turb.st_in.p, nominal = cfg_turb.st_in.p),
    in_w1(start = cfg_turb.st_in.mdot, nominal = cfg_turb.st_in.mdot)
  ) "Fixed pressure controller for main compressor"
  annotation(
    Placement(visible = false, transformation(origin = {-90, -41}, extent = {{-5, -5}, {5, 5}}, rotation = 0)));  

  Modelica.Mechanics.Rotational.Sources.Speed speed_turb(
    exact = false,
    w_ref(start = cfg_turb.N, nominal= cfg_turb.N),
    w(start = cfg_turb.N, nominal = cfg_turb.N)
  )
  annotation(
    Placement(visible = false, transformation(origin = {-73, -41}, extent = {{-5, -5}, {5, 5}}, rotation = 0)));   

  ThermoPower.Water.SourceMassFlow source_cooler(
    redeclare package Medium = medium_cooler,
    w0 = cfg_cooler.cfg_cold.st_in.mdot,
    p0 = cfg_cooler.cfg_cold.st_in.p,
    h  = cfg_cooler.cfg_cold.st_in.h,
    T  = cfg_cooler.cfg_cold.st_in.T) 
  annotation(
    Placement(visible = true, transformation(origin = {56, -38}, extent = {{-4, -4}, {4, 4}}, rotation = 270)));    
    
  ThermoPower.Water.SinkPressure sink_cooler(
    redeclare package Medium = medium_cooler, 
    p0    = cfg_cooler.cfg_cold.st_out.p,
    T     = cfg_cooler.cfg_cold.st_out.T,
    h     = cfg_cooler.cfg_cold.st_out.h,
    use_T = false) 
  annotation(
    Placement(visible = true, transformation(origin = {57, -97}, extent = {{-4, -4}, {4, 4}}, rotation = 270)));
    
  Steps.TPComponents.HE cooler(
    redeclare package FluidMedium = medium_cooler, 
    redeclare package FlueGasMedium = medium_main, 
       
    // Ngo 
    redeclare replaceable model HeatTransfer_F = Steps.TPComponents.NgoHeatTransferFV(),
    redeclare replaceable model HeatTransfer_G = Steps.TPComponents.NgoHeatTransferFV(),
    gasFlow(      
      heatTransfer(
        Cf_C1     = Cf_C1,
        Cf_C2     = Cf_C2,
        Cf_C3     = Cf_C3,
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
        Cf_C1     = Cf_C1,
        Cf_C2     = Cf_C2,
        Cf_C3     = Cf_C3,
        gamma_min = 2000,
        gamma_max = 18000
      ), 
      fixedMassFlowSimplified = true,
      hstartin                = cfg_cooler.cfg_fluid.st_in.h,
      hstartout               = cfg_cooler.cfg_fluid.st_out.h
    ),    
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
    Placement(visible = true, transformation(origin = {56, -66}, extent = {{-8, -8}, {8, 8}}, rotation = 0)));  
    
  Steps.TPComponents.GasStateReader r_cooler_hin(redeclare package Medium = medium_main) 
  annotation(
    Placement(visible = true, transformation(origin = {39, -67}, extent = {{-5, -5}, {5, 5}}, rotation = 0)));     
  Steps.TPComponents.GasStateReader r_cooler_hout(redeclare package Medium = medium_main)
  annotation(
    Placement(visible = true, transformation(origin = {75, -67}, extent = {{-5, -5}, {5, 5}}, rotation = 0)));  
  Steps.TPComponents.WaterStateReader r_cooler_cin(redeclare package Medium = medium_cooler) 
  annotation(
    Placement(visible = true, transformation(origin = {56, -50}, extent = {{-4, -4}, {4, 4}}, rotation = -90))); 
  Steps.TPComponents.WaterStateReader r_cooler_cout(redeclare package Medium = medium_cooler) 
  annotation(
    Placement(visible = true, transformation(origin = {56, -84}, extent = {{-4, -4}, {4, 4}}, rotation = -90)));
    
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
    gas_in(
      p(start   = cfg_comp.st_in.p, nominal = cfg_comp.st_in.p), 
      T(start   = cfg_comp.st_in.T, nominal = cfg_comp.st_in.T)),      
    gas_iso(
      p(nominal = cfg_comp.st_in.p), 
      T(nominal = cfg_comp.st_in.T),
      h(start   = cfg_comp.st_out.h)))
  annotation(
    Placement(visible = true, transformation(origin = {162, -4}, extent = {{-8, -8}, {8, 8}}, rotation = 0)));  
    
  Modelica.Mechanics.Rotational.Sources.ConstantSpeed const_speed_mc(
    w_fixed    = cfg_comp.N,
    useSupport = false)
  annotation(
    Placement(visible = false, transformation(origin = {133, 7}, extent = {{5, -5}, {-5, 5}}, rotation = 0)));
  
  Steps.TPComponents.GasStateReader r_comp_in(redeclare package Medium = medium_main) 
  annotation(
    Placement(visible = true, transformation(origin = {117, 19}, extent = {{-5, -5}, {5, 5}}, rotation = 0)));
  Steps.TPComponents.GasStateReader r_comp_out(redeclare package Medium = medium_main) 
  annotation(
    Placement(visible = true, transformation(origin = {168, 24}, extent = {{-4, -4}, {4, 4}}, rotation = 90))); 
     
  Steps.TPComponents.CompressorFixedPController FPC_mc(
    redeclare package Medium   = medium_main,
    tablePhic = tablePhic_comp_mc,
    tablePR   = tablePR_comp_mc,
    Table     = ThermoPower.Choices.TurboMachinery.TableTypes.matrix,
    Ndesign   = cfg_comp.N,
    Tdes_in   = cfg_comp.st_in.T,
    use_in_p2 = true,
    in_T1(start = cfg_comp.st_in.T, nominal = cfg_comp.st_in.T),
    in_p1(start = cfg_comp.st_in.p, nominal = cfg_comp.st_in.p),
    in_w1(start = cfg_comp.st_in.mdot, nominal = cfg_comp.st_in.mdot)
  ) "Fixed pressure controller for main compressor" 
  annotation(
    Placement(visible = true, transformation(origin = {129, -3}, extent = {{-5, -5}, {5, 5}}, rotation = 0)));
  
  Modelica.Mechanics.Rotational.Sources.Speed speed_mc(
    exact = false,
    w_ref(start = cfg_comp.N, nominal= cfg_comp.N),
    w(start = cfg_comp.N, nominal = cfg_comp.N)
  )
  annotation(
    Placement(visible = true, transformation(origin = {145, -3}, extent = {{-5, -5}, {5, 5}}, rotation = 0)));   
 
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
    gas_in(
      p(start = cfg_recomp.st_in.p, nominal = cfg_recomp.st_in.p), 
      T(start = cfg_recomp.st_in.T, nominal = cfg_recomp.st_in.T)),    
    gas_iso(
      p(nominal = cfg_recomp.st_in.p), 
      T(nominal = cfg_recomp.st_in.T),
      h(start   = cfg_recomp.st_out.h)))
  annotation(
    Placement(visible = true, transformation(origin = {97, -1}, extent = {{-9, -9}, {9, 9}}, rotation = 0)));    
         
  Modelica.Mechanics.Rotational.Sources.ConstantSpeed const_speed_rc(
    w_fixed    = cfg_recomp.N,
    useSupport = false)
  annotation(
    Placement(visible = false, transformation(origin = {91, 7}, extent = {{5, -5}, {-5, 5}}, rotation = 0)));  
    
  Steps.TPComponents.GasStateReader r_recomp_in(redeclare package Medium = medium_main) 
  annotation(
    Placement(visible = true, transformation(origin = {47, 19}, extent = {{5, 5}, {-5, -5}}, rotation = 180)));
  Steps.TPComponents.GasStateReader r_recomp_out(redeclare package Medium = medium_main) 
  annotation(
    Placement(visible = true, transformation(origin = {104, 24}, extent = {{-4, -4}, {4, 4}}, rotation = 90)));
    
  Steps.TPComponents.CompressorFixedPController FPC_rc(
    redeclare package Medium   = medium_main,
    tablePhic = tablePhic_comp_rc,
    tablePR   = tablePR_comp_rc,
    Table     = ThermoPower.Choices.TurboMachinery.TableTypes.matrix,
    Ndesign   = cfg_recomp.N,
    Tdes_in   = cfg_recomp.st_in.T,
    use_in_p2 = true,
    in_T1(start = cfg_recomp.st_in.T, nominal = cfg_recomp.st_in.T),
    in_p1(start = cfg_recomp.st_in.p, nominal = cfg_recomp.st_in.p),
    in_w1(start = cfg_recomp.st_in.mdot, nominal = cfg_recomp.st_in.mdot)
  ) "Fixed pressure controller for recompressor"
  annotation(
    Placement(visible = true, transformation(origin = {65, -3}, extent = {{-5, -5}, {5, 5}}, rotation = 0)));
  
  Modelica.Mechanics.Rotational.Sources.Speed speed_rc(
    exact = false,
    w_ref(start = cfg_recomp.N, nominal= cfg_recomp.N),
    w(start = cfg_recomp.N, nominal = cfg_recomp.N)
  )
  annotation(
    Placement(visible = true, transformation(origin = {81, -3}, extent = {{-5, -5}, {5, 5}}, rotation = 0))); 

  ThermoPower.Gas.FlowSplit splitter(redeclare package Medium = medium_main) 
  annotation(
    Placement(visible = true, transformation(origin = {37, 21}, extent = {{-5, 5}, {5, -5}}, rotation = -90)));
    
  // value to accelerate the simulation
  constant Integer SEC_PER_HOUR = integer(60 * 60 / time_scaling_factor); 
  // constant Real time_scaling_factor = 12; // 5 min = 1 hour
  constant Real time_scaling_factor = 1; // 1 hour = 1 hour
  
  // boundary conditions
 /*  
  // ramp input to simulate the ramp change in ASPEN+ simulation
  // ramp change for case I in ASPEN ASTRI report
  Steps.TPComponents.RampSeq_W ramp_mdot(    
    time_start = 3 * SEC_PER_HOUR,
    interval   = 1 * SEC_PER_HOUR,
    duration   = 0.15 * SEC_PER_HOUR,
    offset     = st_source.mdot
  )
  annotation(
    Placement(visible = true, transformation(origin = {-146, 91}, extent = {{-5, -5}, {5, 5}}, rotation = 0)));

  Modelica.Blocks.Sources.Constant p_const_TIP(k(start = cfg_turb.st_in.p, nominal = cfg_turb.st_in.p))
  annotation(
    Placement(visible = true, transformation(origin = {-146, 63}, extent = {{-5, -5}, {5, 5}}, rotation = 0)));  
  Modelica.Blocks.Sources.Constant p_const_CIP(k(start = cfg_comp.st_in.p, nominal = cfg_comp.st_in.p))
  annotation(
    Placement(visible = true, transformation(origin = {-146, 79}, extent = {{-5, -5}, {5, 5}}, rotation = 0)));
   */ 

 
  // ramp input to simulate the ramp change in ASPEN+ simulation
  // ramp change for case II in the ASPEN ASTRI report
  Steps.TPComponents.RampSeq_W ramp_p(    
    // time_start = 0.01 * SEC_PER_HOUR, 
    // interval = 1 * SEC_PER_HOUR, 
    time_start = 3 * SEC_PER_HOUR,
    interval   = 2 * SEC_PER_HOUR,
    duration   = 0.15 * SEC_PER_HOUR,
      
    offset  = st_source.p,
    ratio_1 = 0.05,
    ratio_2 = 0.05
  )  
  annotation(
    Placement(visible = true, transformation(origin = {-146, 79}, extent = {{-5, -5}, {5, 5}}, rotation = 0)));


  // ramp input for case III in ASPEN ASTRI report 
  /*
  Modelica.Blocks.Sources.Ramp ramp_mdot(
    final duration  = 0.15 * SEC_PER_HOUR,
    final startTime = 0.01 * SEC_PER_HOUR, // start after 36s (for test), should be at 1H. 
    
    final height    = -st_sink.mdot * 0.15 / time_scaling_factor,
    final offset    = st_sink.mdot
  ) 
  annotation(
    Placement(visible = true, transformation(origin = {-119, 21}, extent = {{-5, -5}, {5, 5}}, rotation = 0)));

  Steps.TPComponents.RampSeq_TwoStages ramp_p(
    final time_start = 2 * SEC_PER_HOUR,
    final interval   = 1 * SEC_PER_HOUR,
    final duration_1 = 0.15 * SEC_PER_HOUR,
    final duration_2 = 0.15 * SEC_PER_HOUR,
    
    final height_1 = -1e6 / time_scaling_factor,
    final height_2 = -1e6 / time_scaling_factor,
    final offset   = st_source.p
  )
  annotation(
    Placement(visible = true, transformation(origin = {-120, -32}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
*/
  
  /*
  // ramp component for transient test
  Modelica.Blocks.Sources.Ramp ramp_mdot(
  final height    = st_source.mdot * 0.1,
  final duration  = 9 * 60, // 9 mins
  final startTime = 3,
  final offset    = st_source.mdot);

  // ramp component for transient test
  Modelica.Blocks.Sources.Ramp ramp_pressure(
  final height    = st_sink.p * 0.1,
  final duration  = 9 * 60, // 9 mins
  final startTime = 3,
  final offset    = st_sink.p);
  */    
    
  // calculated variables
  
  // for turbine
  Modelica.SIunits.Power      W_turb   = (r_turb_in.h - r_turb_out.h) * r_turb_in.w / 1e6 "W->MW, net power for turbine";
  Modelica.SIunits.Efficiency eta_turb = turbine.eta * 100;
  // main compressor
  Modelica.SIunits.Power      W_MC   = (r_comp_out.h - r_comp_in.h) * r_comp_in.w / 1e6 "W->MW, net power for main compressor";
  Modelica.SIunits.Efficiency eta_MC = compressor.eta * 100;
  // re compressor
  Modelica.SIunits.Power      W_RC   = (r_recomp_out.h - r_recomp_in.h) * r_recomp_in.w / 1e6 "W->MW, net power for recompressor";
  Modelica.SIunits.Efficiency eta_RC = recompressor.eta * 100;
  
  // power block
  Modelica.SIunits.Power      W_net    = W_turb - W_MC - W_RC "net power generated";
  Modelica.SIunits.Power      Q_heater = (r_heater_hin.h - r_heater_hout.h) * r_heater_hin.w / 1e6 "heat input for the power block";
  Modelica.SIunits.Power      Q_cooler = (r_cooler_cout.h - r_cooler_cin.h) * r_cooler_cin.w / 1e6 "heat input for the power block";
  Modelica.SIunits.Efficiency eta_pb   = W_net / Q_heater * 100 "power block efficiency";
  Real                        SR       = r_recomp_in.w / r_comp_in.w * 100 "split ratio";

  // heat transfer coefficient for HTR and LTR
  // HTR
  Modelica.SIunits.Power                 Q_HTR      = (r_HTR_cout.h - r_HTR_cin.h) * r_HTR_cin.w / 1e6 "W->MW, heat input for HTR";
  Modelica.SIunits.TemperatureDifference dT1_HTR    = (r_HTR_hin.T - r_HTR_cout.T);
  Modelica.SIunits.TemperatureDifference dT2_HTR    = (r_HTR_hout.T - r_HTR_cin.T);
  Real                                   T_ltmd_HTR = if dT1_HTR  > 0 and dT2_HTR > 0 then (dT2_HTR - dT1_HTR) / Modelica.Math.log(abs(dT2_HTR / dT1_HTR)) else -1;
  Real                                   UA_HTR     = if T_ltmd_HTR > 0 then Q_HTR / T_ltmd_HTR else 0.0;
  // LTR
  Modelica.SIunits.Power                 Q_LTR      = (r_LTR_cout.h - r_LTR_cin.h) * r_LTR_cin.w / 1e6 "W->MW";
  Modelica.SIunits.TemperatureDifference dT1_LTR    = (r_LTR_hin.T - r_LTR_cout.T);
  Modelica.SIunits.TemperatureDifference dT2_LTR    = (r_LTR_hout.T - r_LTR_cin.T);
  Real                                   T_ltmd_LTR = if dT1_LTR > 0 and dT2_LTR > 0 then (dT2_LTR - dT1_LTR) / Modelica.Math.log(abs(dT2_LTR / dT1_LTR)) else -1;
  Real                                   UA_LTR     = if T_ltmd_LTR > 0 then Q_LTR / T_ltmd_LTR else 0.0;
  
  // Liquid Na exit temperature
  Modelica.SIunits.Temperature T_heater_hot_out = r_HTR_hout.T;    

protected

  // performance map for main compressor
  parameter Real tableEta_comp_mc[5, 4]  = [0, 95, 100, 105; 1, 0.85310219, 0.837591241, 0.832420925; 2, 0.868613139, 0.857238443, 0.851034063; 3, 0.860340633, 0.85, 0.842761557; 4, 0.85310219, 0.839659367, 0.816909976];
  parameter Real tablePhic_comp_mc[5, 4] = [0, 95, 100, 105; 1, 0.000134346, 0.000150832, 0.000164161; 2, 0.000137854, 0.000153638, 0.00016802; 3, 0.000142414, 0.000158549, 0.000169774; 4, 0.000145921, 0.000161706, 0.000171528];
  parameter Real tablePR_comp_mc[5, 4]   = [0, 95, 100, 105; 1, 1.967529638, 2.350588505, 2.785882673; 2, 1.915294338, 2.315764972, 2.681412073; 3, 1.810823737, 2.220000255, 2.524706172; 4, 1.654117837, 2.115529655, 2.359294389];

  // performance map for re compressor      
  parameter Real tableEta_comp_rc[5, 4]  = [0, 95, 100, 105; 1, 0.85310219, 0.837591241, 0.832420925; 2, 0.868613139, 0.857238443, 0.851034063; 3, 0.860340633, 0.85, 0.842761557; 4, 0.85310219, 0.839659367, 0.816909976];
  parameter Real tablePhic_comp_rc[5, 4] = [0, 95, 100, 105; 1, 7.17663E-05, 8.05731E-05, 8.76935E-05; 2, 7.36401E-05, 8.20721E-05, 8.97547E-05; 3, 7.6076E-05, 8.46954E-05, 9.06916E-05; 4, 7.79498E-05, 8.63819E-05, 9.16285E-05];
  parameter Real tablePR_comp_rc[5, 4]   = [0, 95, 100, 105; 1, 1.967529638, 2.350588505, 2.785882673; 2, 1.915294338, 2.315764972, 2.681412073; 3, 1.810823737, 2.220000255, 2.524706172; 4, 1.654117837, 2.115529655, 2.359294389];
  
equation  
  connect(source.flange, r_turb_in.inlet) 
  annotation(
    Line(points = {{-102, -7}, {-101, -7}, {-101, -28}, {-100, -28}}, color = {0, 0, 255}, thickness = 1));    
  connect(r_turb_in.outlet, turbine.inlet) 
  annotation(
    Line(points = {{-96, -28}, {-63, -28}, {-63, -32}, {-62, -32}}, color = {0, 0, 255}, thickness = 1)); 
    /*
    // Case I - variable speed - let FPC determine the shaft speed
    connect(r_turb_in.T, FPC_turb.in_T1)
    annotation(
      Line(points = {{-96, -28}, {-100, -28}, {-100, -41}, {-95, -41}}, color = {159, 159, 223}, pattern = LinePattern.Dash));
    connect(r_turb_in.w, FPC_turb.in_w1)
    annotation(
      Line(points = {{-98, -28}, {-98, -38},{-95, -38}}, color = {159, 159, 223}, pattern = LinePattern.Dash));
    connect(p_const_TIP.y, FPC_turb.in_p1)
    annotation(
      Line(points = {{-140, 64}, {-135, 64}, {-135, 72}, {-135, -44}, {-94, -44}}, color = {0, 0, 127}, pattern = LinePattern.Dash)); // set the target inlet pressure 
    connect(p_const_CIP.y, FPC_turb.in_p2)
    annotation(
      Line(points = {{-140.5, 79}, {-132, 79}, {-132, -32}, {-90, -32}, {-90, -36}}, color = {0, 0, 127}, pattern = LinePattern.Dash)); // set the target outlet pressure    
    connect(FPC_turb.omega, speed_turb.w_ref)
    annotation(
      Line(points = {{-85, -41}, {-79, -41}}, color = {0, 0, 127}));     
    connect(speed_turb.flange, turbine.shaft_a)
    annotation(
      Line(points = {{-68, -41}, {-68, -40}, {-60, -40}}));
    */
    
    // Test case II or common case - constant speed
    connect(turbine.shaft_b, const_speed_turb.flange)     
    annotation(
      Line(points = {{-68, -41}, {-68, -40}, {-60, -40}}));
    
    
  connect(turbine.outlet, r_turb_out.inlet) 
  annotation(
    Line(points = {{-46, -32}, {-42, -32}}, color = {170, 170, 255}, thickness = 1));      
  connect(r_turb_out.outlet, r_HTR_hin.inlet) 
  annotation(
    Line(points = {{-42, -28}, {-42, 32}, {-38, 32}}, color = {170, 170, 255}, thickness = 1));
  connect(r_HTR_hin.outlet, HTR.gasIn) 
  annotation(
    Line(points = {{-34, 32}, {-30, 32}}, color = {170, 170, 255}, thickness = 0.5));
  connect(HTR.gasOut, r_HTR_hout.inlet) 
  annotation(
    Line(points = {{-14, 32}, {-10, 32}}, color = {170, 170, 255}, thickness = 1));
  connect(r_HTR_hout.outlet, r_LTR_hin.inlet) 
  annotation(
    Line(points = {{-6, 32}, {4, 32}}, color = {170, 170, 255}, thickness = 1));
  connect(r_LTR_hin.outlet, LTR.gasIn) annotation(
    Line(points = {{8, 32}, {12, 32}}, color = {170, 170, 255}, thickness = 1));
  connect(LTR.gasOut, r_LTR_hout.inlet) 
  annotation(
    Line(points = {{28, 32}, {32, 32}}, color = {170, 170, 255}, thickness = 1));
  connect(r_LTR_hout.outlet, splitter.inlet)
  annotation(
    Line(points = {{36, 32}, {37, 32}, {37, 24}}, color = {170, 170, 255}, thickness = 1));
  connect(splitter.outlet1, r_cooler_hin.inlet) 
  annotation(
    Line(points = {{35, 18}, {35, -67}}, color = {170, 170, 255}, thickness = 0.5));
        
    connect(splitter.outlet2, r_recomp_in.inlet) 
    annotation(
      Line(points = {{39, 18}, {45, 18}, {45, 19}}, color = {170, 170, 255}, thickness = 0.5)); 
     
    connect(r_recomp_in.outlet, recompressor.inlet) 
    annotation(
      Line(points = {{50, 19}, {50, 19.5}, {90, 19.5}, {90, 6}}, color = {170, 170, 255}, thickness = 0.5)); 
      /*
      // for Case I - variable speed, fixed inlet and outlet pressure
      connect(p_const_CIP.y, FPC_rc.in_p1)
      annotation(
        Line(points = {{-140.5, 79}, {44, 79}, {44, -6}, {60, -6}}, color = {0, 0, 127}, pattern = LinePattern.Dash));
      connect(r_recomp_in.T, FPC_rc.in_T1) 
      annotation(
        Line(points = {{50, 19}, {47.5, 19}, {47.5, -3}, {60, -3}}, color = {159, 159, 223}, pattern = LinePattern.Dash));
      connect(r_recomp_in.w, FPC_rc.in_w1) 
      annotation(
        Line(points = {{50, 19}, {50.5, 19}, {50.5, 0}, {60, 0}}, color = {159, 159, 223}, pattern = LinePattern.Dash));    
      connect(p_const_TIP.y, FPC_rc.in_p2)
      annotation(
        Line(points = {{-140, 64}, {-135, 64}, {-135, 72}, {64, 72}, {64, 2}, {66, 2}}, color = {0, 0, 127}, pattern = LinePattern.Dash));  // set the target outlet pressure         
      */
      
      // Case II - variable speed, fixed outlet pressure only
      connect(r_recomp_in.p, FPC_rc.in_p1)
      annotation(
        Line(points = {{44.5, 19}, {44.5, -6}, {60, -6}}, color = {159, 159, 223}, pattern = LinePattern.Dash));
      connect(r_recomp_in.T, FPC_rc.in_T1)
      annotation(
        Line(points = {{47.5, 19}, {47.5, -3}, {60, -3}}, color = {159, 159, 223}, pattern = LinePattern.Dash));
      connect(r_recomp_in.w, FPC_rc.in_w1)
      annotation(
        Line(points = {{50.5, 19}, {50.5, 0}, {60, 0}}, color = {159, 159, 223}, pattern = LinePattern.Dash));          
      connect(ramp_p.y, FPC_rc.in_p2)
      annotation(
        Line(points = {{-140.5, 79}, {44, 79}, {65, 79}, {65, 2}}, color = {0, 0, 127}, pattern = LinePattern.Dash));   // set the target outlet pressure following the source pressure            

      connect(FPC_rc.omega, speed_rc.w_ref)
      annotation(
        Line(points = {{70, -3}, {68, -3}, {68, -2}, {75, -2}, {75, -3}}, color = {0, 0, 127}, pattern = LinePattern.Dash));            
      connect(speed_rc.flange, recompressor.shaft_a) 
      annotation(
        Line(points = {{86, -3}, {82.75, -3}, {82.75, -2}, {92, -2}, {92, -1}}));   

    connect(recompressor.outlet, r_recomp_out.inlet) 
    annotation(
      Line(points = {{104, 6}, {104, 22}}, color = {0, 0, 255}, thickness = 0.5));
      
      /*
      // no speed control (neither Case I nor Case II), fixed re-compressor speed               
      connect(recompressor.shaft_b, const_speed_rc.flange) annotation(
        Line(points = {{86, 7}, {74, 7}}));  
      */
      
    connect(r_recomp_out.outlet, mixer.inlet2) 
    annotation(
      Line(points = {{104, 26}, {104, 56}, {-2, 56}}, color = {0, 0, 255}, thickness = 0.5));
      
  connect(r_cooler_hin.outlet, cooler.gasIn) 
  annotation(
    Line(points = {{42, -67}, {48, -67}, {48, -66}}, color = {170, 170, 255}, thickness = 0.5));
  connect(cooler.gasOut, r_cooler_hout.inlet) 
  annotation(
    Line(points = {{64, -66}, {68, -66}, {68, -67}, {72, -67}}, color = {170, 170, 255}, thickness = 0.5));
  connect(r_cooler_hout.outlet, r_comp_in.inlet) 
  annotation(
    Line(points = {{78, -67}, {110, -67}, {110, 18}, {114, 18}, {114, 19}}, color = {170, 170, 255}, thickness = 0.5));
  connect(r_comp_in.outlet, compressor.inlet) 
  annotation(
    Line(points = {{120, 19}, {120, 20}, {156, 20}, {156, 2}}, color = {170, 170, 255}, thickness = 0.5));
    
    /*
    // Case I - ramp mdot, controller determined varialbe speed, fixed in_p1 of controller
    connect(p_const_CIP.y, FPC_mc.in_p1)
    annotation(
      Line(points = {{-140.5, 79}, {114, 79}, {114, -6}, {124, -6}}, color = {0, 0, 127}, pattern = LinePattern.Dash));
    connect(r_comp_in.T, FPC_mc.in_T1)
    annotation(
      Line(points = {{120, 19}, {117, 19}, {117, -3}, {124, -3}}, color = {159, 159, 223}, pattern = LinePattern.Dash));
    connect(r_comp_in.w, FPC_mc.in_w1)
    annotation(
      Line(points = {{120, 19}, {119.5, 19}, {119.5, 0}, {124, 0}}, color = {159, 159, 223}, pattern = LinePattern.Dash));     
    connect(p_const_TIP.y, FPC_mc.in_p2)
    annotation(
      Line(points = {{-140, 64}, {-135, 64}, {-135, 72}, {130, 72}, {130, 2}}, color = {0, 0, 127}, pattern = LinePattern.Dash));  // set the target outlet pressure       
    */
    
    // Case II - ramp source pressure, controller determined variable speed
    connect(r_comp_in.p, FPC_mc.in_p1)    
    annotation(
      Line(points = {{114.5, 19}, {114.5, -6}, {124, -6}}, color = {159, 159, 223}, pattern = LinePattern.Dash)); 
    connect(r_comp_in.T, FPC_mc.in_T1)    
    annotation(
      Line(points = {{117, 19}, {117, -3}, {124, -3}}, color = {159, 159, 223}, pattern = LinePattern.Dash));
    connect(r_comp_in.w, FPC_mc.in_w1)    
    annotation(
      Line(points = {{119.5, 19}, {119.5, 0}, {124, 0}}, color = {159, 159, 223}, pattern = LinePattern.Dash));     
     
    connect(ramp_p.y, FPC_mc.in_p2)      
    annotation(
        Line(points = {{-140.5, 79}, {44, 79}, {129, 79}, {129, 2}}, color = {0, 0, 127}, pattern = LinePattern.Dash));  // set the target outlet pressure following the source pressure      
    

    connect(FPC_mc.omega, speed_mc.w_ref)
    annotation(
      Line(points = {{134, -3}, {134, -2}, {140, -2}, {140, -3}, {139, -3}}, color = {0, 0, 127}, pattern = LinePattern.Dash));         
    connect(speed_mc.flange, compressor.shaft_a)
    annotation(
      Line(points = {{150, -3}, {157, -3}, {157, -4}}));      
  
  connect(compressor.outlet, r_comp_out.inlet) 
  annotation(
    Line(points = {{168, 2}, {168, 22}}, color = {0, 0, 255}, thickness = 0.5));
    
    /*
    // neither Case I nor Case II - constant speed for common situation
    connect(compressor.shaft_b, const_speed_mc.flange) annotation(
      Line(points = {{128, 7}, {122.5, 7}, {122.5, 6}, {117, 6}}));
    */
    
  connect(r_comp_out.outlet, r_LTR_cin.inlet) 
  annotation(
    Line(points = {{168, 26}, {168, 64}, {20, 64}, {20, 48}}, color = {0, 0, 255}, thickness = 0.5));
  connect(r_LTR_cin.outlet, LTR.waterIn) 
  annotation(
    Line(points = {{20, 44}, {20, 40}}, color = {0, 0, 255}, thickness = 0.5));
  connect(LTR.waterOut, r_LTR_cout.inlet) 
  annotation(
    Line(points = {{20, 24}, {20, 18}}, color = {0, 0, 255}, thickness = 0.5));
  connect(r_LTR_cout.outlet, mixer.inlet1) 
  annotation(
    Line(points = {{20, 14}, {-2, 14}, {-2, 52}}, color = {0, 0, 255}, thickness = 0.5));
  connect(mixer.outlet, r_HTR_cin.inlet) annotation(
    Line(points = {{-6, 54}, {-22, 54}, {-22, 48}}, color = {0, 0, 255}, thickness = 1));
  connect(r_HTR_cin.outlet, HTR.waterIn) 
  annotation(
    Line(points = {{-22, 44}, {-22, 40}}, color = {0, 0, 255}, thickness = 1));
  connect(HTR.waterOut, r_HTR_cout.inlet) 
  annotation(
    Line(points = {{-22, 24}, {-22, 18}}, color = {0, 0, 255}, thickness = 1));
  connect(r_HTR_cout.outlet, r_heater_cin.inlet) 
  annotation(
    Line(points = {{-22, 14}, {-48, 14}, {-48, 32}, {-54, 32}}, color = {0, 0, 255}, thickness = 1));
  connect(r_heater_cin.outlet, heater.gasIn) 
  annotation(
    Line(points = {{-58, 32}, {-66, 32}}, color = {0, 0, 255}, thickness = 1));
  connect(heater.gasOut, r_heater_cout.inlet) 
  annotation(
    Line(points = {{-82, 32}, {-88, 32}}, color = {0, 0, 255}, thickness = 1));
  connect(r_heater_cout.outlet, sink.flange) annotation(
    Line(points = {{-92, 32}, {-95, 32}, {-95, 33}, {-100, 33}}, color = {0, 0, 255}, thickness = 1));
    
  // hot stream for heater
  connect(source_heater_hot.flange, r_heater_hin.inlet) annotation(
    Line(points = {{-74, 8}, {-74, 14}}, color = {255, 0, 0}, thickness = 0.75));
  connect(r_heater_hin.outlet, heater.waterIn) annotation(
    Line(points = {{-74, 18}, {-74, 24}}, color = {255, 0, 0}, thickness = 0.75));
  connect(heater.waterOut, r_heater_hout.inlet) annotation(
    Line(points = {{-74, 40}, {-74, 48}}, color = {255, 0, 0}, thickness = 0.75));
  connect(r_heater_hout.outlet, sink_heater_hot.flange) annotation(
    Line(points = {{-74, 52}, {-74, 58}}, color = {255, 0, 0}, thickness = 0.75));
    
  // cold side of cooler
  connect(source_cooler.flange, r_cooler_cin.inlet) annotation(
    Line(points = {{56, -42}, {56, -48}}, color = {0, 255, 0}, thickness = 0.75));    
  connect(r_cooler_cin.outlet, cooler.waterIn) annotation(
    Line(points = {{56, -52}, {56, -58}}, color = {0, 255, 0}, thickness = 0.75));
  connect(cooler.waterOut, r_cooler_cout.inlet) annotation(
    Line(points = {{56, -74}, {56, -82}}, color = {0, 255, 0}, thickness = 0.75));
  connect(r_cooler_cout.outlet, sink_cooler.flange) annotation(
    Line(points = {{56, -86}, {56, -93}, {57, -93}}, color = {0, 255, 0}, thickness = 0.75));
  /*
  // ramp input  
  // for case I  
  connect(ramp_mdot.y, source.in_w0)
  annotation(
    Line(points = {{-140.5, 91}, {-128, 91}, {-128, -0.75}, {-110, -0.75}, {-110, -4.5}}, color = {0, 0, 127}, pattern = LinePattern.Dash));
    
  // set target outlet pressure
  connect(p_const_TIP.y, sink.in_p0)
  annotation(
    Line(points = {{-140, 64}, {-135, 64}, {-135, 72}, {-102, 72}, {-102, 36}}, color = {0, 0, 127}, pattern = LinePattern.Dash));
  */ 
  
  // for case II 
  connect(ramp_p.y, source.in_p0)   
  annotation(
    Line(points = {{-140.5, 79}, {-128, 79}, {-128, -0.75}, {-110, -0.75}, {-110, -4.5}}, color = {0, 0, 127}, pattern = LinePattern.Dash));  
  
  
  /*
  connect(ramp_mdot.y, sink.in_w0) annotation(
    Line(points = {{-113.5, 21}, {-102, 21}, {-102, 30.5}}, color = {0, 0, 127}));
  */
  
  
  annotation(
    Diagram(coordinateSystem(extent = {{-160, -100}, {180, 100}})),
    experiment(StartTime = 0, StopTime = 10, Tolerance = 1e-3, Interval = 1),
    __OpenModelica_commandLineOptions = "--matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection -d=initialization,NLSanalyticJacobian,bltdump",
    // __OpenModelica_simulationFlags(lv = "LOG_DEBUG,LOG_NLS,LOG_NLS_V,LOG_STATS,LOG_INIT,LOG_STDOUT, -w", outputFormat = "mat", s = "dassl", nls = "homotopy"));
    __OpenModelica_simulationFlags(lv = "LOG_DEBUG,LOG_NLS,LOG_STATS,LOG_INIT,LOG_STDOUT -w", newtonFTol = "1e-6", newtonXTol= "1e-6", outputFormat = "mat", s = "dassl", nls = "homotopy"));
end TPDyn_RCBCycle;
