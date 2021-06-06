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
    table_k_LTR_wall = table_k_metalwall,
    table_k_HTR_wall = table_k_metalwall,
    
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
  
  // set the values of parameters accordingly
  parameter Model.HeatExchangerConfig cfg_heater = cfg.cfg_heater;
  parameter Model.HeatExchangerConfig cfg_LTR    = cfg.cfg_LTR;
  parameter Model.HeatExchangerConfig cfg_HTR    = cfg.cfg_HTR;
  parameter Model.TurbomachineryConfig cfg_turb  = cfg.cfg_turb;
  parameter Model.TurbomachineryConfig cfg_comp   = cfg.cfg_comp;
  parameter Model.TurbomachineryConfig cfg_recomp = cfg.cfg_recomp;
  parameter Model.HeatExchangerConfig cfg_cooler  = cfg.cfg_cooler;
  parameter Model.SplitterConfig cfg_splitter     = cfg.cfg_splitter;
  parameter Model.SplitterConfig cfg_merger      = cfg.cfg_merger;
  
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
    Placement(visible = true, transformation(origin = {-74, 2}, extent = {{-4, -4}, {4, 4}}, rotation = 90)));

  ThermoPower.Water.SinkPressure sink_heater_hot(
    redeclare package Medium = medium_heater, 
    p0    = cfg_heater.cfg_hot.st_out.p,
    T     = cfg_heater.cfg_hot.st_out.T,
    use_T = false) 
  annotation(
    Placement(visible = true, transformation(origin = {-74, 66}, extent = {{-4, -4}, {4, 4}}, rotation = 90)));

  Steps.TPComponents.HE heater(
    redeclare package FluidMedium = medium_heater, 
    redeclare package FlueGasMedium = medium_main, 
    
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
      Tstartin  = cfg_heater.cfg_gas.st_in.T,
      Tstartout = cfg_heater.cfg_gas.st_out.T,
      Nt        = cfg_heater.cfg_gas.N_ch,
      Dhyd      = cfg_heater.cfg_gas.geo_area.d_hyd
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
      hstartin                = cfg_heater.cfg_fluid.st_in.h,
      hstartout               = cfg_heater.cfg_fluid.st_out.h
    ),  
    
    // // Dittus Boelter
    // redeclare replaceable model HeatTransfer_F = ThermoPower.Thermal.HeatTransferFV.DittusBoelter,    
    // redeclare replaceable model HeatTransfer_G = ThermoPower.Thermal.HeatTransferFV.DittusBoelter,
    // fluidFlow(
    //   heatTransfer(heating=true), 
    //   fixedMassFlowSimplified = true,
    //   hstartin                = cfg_heater.cfg_fluid.st_in.h,
    //   hstartout               = cfg_heater.cfg_fluid.st_out.h),   // set the fluid flow as fixed mdot for simplarity
    // gasFlow(
    //   heatTransfer(heating=false), 
    //   Tstartin  = cfg_heater.cfg_gas.st_in.T,
    //   Tstartout = cfg_heater.cfg_gas.st_out.T,
    //   Nt        = cfg_heater.cfg_gas.N_ch,
    //   Dhyd      = cfg_heater.cfg_gas.geo_area.d_hyd),
    
    
    // redeclare replaceable model HeatTransfer_F = ThermoPower.Thermal.HeatTransferFV.IdealHeatTransfer,    
    // redeclare replaceable model HeatTransfer_G = ThermoPower.Thermal.HeatTransferFVH.IdealHeatTransfer, 
    // fluidFlow(
    //   heatTransfer(gamma = cfg_heater.cfg_hot.gamma_HE),
    //   fixedMassFlowSimplified = true,
    //   hstartin                = cfg_heater.cfg_hot.st_in.h,
    //   hstartout               = cfg_heater.cfg_hot.st_out.h),   // set the fluid flow as fixed mdot for simplarity
    // gasFlow(
    //   heatTransfer(gamma = cfg_heater.cfg_cold.gamma_HE),
    //   Tstartin    = cfg_heater.cfg_cold.st_in.T,
    //   Tstartout   = cfg_heater.cfg_cold.st_out.T),
    
    // redeclare replaceable model HeatTransfer_F = ThermoPower.Thermal.HeatTransferFV.ConstantHeatTransferCoefficient(),     
    // // redeclare replaceable model HeatTransfer_G = ThermoPower.Thermal.HeatTransferFV.IdealHeatTransfer,      
    // redeclare replaceable model HeatTransfer_G = ThermoPower.Thermal.HeatTransferFV.ConstantHeatTransferCoefficient(),       

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
  
  Steps.TPComponents.WaterStateReader r_heater_hin(redeclare package Medium = medium_heater) annotation(
    Placement(visible = true, transformation(origin = {-74, 14}, extent = {{-4, -4}, {4, 4}}, rotation = 90)));
  Steps.TPComponents.WaterStateReader r_heater_hout(redeclare package Medium = medium_heater) annotation(
    Placement(visible = true, transformation(origin = {-74, 52}, extent = {{-4, -4}, {4, 4}}, rotation = 90)));
  Steps.TPComponents.GasStateReader r_heater_cin(redeclare package Medium = medium_main) annotation(
    Placement(visible = true, transformation(origin = {-56, 32}, extent = {{-4, -4}, {4, 4}}, rotation = 180)));
  Steps.TPComponents.GasStateReader r_heater_cout(redeclare package Medium = medium_main) annotation(
    Placement(visible = true, transformation(origin = {-90, 32}, extent = {{-4, -4}, {4, 4}}, rotation = 180)));

/*
  ThermoPower.Gas.SourceMassFlow source(
    redeclare package Medium = medium_main, 
    T         = st_source.T,
    p0        = st_source.p,
    w0        = st_source.mdot,
    use_in_T  = false,
    use_in_w0 = false,
    gas(
      p(start = st_source.p, nominal = st_source.p), 
      T(start = st_source.T, nominal = st_source.T)))
      //h(start = st_source.h, nominal = st_source.h)))
  annotation(
    Placement(visible = true, transformation(origin = {-105, -39}, extent = {{-5, -5}, {5, 5}}, rotation = 0)));
*/

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
    Placement(visible = true, transformation(origin = {-105, -39}, extent = {{-5, -5}, {5, 5}}, rotation = 0)));


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
  ThermoPower.Gas.SinkPressure sink(
    redeclare package Medium = medium_main, 
    p0 = st_sink.p,
    T  = st_sink.T,
    use_in_p0 = true)
  annotation(
    Placement(transformation(origin = {0, -80}, extent = {{-10, -10}, {10, 10}}, rotation = 270)));
*/


  ThermoPower.Gas.SinkMassFlow sink(
    redeclare package Medium = medium_main, 
    T        = st_sink.T,
    p0       = st_sink.p,
    use_in_T = false,
    use_in_w0 = true,
    w0       = st_sink.mdot
  ) 
  annotation(
    Placement(visible = true, transformation(origin = {-105, 33}, extent = {{-5, -5}, {5, 5}}, rotation = 180)));

/*
  ThermoPower.Gas.SourceMassFlow source_temp(
    redeclare package Medium = medium_main, 
    T        = st_source_temp.T,
    p0       = st_source_temp.p,
    w0       = st_source_temp.mdot,
    use_in_T = false,  
    use_in_w0 = false,  
    gas(
      p(start = st_source_temp.p, nominal = st_source_temp.p), 
      T(start = st_source_temp.T, nominal = st_source_temp.T)))
      // h(start = st_source_temp.h, nominal = st_source_temp.h)))
  annotation(
    Placement(transformation(origin = {0, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 270)));
*/
/*
  ThermoPower.Gas.SourcePressure source_temp(
    redeclare package Medium = medium_main, 
    T        = st_source_temp.T,
    p0       = st_source_temp.p,
    use_in_T = false,
    gas(
      p(start = st_source_temp.p, nominal = st_source_temp.p), 
      T(start = st_source_temp.T, nominal = st_source_temp.T),
      h(start = st_source_temp.h, nominal = st_source_temp.h))) 
  annotation(
    Placement(transformation(origin = {0, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 270)));
*/
/*
  ThermoPower.Gas.SinkPressure sink_temp(
    redeclare package Medium = medium_main,
    T  = st_sink_temp.T,
    p0 = st_sink_temp.p)
  annotation(
    Placement(transformation(extent = {{60, -10}, {80, 10}}, rotation = 0)));
*/
/*
  ThermoPower.Gas.SinkMassFlow sink_temp(
    redeclare package Medium = medium_main, 
    T        = st_sink_temp.T,
    p0       = st_sink_temp.p,
    use_in_T = false,
    w0       = st_sink_temp.mdot
  );
*/

  // use FlowJoin to mix flow
  Gas.FlowJoin mixer(redeclare package Medium = medium_main) 
  annotation(
    Placement(visible = true, transformation(origin = {-4, 54}, extent = {{-4, -4}, {4, 4}}, rotation = 180))); 

  Steps.TPComponents.PCHE LTR(
    redeclare package FluidMedium = medium_main, 
    redeclare package FlueGasMedium = medium_main, 
     
    // use Marchionni PCHE HeatTransfer
    // slow but can have a result - set a_phi = 0 to use Gnielinski's correlation 
    redeclare replaceable model HeatTransfer_F = Steps.TPComponents.NgoHeatTransferFV(),
    redeclare replaceable model HeatTransfer_G = Steps.TPComponents.NgoHeatTransferFV(),
    gasFlow(      
      heatTransfer(
        pitch = cfg_LTR.cfg_hot.l_pitch,
        phi   = cfg_LTR.cfg_hot.a_phi,
        Cf_C1 = Cf_C1,
        Cf_C2 = Cf_C2,
        Cf_C3 = Cf_C3,
        gamma_min = 2000,
        gamma_max = 5000        
      )
    ),
    fluidFlow(      
      heatTransfer(
        pitch = cfg_LTR.cfg_cold.l_pitch,
        phi   = cfg_LTR.cfg_cold.a_phi,
        Cf_C1 = Cf_C1,
        Cf_C2 = Cf_C2,
        Cf_C3 = Cf_C3,
        gamma_min = 2000,
        gamma_max = 5000        
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
  )
  annotation(
    Placement(visible = true, transformation(origin = {20, 32}, extent = {{-8, -8}, {8, 8}}, rotation = 0)));

  Steps.TPComponents.GasStateReader r_LTR_hin(redeclare package Medium = medium_main) annotation(
    Placement(visible = true, transformation(origin = {6, 32}, extent = {{-4, -4}, {4, 4}}, rotation = 0)));
  Steps.TPComponents.GasStateReader r_LTR_hout(redeclare package Medium = medium_main) annotation(
    Placement(visible = true, transformation(origin = {34, 32}, extent = {{-4, -4}, {4, 4}}, rotation = 0)));
  Steps.TPComponents.GasStateReader r_LTR_cin(redeclare package Medium = medium_main) annotation(
    Placement(visible = true, transformation(origin = {20, 46}, extent = {{-4, -4}, {4, 4}}, rotation = -90)));
  Steps.TPComponents.GasStateReader r_LTR_cout(redeclare package Medium = medium_main) annotation(
    Placement(visible = true, transformation(origin = {20, 16}, extent = {{-4, -4}, {4, 4}}, rotation = -90)));

  Steps.TPComponents.PCHE HTR(
    redeclare package FluidMedium = medium_main, 
    redeclare package FlueGasMedium = medium_main, 
     
    // use Marchionni PCHE HeatTransfer
    // slow but can have a result - set a_phi = 0 to use Gnielinski's correlation 
    
    redeclare replaceable model HeatTransfer_F = Steps.TPComponents.GnielinskiHeatTransferFV(),
    redeclare replaceable model HeatTransfer_G = Steps.TPComponents.GnielinskiHeatTransferFV(),
    gasFlow(      
      heatTransfer(
        pitch = cfg_HTR.cfg_hot.l_pitch,
        phi   = cfg_HTR.cfg_hot.a_phi,
        Cf_C1 = Cf_C1,
        Cf_C2 = Cf_C2,
        Cf_C3 = Cf_C3,
        gamma_min = 2000,
        gamma_max = 5000        
      )
    ),
    fluidFlow(      
      heatTransfer(
        pitch = cfg_HTR.cfg_cold.l_pitch,
        phi   = cfg_HTR.cfg_cold.a_phi,
        Cf_C1 = Cf_C1,
        Cf_C2 = Cf_C2,
        Cf_C3 = Cf_C3,
        gamma_min = 2000,
        gamma_max = 5000
      )
    ),        

    // redeclare replaceable model HeatTransfer_F = ThermoPower.Thermal.HeatTransferFV.DittusBoelter,    
    // redeclare replaceable model HeatTransfer_G = ThermoPower.Thermal.HeatTransferFV.DittusBoelter,
    // gasFlow(heatTransfer(heating=true)),
    // fluidFlow(heatTransfer(heating=false)),    
    
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
  )
  annotation(
    Placement(visible = true, transformation(origin = {-22, 32}, extent = {{-8, -8}, {8, 8}}, rotation = 0)));
  
  Steps.TPComponents.GasStateReader r_HTR_hin(redeclare package Medium = medium_main) annotation(
    Placement(visible = true, transformation(origin = {-36, 32}, extent = {{-4, -4}, {4, 4}}, rotation = 0)));
  Steps.TPComponents.GasStateReader r_HTR_hout(redeclare package Medium = medium_main) annotation(
    Placement(visible = true, transformation(origin = {-8, 32}, extent = {{-4, -4}, {4, 4}}, rotation = 0)));
  Steps.TPComponents.GasStateReader r_HTR_cin(redeclare package Medium = medium_main) annotation(
    Placement(visible = true, transformation(origin = {-22, 46}, extent = {{-4, -4}, {4, 4}}, rotation = -90)));
  Steps.TPComponents.GasStateReader r_HTR_cout(redeclare package Medium = medium_main) annotation(
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
    Placement(visible = true, transformation(origin = {-72, -46}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));

  Modelica.Mechanics.Rotational.Sources.ConstantSpeed const_speed_turb(
      w_fixed=cfg_turb.N, useSupport=false) annotation(
    Placement(visible = true, transformation(origin = {-49, -47}, extent = {{5, -5}, {-5, 5}}, rotation = 0)));
 
  Steps.TPComponents.GasStateReader r_turb_in(redeclare package Medium = medium_main) annotation(
    Placement(visible = true, transformation(origin = {-92, -38}, extent = {{-4, -4}, {4, 4}}, rotation = 0)));

  Steps.TPComponents.GasStateReader r_turb_out(redeclare package Medium = medium_main) annotation(
    Placement(visible = true, transformation(origin = {-42, -30}, extent = {{-4, -4}, {4, 4}}, rotation = 90)));
 
  ThermoPower.Water.SourceMassFlow source_cooler(
    redeclare package Medium = medium_cooler,
    w0 = cfg_cooler.cfg_cold.st_in.mdot,
    p0 = cfg_cooler.cfg_cold.st_in.p,
    h  = cfg_cooler.cfg_cold.st_in.h,
    T  = cfg_cooler.cfg_cold.st_in.T) 
  annotation(
    Placement(visible = true, transformation(origin = {56, -18}, extent = {{-6, -6}, {6, 6}}, rotation = 270)));
    
  ThermoPower.Water.SinkPressure sink_cooler(
    redeclare package Medium = medium_cooler, 
    p0    = cfg_cooler.cfg_cold.st_out.p,
    T     = cfg_cooler.cfg_cold.st_out.T,
    use_T = false) 
  annotation(
    Placement(visible = true, transformation(origin = {57, -77}, extent = {{-5, -5}, {5, 5}}, rotation = 270)));

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
    SSInit          = true,
    gasQuasiStatic  = false,
    FluidPhaseStart = ThermoPower.Choices.FluidPhase.FluidPhases.Liquid,
    metalTube(WallRes=false))
  annotation(
    Placement(visible = true, transformation(origin = {56, -46}, extent = {{-8, -8}, {8, 8}}, rotation = 0)));
  
  Steps.TPComponents.GasStateReader r_cooler_hin(redeclare package Medium = medium_main) annotation(
    Placement(visible = true, transformation(origin = {39, -47}, extent = {{-5, -5}, {5, 5}}, rotation = 0)));
  Steps.TPComponents.GasStateReader r_cooler_hout(redeclare package Medium = medium_main) annotation(
    Placement(visible = true, transformation(origin = {75, -47}, extent = {{-5, -5}, {5, 5}}, rotation = 0)));
  Steps.TPComponents.WaterStateReader r_cooler_cin(redeclare package Medium = medium_cooler) annotation(
    Placement(visible = true, transformation(origin = {56, -30}, extent = {{-4, -4}, {4, 4}}, rotation = -90)));
  Steps.TPComponents.WaterStateReader r_cooler_cout(redeclare package Medium = medium_cooler) annotation(
    Placement(visible = true, transformation(origin = {56, -64}, extent = {{-4, -4}, {4, 4}}, rotation = -90)));

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
  annotation(
    Placement(visible = true, transformation(origin = {112, 6}, extent = {{-8, -8}, {8, 8}}, rotation = 0)));
            
  Modelica.Mechanics.Rotational.Sources.ConstantSpeed const_speed_mc(
    w_fixed    = cfg_comp.N,
    useSupport = false)
  annotation(
    Placement(visible = true, transformation(origin = {133, 7}, extent = {{5, -5}, {-5, 5}}, rotation = 0)));

  Steps.TPComponents.GasStateReader r_comp_in(redeclare package Medium = medium_main) annotation(
    Placement(visible = true, transformation(origin = {101, 3}, extent = {{-5, -5}, {5, 5}}, rotation = 90)));

  Steps.TPComponents.GasStateReader r_comp_out(redeclare package Medium = medium_main) annotation(
    Placement(visible = true, transformation(origin = {118, 22}, extent = {{-4, -4}, {4, 4}}, rotation = 90)));
 
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
  annotation(
    Placement(visible = true, transformation(origin = {69, 7}, extent = {{-9, -9}, {9, 9}}, rotation = 0)));
            
  Modelica.Mechanics.Rotational.Sources.ConstantSpeed const_speed_rc(
    w_fixed    = cfg_recomp.N,
    useSupport = false)
  annotation(
    Placement(visible = true, transformation(origin = {91, 7}, extent = {{5, -5}, {-5, 5}}, rotation = 0)));

  Steps.TPComponents.GasStateReader r_recomp_in(redeclare package Medium = medium_main) annotation(
    Placement(visible = true, transformation(origin = {53, 3}, extent = {{5, 5}, {-5, -5}}, rotation = -90)));

  Steps.TPComponents.GasStateReader r_recomp_out(redeclare package Medium = medium_main) annotation(
    Placement(visible = true, transformation(origin = {76, 22}, extent = {{-4, -4}, {4, 4}}, rotation = 90))); 

  ThermoPower.Gas.FlowSplit splitter(redeclare package Medium = medium_main) annotation(
    Placement(visible = true, transformation(origin = {37, -3}, extent = {{-5, 5}, {5, -5}}, rotation = -90)));

  // calculated variables
  
  // for turbine
  Modelica.SIunits.Power W_turb = (r_turb_in.h - r_turb_out.h) * r_turb_in.w / 1e6 "W->MW, net power for turbine";
  Modelica.SIunits.Efficiency eta_turb = turbine.eta * 100;
  // main compressor
  Modelica.SIunits.Power W_MC = (r_comp_out.h - r_comp_in.h) * r_comp_in.w / 1e6 "W->MW, net power for main compressor";
  Modelica.SIunits.Efficiency eta_MC = compressor.eta * 100;
  // re compressor
  Modelica.SIunits.Power W_RC = (r_recomp_out.h - r_recomp_in.h) * r_recomp_in.w / 1e6 "W->MW, net power for recompressor";
  Modelica.SIunits.Efficiency eta_RC = recompressor.eta * 100;
  
  // power block
  Modelica.SIunits.Power W_net = W_turb - W_MC - W_RC "net power generated";
  Modelica.SIunits.Power Q_heater = (r_heater_hin.h - r_heater_hout.h) * r_heater_hin.w / 1e6 "heat input for the power block";
  Modelica.SIunits.Power Q_cooler = (r_cooler_cout.h - r_cooler_cin.h) * r_cooler_cin.w / 1e6 "heat input for the power block";
  Modelica.SIunits.Efficiency eta_pb = W_net / Q_heater * 100 "power block efficiency";
  Real SR = r_recomp_in.w / r_comp_in.w * 100 "split ratio";

  // heat transfer coefficient for HTR and LTR
  // HTR
  Modelica.SIunits.Power Q_HTR = (r_HTR_cout.h - r_HTR_cin.h) * r_HTR_cin.w / 1e6 "W->MW, heat input for HTR";
  Modelica.SIunits.TemperatureDifference dT1_HTR = (r_HTR_hin.T - r_HTR_cout.T);
  Modelica.SIunits.TemperatureDifference dT2_HTR = (r_HTR_hout.T - r_HTR_cin.T);
  Real T_ltmd_HTR = if dT1_HTR  > 0 and dT2_HTR > 0 then (dT2_HTR - dT1_HTR) / Modelica.Math.log(abs(dT2_HTR / dT1_HTR)) else -1;
  Real UA_HTR = if T_ltmd_HTR > 0 then Q_HTR / T_ltmd_HTR else 0.0;
  // LTR
  Modelica.SIunits.Power Q_LTR = (r_LTR_cout.h - r_LTR_cin.h) * r_LTR_cin.w / 1e6 "W->MW";
  Modelica.SIunits.TemperatureDifference dT1_LTR = (r_LTR_hin.T - r_LTR_cout.T);
  Modelica.SIunits.TemperatureDifference dT2_LTR = (r_LTR_hout.T - r_LTR_cin.T);
  Real T_ltmd_LTR = if dT1_LTR > 0 and dT2_LTR > 0 then (dT2_LTR - dT1_LTR) / Modelica.Math.log(abs(dT2_LTR / dT1_LTR)) else -1;
  Real UA_LTR = if T_ltmd_LTR > 0 then Q_LTR / T_ltmd_LTR else 0.0;  
  
  // Liquid Na exit temperature
  Modelica.SIunits.Temperature T_heater_hot_out = r_HTR_hout.T;  
  
  // value to accelerate the simulation
  constant Integer SEC_PER_HOUR = integer(60 * 60 / time_scaling_factor); 
  // constant Real time_scaling_factor = 12; // 5 min = 1 hour
  constant Real time_scaling_factor = 1; // 1 hour = 1 hour
  // constant Integer SEC_PER_HOUR = 60 * 60;
  /*
  // ramp input to simulate the ramp change in ASPEN+ simulation
  // ramp change for case I and case II in Guan's report
  Steps.TPComponents.RampSeq_W ramp_mdot(    
    time_start = 3 * SEC_PER_HOUR,
    interval = 1 * SEC_PER_HOUR,
    duration = 0.15 * SEC_PER_HOUR,       
    offset = st_source.mdot
  );
  */
  
  // ramp input for case III in Guan's report
  
  Modelica.Blocks.Sources.Ramp ramp_mdot(
    final duration  = 0.15 * SEC_PER_HOUR,
    final startTime = 1 * SEC_PER_HOUR,
    
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
  connect(source.flange, r_turb_in.inlet) annotation(
    Line(points = {{-100, -39}, {-97, -39}, {-97, -38}, {-94, -38}}, color = {0, 0, 255}, thickness = 1));
  connect(r_turb_in.outlet, turbine.inlet) annotation(
    Line(points = {{-90, -38}, {-80, -38}}, color = {0, 0, 255}, thickness = 1));  
  connect(turbine.outlet, r_turb_out.inlet) annotation(
    Line(points = {{-64, -38}, {-42, -38}, {-42, -32}}, color = {170, 170, 255}, thickness = 1));
    
    connect(turbine.shaft_b, const_speed_turb.flange) annotation(
      Line(points = {{-54, -47}, {-60, -47}, {-60, -46}, {-66, -46}}));      

  connect(r_turb_out.outlet, r_HTR_hin.inlet) annotation(
    Line(points = {{-42, -28}, {-42, 32}, {-38, 32}}, color = {170, 170, 255}, thickness = 1));   
   connect(r_HTR_hin.outlet, HTR.gasIn) annotation(
    Line(points = {{-34, 32}, {-30, 32}}, color = {170, 170, 255}, thickness = 0.5));  

  connect(HTR.gasOut, r_HTR_hout.inlet) annotation(
    Line(points = {{-14, 32}, {-10, 32}}, color = {170, 170, 255}, thickness = 1));
  connect(r_HTR_hout.outlet, r_LTR_hin.inlet) annotation(
    Line(points = {{-6, 32}, {4, 32}}, color = {170, 170, 255}, thickness = 1));
  connect(r_LTR_hin.outlet, LTR.gasIn) annotation(
    Line(points = {{8, 32}, {12, 32}}, color = {170, 170, 255}, thickness = 1));
  connect(LTR.gasOut, r_LTR_hout.inlet) annotation(
    Line(points = {{28, 32}, {32, 32}}, color = {170, 170, 255}, thickness = 1));
  connect(r_LTR_hout.outlet, splitter.inlet) annotation(
    Line(points = {{36, 32}, {37, 32}, {37, 0}}, color = {170, 170, 255}, thickness = 1));
  connect(splitter.outlet1, r_cooler_hin.inlet) annotation(
    Line(points = {{35, -6}, {35, -5.5}, {36, -5.5}, {36, -47}}, color = {170, 170, 255}, thickness = 0.5));

    connect(splitter.outlet2, r_recomp_in.inlet) annotation(
      Line(points = {{39, -6}, {52, -6}, {52, 0}, {53, 0}}, color = {170, 170, 255}, thickness = 0.5));
    connect(r_recomp_in.outlet, recompressor.inlet) annotation(
      Line(points = {{53, 6}, {52, 6}, {52, 14}, {62, 14}}, color = {170, 170, 255}, thickness = 0.5));
    connect(recompressor.outlet, r_recomp_out.inlet) annotation(
      Line(points = {{76, 14}, {76, 20}}, color = {0, 0, 255}, thickness = 0.5));

      connect(recompressor.shaft_b, const_speed_rc.flange) annotation(
        Line(points = {{86, 7}, {74, 7}}));

    connect(r_recomp_out.outlet, mixer.inlet2) annotation(
      Line(points = {{76, 24}, {76, 56}, {-2, 56}}, color = {0, 0, 255}, thickness = 0.5));

  connect(r_cooler_hin.outlet, cooler.gasIn) annotation(
    Line(points = {{42, -47}, {48, -47}, {48, -46}}, color = {170, 170, 255}, thickness = 0.5));
  connect(cooler.gasOut, r_cooler_hout.inlet) annotation(
    Line(points = {{64, -46}, {68, -46}, {68, -47}, {72, -47}}, color = {170, 170, 255}, thickness = 0.5));
  connect(r_cooler_hout.outlet, r_comp_in.inlet) annotation(
    Line(points = {{78, -47}, {101, -47}, {101, 0}}, color = {170, 170, 255}, thickness = 0.5));
  connect(r_comp_in.outlet, compressor.inlet) annotation(
    Line(points = {{101, 6}, {100, 6}, {100, 12}, {106, 12}}, color = {170, 170, 255}, thickness = 0.5));

    connect(compressor.shaft_b, const_speed_mc.flange) annotation(
      Line(points = {{128, 7}, {122.5, 7}, {122.5, 6}, {117, 6}}));

  connect(compressor.outlet, r_comp_out.inlet) annotation(
    Line(points = {{118, 12}, {118, 20}}, color = {0, 0, 255}, thickness = 0.5));
  connect(r_comp_out.outlet, r_LTR_cin.inlet) annotation(
    Line(points = {{118, 24}, {118, 64}, {20, 64}, {20, 48}}, color = {0, 0, 255}, thickness = 0.5));
  connect(r_LTR_cin.outlet, LTR.waterIn) annotation(
    Line(points = {{20, 44}, {20, 40}}, color = {0, 0, 255}, thickness = 0.5));
  connect(LTR.waterOut, r_LTR_cout.inlet) annotation(
    Line(points = {{20, 24}, {20, 18}}, color = {0, 0, 255}, thickness = 0.5));
  connect(r_LTR_cout.outlet, mixer.inlet1) annotation(
    Line(points = {{20, 14}, {-2, 14}, {-2, 52}}, color = {0, 0, 255}, thickness = 0.5));
  connect(mixer.outlet, r_HTR_cin.inlet) annotation(
    Line(points = {{-6, 54}, {-22, 54}, {-22, 48}}, color = {0, 0, 255}, thickness = 1));
  connect(r_HTR_cin.outlet, HTR.waterIn) annotation(
    Line(points = {{-22, 44}, {-22, 40}}, color = {0, 0, 255}, thickness = 1));
  connect(HTR.waterOut, r_HTR_cout.inlet) annotation(
    Line(points = {{-22, 24}, {-22, 18}}, color = {0, 0, 255}, thickness = 1));
  connect(r_HTR_cout.outlet, r_heater_cin.inlet) annotation(
    Line(points = {{-22, 14}, {-48, 14}, {-48, 32}, {-54, 32}}, color = {0, 0, 255}, thickness = 1));
  connect(r_heater_cin.outlet, heater.gasIn) annotation(
    Line(points = {{-58, 32}, {-66, 32}}, color = {0, 0, 255}, thickness = 1));
  connect(heater.gasOut, r_heater_cout.inlet) annotation(
    Line(points = {{-82, 32}, {-88, 32}}, color = {0, 0, 255}, thickness = 1));
  connect(r_heater_cout.outlet, sink.flange) annotation(
    Line(points = {{-92, 32}, {-95, 32}, {-95, 33}, {-100, 33}}, color = {0, 0, 255}, thickness = 1));

  // hot stream for heater
  connect(source_heater_hot.flange, r_heater_hin.inlet) annotation(
    Line(points = {{-74, 6}, {-74, 12}}, color = {255, 0, 0}, thickness = 0.75));
  connect(r_heater_hin.outlet, heater.waterIn) annotation(
    Line(points = {{-74, 16}, {-74, 24}}, color = {255, 0, 0}, thickness = 0.75));
  connect(heater.waterOut, r_heater_hout.inlet) annotation(
    Line(points = {{-74, 40}, {-74, 50}}, color = {255, 0, 0}, thickness = 0.75));
  connect(r_heater_hout.outlet, sink_heater_hot.flange) annotation(
    Line(points = {{-74, 54}, {-74, 62}}, color = {255, 0, 0}, thickness = 0.75));

  // cold side of cooler
  connect(source_cooler.flange, r_cooler_cin.inlet) annotation(
    Line(points = {{56, -24}, {56, -28}}, color = {0, 255, 0}, thickness = 0.75));
  connect(r_cooler_cin.outlet, cooler.waterIn) annotation(
    Line(points = {{56, -32}, {56, -38}}, color = {0, 255, 0}, thickness = 0.75));
  connect(cooler.waterOut, r_cooler_cout.inlet) annotation(
    Line(points = {{56, -54}, {56, -62}}, color = {0, 255, 0}, thickness = 0.75));
  connect(r_cooler_cout.outlet, sink_cooler.flange) annotation(
    Line(points = {{56, -66}, {56, -72}, {57, -72}}, color = {0, 255, 0}, thickness = 0.75));

  // ramp input
  connect(ramp_p.y, source.in_p0) annotation(
    Line(points = {{-113, -32}, {-107.5, -32}, {-107.5, -36}, {-108, -36}}, color = {0, 0, 127}));
  connect(ramp_mdot.y, sink.in_w0) annotation(
    Line(points = {{-113.5, 21}, {-102, 21}, {-102, 30.5}}, color = {0, 0, 127}));

  annotation(
    Diagram(coordinateSystem(extent = {{-140, -100}, {140, 100}})),
    experiment(StartTime = 0, StopTime = 50, Tolerance = 1e-3, Interval = 1),
    __OpenModelica_commandLineOptions = "--matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection -d=initialization,NLSanalyticJacobian",
    // __OpenModelica_simulationFlags(lv = "LOG_DEBUG,LOG_NLS,LOG_NLS_V,LOG_STATS,LOG_INIT,LOG_STDOUT, -w", outputFormat = "mat", s = "dassl", nls = "homotopy"));
    __OpenModelica_simulationFlags(lv = "LOG_DEBUG,LOG_NLS,LOG_STATS,LOG_INIT,LOG_STDOUT -w", newtonFTol = "1e-6", newtonXTol= "1e-6", outputFormat = "mat", s = "dassl", nls = "homotopy"));
end TPDyn_RCBCycle;
