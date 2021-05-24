within Steps.Cycle.Transient;

model TP_RCBCycle

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
  parameter Model.HeatExchangerConfig cfg_heater = cfg.cfg_heater;
  parameter Model.HeatExchangerConfig cfg_LTR    = cfg.cfg_LTR;
  parameter Model.HeatExchangerConfig cfg_HTR    = cfg.cfg_HTR;
  parameter Model.TurbomachineryConfig cfg_turb  = cfg.cfg_turb;
  parameter Model.TurbomachineryConfig cfg_comp   = cfg.cfg_comp;
  parameter Model.TurbomachineryConfig cfg_recomp = cfg.cfg_recomp;
  parameter Model.HeatExchangerConfig cfg_cooler  = cfg.cfg_cooler;
  parameter Model.SplitterConfig cfg_splitter     = cfg.cfg_splitter;
  parameter Model.SplitterConfig cfg_merger      = cfg.cfg_merger;
  
  parameter Model.ThermoState st_bypass      = cfg.st_recomp_out;
  parameter Model.ThermoState st_source_temp = cfg_recomp.st_out;
  parameter Model.ThermoState st_sink_temp    = cfg_recomp.st_in;
  parameter Model.ThermoState st_source  = cfg_heater.cfg_cold.st_out;
  parameter Model.ThermoState st_sink   = cfg_heater.cfg_cold.st_out;  

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
    
    // DittusBoelter
    redeclare replaceable model HeatTransfer_F = ThermoPower.Thermal.HeatTransferFV.DittusBoelter,    
    redeclare replaceable model HeatTransfer_G = ThermoPower.Thermal.HeatTransferFV.DittusBoelter,
    fluidFlow(
      heatTransfer(heating=true), 
      fixedMassFlowSimplified = true,
      hstartin                = cfg_heater.cfg_hot.st_in.h,
      hstartout               = cfg_heater.cfg_hot.st_out.h),   // set the fluid flow as fixed mdot for simplarity
    gasFlow(
      heatTransfer(heating=false), 
      Tstartin  = cfg_heater.cfg_cold.st_in.T,
      Tstartout = cfg_heater.cfg_cold.st_out.T),
    
    
    // redeclare replaceable model HeatTransfer_F = ThermoPower.Thermal.HeatTransferFV.IdealHeatTransfer,    
    // redeclare replaceable model HeatTransfer_G = ThermoPower.Thermal.HeatTransferFVH.IdealHeatTransfer, 

    // redeclare replaceable model HeatTransfer_F = ThermoPower.Thermal.HeatTransferFV.ConstantHeatTransferCoefficient(),     
    // // redeclare replaceable model HeatTransfer_G = ThermoPower.Thermal.HeatTransferFV.IdealHeatTransfer,      
    // redeclare replaceable model HeatTransfer_G = ThermoPower.Thermal.HeatTransferFV.ConstantHeatTransferCoefficient(),      
    // fluidFlow(
    //   heatTransfer(gamma = cfg_heater.cfg_hot.gamma_HE),
    //   fixedMassFlowSimplified = true,
    //   hstartin                = cfg_heater.cfg_hot.st_in.h,
    //   hstartout               = cfg_heater.cfg_hot.st_out.h),   // set the fluid flow as fixed mdot for simplarity
    // gasFlow(
    //   heatTransfer(gamma = cfg_heater.cfg_cold.gamma_HE),
    //   Tstartin    = cfg_heater.cfg_cold.st_in.T,
    //   Tstartout   = cfg_heater.cfg_cold.st_out.T),     

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
  ThermoPower.Gas.SourceMassFlow source_heater_hot(
    redeclare package Medium = medium_main,
    T        = cfg_heater.cfg_hot.st_in.T,
    p0       = cfg_heater.cfg_hot.st_in.p,
    use_in_T = false,
    w0       = cfg_heater.cfg_hot.st_in.mdot,
    gas(
      p(start = cfg_heater.cfg_hot.st_in.p, nominal = cfg_heater.cfg_hot.st_in.p), 
      T(start = cfg_heater.cfg_hot.st_in.T, nominal = cfg_heater.cfg_hot.st_in.T))     
  );
  
  
  ThermoPower.Gas.SinkPressure sink_heater_hot(
    redeclare package Medium = medium_heater,
    T  = cfg_heater.cfg_hot.st_out.T,
    p0 = cfg_heater.cfg_hot.st_out.p) 
  annotation(
    Placement(transformation(extent = {{60, -10}, {80, 10}}, rotation = 0)));


  Steps.TPComponents.PCHE heater(
    redeclare package FluidMedium = medium_heater, 
    redeclare package FlueGasMedium = medium_main, 
     
    // use Marchionni PCHE HeatTransfer
    // slow but can have a result - set a_phi = 0 to use Gnielinski's correlation 
    
    redeclare replaceable model HeatTransfer_F = Steps.TPComponents.MarchionniPCHEHeatTransferFV(),
    redeclare replaceable model HeatTransfer_G = Steps.TPComponents.MarchionniPCHEHeatTransferFV(),
    gasFlow(
      avoidInletEnthalpyDerivative = false,
      heatTransfer(
        pitch = cfg_heater.cfg_hot.l_pitch,
        phi   = cfg_heater.cfg_hot.a_phi,
        Cf_C1 = Cf_C1,
        Cf_C2 = Cf_C2,
        Cf_C3 = Cf_C3
      )
    ),
    fluidFlow(
      avoidInletEnthalpyDerivative = true,
      heatTransfer(
        pitch = cfg_heater.cfg_cold.l_pitch,
        phi   = cfg_heater.cfg_cold.a_phi,
        Cf_C1 = Cf_C1,
        Cf_C2 = Cf_C2,
        Cf_C3 = Cf_C3
      )
    ),        

    // redeclare replaceable model HeatTransfer_F = ThermoPower.Thermal.HeatTransferFV.DittusBoelter,    
    // redeclare replaceable model HeatTransfer_G = ThermoPower.Thermal.HeatTransferFV.DittusBoelter,
    // gasFlow(heatTransfer(heating=true)),
    // fluidFlow(heatTransfer(heating=false)),   
    
    // fast and works fine for now. Error occurs when mass flow rate is zero, i.e. one flow is shut down. 
    // redeclare replaceable model HeatTransfer_F = ThermoPower.Thermal.HeatTransferFV.IdealHeatTransfer,      
    // redeclare replaceable model HeatTransfer_G = ThermoPower.Thermal.HeatTransferFV.IdealHeatTransfer,  
    cfg              = cfg_heater,
    N_G              = N_seg_heater,
    N_F              = N_seg_heater,
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
  );

  Steps.TPComponents.GasStateReader r_heater_hin(redeclare package Medium = medium_heater) annotation(
    Placement(visible = true, transformation(origin = {66, -56}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
  Steps.TPComponents.GasStateReader r_heater_hout(redeclare package Medium = medium_heater) annotation(
    Placement(visible = true, transformation(origin = {96, 44}, extent = {{-6, -6}, {6, 6}}, rotation = 180)));
  Steps.TPComponents.GasStateReader r_heater_cin(redeclare package Medium = medium_main) annotation(
    Placement(visible = true, transformation(origin = {-22, 50}, extent = {{-6, -6}, {6, 6}}, rotation = -90)));
  Steps.TPComponents.GasStateReader r_heater_cout(redeclare package Medium = medium_main) annotation(
    Placement(visible = true, transformation(origin = {-32, 20}, extent = {{-6, -6}, {6, 6}}, rotation = 180)));
*/

  ThermoPower.Gas.SourceMassFlow source(
    redeclare package Medium = medium_main, 
    T        = st_source.T,
    p0       = st_source.p,
    w0       = st_source.mdot,
    use_in_T = false,
    gas(
      p(start = st_source.p, nominal = st_source.p), 
      T(start = st_source.T, nominal = st_source.T)))
  annotation(
    Placement(transformation(extent = {{-70, -10}, {-50, 10}}, rotation = 0))); 

/*
  ThermoPower.Gas.SourcePressure source(
    redeclare package Medium = medium_main, 
    T        = st_source.T,
    p0       = st_source.p,
    use_in_T = true,
    gas(
      p(start = st_source.p, nominal = st_source.p), 
      T(start = st_source.T, nominal = st_source.T))) 
  annotation(
    Placement(transformation(origin = {0, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 270)));


  ThermoPower.Gas.SourcePressure source_temp(
    redeclare package Medium = medium_main, 
    T        = st_source_temp.T,
    p0       = st_source_temp.p,
    use_in_T = true,
    gas(
      p(start = st_source_temp.p, nominal = st_source_temp.p), 
      T(start = st_source_temp.T, nominal = st_source_temp.T))) 
  annotation(
    Placement(transformation(origin = {0, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 270)));
*/

/*
  ThermoPower.Gas.SourceMassFlow source_mixer_in(
    redeclare package Medium = medium_main,
    T        = st_bypass.T,
    p0       = st_bypass.p,
    use_in_T = false,
    w0       = st_bypass.mdot
  );  
*/
 
  ThermoPower.Gas.SinkPressure sink(
    redeclare package Medium = medium_main, 
    p0 = st_sink.p,
    T  = st_sink.T)
  annotation(
    Placement(transformation(origin = {0, -80}, extent = {{-10, -10}, {10, 10}}, rotation = 270)));
  /*
  ThermoPower.Gas.SourceMassFlow source_temp(
    redeclare package Medium = medium_main, 
    T        = st_source_temp.T,
    p0       = st_source_temp.p,
    use_in_T = false,
    w0       = st_source_temp.mdot)
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
  Gas.FlowJoin mixer(redeclare package Medium = medium_main);  

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
      T(nominal = turbine.Tstart_in)),
    gas_iso(
      p(nominal = turbine.pstart_out), 
      T(nominal = turbine.Tstart_out),
      h(start   = cfg_turb.st_out.h)))    
    annotation(
      Placement(transformation(extent = {{-40, -20}, {0, 20}}, rotation = 0)));

  Modelica.Mechanics.Rotational.Sources.ConstantSpeed const_speed_turb(
      w_fixed=cfg_turb.N, useSupport=false) annotation(
    Placement(visible = true, transformation(origin = {81, -7}, extent = {{-5, -5}, {5, 5}}, rotation = 0)));  
 
  ThermoPower.Gas.SensT sens_turbine(redeclare package Medium = medium_main) annotation(
    Placement(visible = true, transformation(origin = {20, 44}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    

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
       
    // // DittusBoelter
    // redeclare replaceable model HeatTransfer_F = ThermoPower.Thermal.HeatTransferFV.DittusBoelter,    
    // redeclare replaceable model HeatTransfer_G = ThermoPower.Thermal.HeatTransferFV.DittusBoelter,
    // fluidFlow(
    //   heatTransfer(heating=false), 
    //   fixedMassFlowSimplified = true,
    //   hstartin                = cfg_cooler.cfg_cold.st_in.h,
    //   hstartout               = cfg_cooler.cfg_cold.st_out.h),   // set the fluid flow as fixed mdot for simplarity
    // gasFlow(
    //   heatTransfer(heating=true), 
    //   Tstartin  = cfg_cooler.cfg_hot.st_in.T,
    //   Tstartout = cfg_cooler.cfg_hot.st_out.T),
    
    // redeclare replaceable model HeatTransfer_F = ThermoPower.Thermal.HeatTransferFV.IdealHeatTransfer,    
    // redeclare replaceable model HeatTransfer_G = ThermoPower.Thermal.HeatTransferFV.IdealHeatTransfer, 
    // fluidFlow(
    //   fixedMassFlowSimplified = true,
    //   hstartin                = cfg_cooler.cfg_cold.st_in.h,
    //   hstartout               = cfg_cooler.cfg_cold.st_out.h),   // set the fluid flow as fixed mdot for simplarity
    // gasFlow(
    //   Tstartin    = cfg_cooler.cfg_hot.st_in.T,
    //   Tstartout   = cfg_cooler.cfg_hot.st_out.T),
    
    redeclare replaceable model HeatTransfer_F = ThermoPower.Thermal.HeatTransferFV.ConstantHeatTransferCoefficient(),      
    redeclare replaceable model HeatTransfer_G = ThermoPower.Thermal.HeatTransferFV.ConstantHeatTransferCoefficient(), 
        
    fluidFlow(
      heatTransfer(gamma = cfg_cooler.cfg_cold.gamma_HE),
      fixedMassFlowSimplified = true,
      hstartin                = cfg_cooler.cfg_cold.st_in.h,
      hstartout               = cfg_cooler.cfg_cold.st_out.h),   // set the fluid flow as fixed mdot for simplarity
    gasFlow(
      heatTransfer(gamma = cfg_cooler.cfg_hot.gamma_HE),
      Tstartin    = cfg_cooler.cfg_hot.st_in.T,
      Tstartout   = cfg_cooler.cfg_hot.st_out.T),   
           
    redeclare model HeatExchangerTopology = ThermoPower.Thermal.HeatExchangerTopologies.CounterCurrentFlow,
    cfg             = cfg_cooler,
    N_G             = N_seg_cooler,
    N_F             = N_seg_cooler,
    SSInit          = true,
    gasQuasiStatic  = true,
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
            
  Modelica.Mechanics.Rotational.Sources.ConstantSpeed const_speed_mc(
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
            
  Modelica.Mechanics.Rotational.Sources.ConstantSpeed const_speed_rc(
    w_fixed    = cfg_recomp.N,
    useSupport = false)
  annotation (Placement(transformation(
      extent={{-50,-10},{-30,10}}, rotation=0)));

  
  ThermoPower.Gas.FlowSplit splitter(redeclare package Medium = medium_main); 

  /*
  ThermoPower.Gas.Utility.ClosedSystemInit sys_init(
    redeclare package Medium = medium_main,
    //w_b(nominal = bc_heater.st_cold_out.mdot),
    pstart = bc_heater.st_cold_out.p);
  */
  /*
  ThermoPower.Gas.ThroughMassFlow source_through_cooler(
    redeclare package Medium = medium_main,
    w0 = bc_cooler.st_hot_out.mdot) annotation(
    Placement(visible = true, transformation(origin = {88, 12}, extent = {{-6, -6}, {6, 6}}, rotation = 90)));

  
  ThermoPower.Gas.ThroughMassFlow source_through_heater(
    redeclare package Medium = medium_main,
    w0 = bc_heater.st_cold_out.mdot);
*/ 


  // calculated variables
  /*
  // for turbine
  Modelica.SIunits.Power W_turb = (r05.h - r06.h) * r05.w / 1e6 "W->MW, net power for turbine";
  Modelica.SIunits.Efficiency eta_turb = turbine.eta * 100;
  // main compressor
  Modelica.SIunits.Power W_MC = (r02.h - r01.h) * r01.w / 1e6 "W->MW, net power for main compressor";
  Modelica.SIunits.Efficiency eta_MC = compressor.eta * 100;
  // re compressor
  Modelica.SIunits.Power W_RC = (r09.h - r08b.h) * r08b.w / 1e6 "W->MW, net power for recompressor";
  Modelica.SIunits.Efficiency eta_RC = recompressor.eta * 100;
  
  // power block
  Modelica.SIunits.Power W_net = W_turb - W_MC - W_RC "net power generated";
  Modelica.SIunits.Power Q_heater = (rh1.h - rh2.h) * rh1.w / 1e6 "heat input for the power block";
  Modelica.SIunits.Power Q_cooler = (rc2.h - rc1.h) * rc1.w / 1e6 "heat input for the power block";
  Modelica.SIunits.Efficiency eta_pb = W_net / Q_heater * 100 "power block efficiency";
  Real SR = r08b.w / r08.w * 100 "split ratio";

  // heat transfer coefficient for HTR and LTR
  // HTR
  Modelica.SIunits.Power Q_HTR = (r04.h - r03.h) * r03.w / 1e6 "W->MW, heat input for HTR";
  Modelica.SIunits.TemperatureDifference dT1_HTR = (r06.T - r04.T);
  Modelica.SIunits.TemperatureDifference dT2_HTR = (r07.T - r03.T);
  Real T_ltmd_HTR = if dT1_HTR  > 0 and dT2_HTR > 0 then (dT2_HTR - dT1_HTR) / Modelica.Math.log(abs(dT2_HTR / dT1_HTR)) else -1;
  Real UA_HTR = if T_ltmd_HTR > 0 then Q_HTR / T_ltmd_HTR else 0.0;
  // LTR
  Modelica.SIunits.Power Q_LTR = (r10.h - r02.h) * r02.w / 1e6 "W->MW";
  Modelica.SIunits.TemperatureDifference dT1_LTR = (r07.T - r10.T);
  Modelica.SIunits.TemperatureDifference dT2_LTR = (r08.T - r02.T);
  Real T_ltmd_LTR = if dT1_LTR > 0 and dT2_LTR > 0 then (dT2_LTR - dT1_LTR) / Modelica.Math.log(abs(dT2_LTR / dT1_LTR)) else -1;
  Real UA_LTR = if T_ltmd_LTR > 0 then Q_LTR / T_ltmd_LTR else 0.0;  
  
  // Liquid Na exit temperature
  Modelica.SIunits.Temperature T_heater_hot_out = rh2.T;
  */
  
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
  Modelica.Blocks.Sources.IntegerConstant const_T_offset_c(k = integer(st_source.T)) annotation(
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
  // performance map for main compressor
  parameter Real tableEta_comp_mc[5, 4]  = [0, 95, 100, 105; 1, 0.85310219, 0.837591241, 0.832420925; 2, 0.868613139, 0.857238443, 0.851034063; 3, 0.860340633, 0.85, 0.842761557; 4, 0.85310219, 0.839659367, 0.816909976];
  parameter Real tablePhic_comp_mc[5, 4] = [0, 95, 100, 105; 1, 0.000134346, 0.000150832, 0.000164161; 2, 0.000137854, 0.000153638, 0.00016802; 3, 0.000142414, 0.000158549, 0.000169774; 4, 0.000145921, 0.000161706, 0.000171528];
  parameter Real tablePR_comp_mc[5, 4]   = [0, 95, 100, 105; 1, 1.967529638, 2.350588505, 2.785882673; 2, 1.915294338, 2.315764972, 2.681412073; 3, 1.810823737, 2.220000255, 2.524706172; 4, 1.654117837, 2.115529655, 2.359294389];

  // performance map for re compressor      
  parameter Real tableEta_comp_rc[5, 4]  = [0, 95, 100, 105; 1, 0.85310219, 0.837591241, 0.832420925; 2, 0.868613139, 0.857238443, 0.851034063; 3, 0.860340633, 0.85, 0.842761557; 4, 0.85310219, 0.839659367, 0.816909976];
  parameter Real tablePhic_comp_rc[5, 4] = [0, 95, 100, 105; 1, 7.17663E-05, 8.05731E-05, 8.76935E-05; 2, 7.36401E-05, 8.20721E-05, 8.97547E-05; 3, 7.6076E-05, 8.46954E-05, 9.06916E-05; 4, 7.79498E-05, 8.63819E-05, 9.16285E-05];
  parameter Real tablePR_comp_rc[5, 4]   = [0, 95, 100, 105; 1, 1.967529638, 2.350588505, 2.785882673; 2, 1.915294338, 2.315764972, 2.681412073; 3, 1.810823737, 2.220000255, 2.524706172; 4, 1.654117837, 2.115529655, 2.359294389];
   
equation
/*
  // HTR alone
  connect(source_temp.flange, HTR.waterIn);

  connect(HTR.waterOut, T_waterOut.inlet) annotation(
    Line(points = {{8.88178e-016, -44}, {8.88178e-016, -20}, {0, -20}}, thickness = 0.5, color = {0, 0, 255}));
  connect(T_waterOut.outlet, sink.flange) annotation(
    Line(points = {{1.83697e-015, -70}, {1.83697e-015, -56}, {-8.88178e-016, -56}}, thickness = 0.5, color = {0, 0, 255}));

  // connect(HTR.waterOut, sink.flange) 
   //annotation(Line(points = {{1.83697e-015, -70}, {1.83697e-015, -56}, {-8.88178e-016, -56}}, thickness = 0.5, color = {0, 0, 255}));   
  
  connect(source.flange, HTR.gasIn) annotation(
    Line(points = {{-50, 0}, {-20, 0}}, color = {159, 159, 223}, thickness = 0.5, smooth = Smooth.None));   
    
  connect(HTR.gasOut, T_gasOut.inlet ) annotation(
    Line(points = {{34, 0}, {34, 0}, {20, 0}}, color = {159, 159, 223}, thickness = 0.5));
  connect(T_gasOut.outlet, sink_temp.flange) annotation(
    Line(points = {{46, 0}, {46, 0}, {60, 0}}, color = {159, 159, 223}, thickness = 0.5));  

  // LTR alone
  connect(source_temp.flange, LTR.waterIn);

  connect(LTR.waterOut, T_waterOut.inlet) annotation(
    Line(points = {{8.88178e-016, -44}, {8.88178e-016, -20}, {0, -20}}, thickness = 0.5, color = {0, 0, 255}));
  connect(T_waterOut.outlet, sink.flange) annotation(
    Line(points = {{1.83697e-015, -70}, {1.83697e-015, -56}, {-8.88178e-016, -56}}, thickness = 0.5, color = {0, 0, 255}));

  // connect(HTR.waterOut, sink.flange) 
   //annotation(Line(points = {{1.83697e-015, -70}, {1.83697e-015, -56}, {-8.88178e-016, -56}}, thickness = 0.5, color = {0, 0, 255}));   
  
  connect(source.flange, LTR.gasIn) annotation(
    Line(points = {{-50, 0}, {-20, 0}}, color = {159, 159, 223}, thickness = 0.5, smooth = Smooth.None));   
    
  connect(LTR.gasOut, T_gasOut.inlet ) annotation(
    Line(points = {{34, 0}, {34, 0}, {20, 0}}, color = {159, 159, 223}, thickness = 0.5));
  connect(T_gasOut.outlet, sink_temp.flange) annotation(
    Line(points = {{46, 0}, {46, 0}, {60, 0}}, color = {159, 159, 223}, thickness = 0.5));  
*/
   
 
/*
  //HTR + mixer

  connect(source_mixer_in.flange, mixer.inlet1);  
  connect(source_temp.flange, mixer.inlet2);  
  connect(mixer.outlet, HTR.waterIn);

  connect(HTR.waterOut, T_waterOut.inlet) annotation(
    Line(points = {{8.88178e-016, -44}, {8.88178e-016, -20}, {0, -20}}, thickness = 0.5, color = {0, 0, 255}));
  connect(sink.flange, T_waterOut.outlet) annotation(
    Line(points = {{1.83697e-015, -70}, {1.83697e-015, -56}, {-8.88178e-016, -56}}, thickness = 0.5, color = {0, 0, 255}));

  // connect(HTR.waterOut, sink.flange) 
   //annotation(Line(points = {{1.83697e-015, -70}, {1.83697e-015, -56}, {-8.88178e-016, -56}}, thickness = 0.5, color = {0, 0, 255}));   
  
  connect(source.flange, HTR.gasIn) annotation(
    Line(points = {{-50, 0}, {-20, 0}}, color = {159, 159, 223}, thickness = 0.5, smooth = Smooth.None));   
    
  connect(HTR.gasOut, T_gasOut.inlet ) annotation(
    Line(points = {{34, 0}, {34, 0}, {20, 0}}, color = {159, 159, 223}, thickness = 0.5));
  connect(T_gasOut.outlet, sink_temp.flange) annotation(
    Line(points = {{46, 0}, {46, 0}, {60, 0}}, color = {159, 159, 223}, thickness = 0.5));
  
  // connect(HTR.gasOut, sink_temp.flange) annotation(
  //   Line(points = {{46, 0}, {46, 0}, {60, 0}}, color = {159, 159, 223}, thickness = 0.5));
*/  
/*
  // mixer + LTR
  // water/cold side  
  connect(source_mixer_in.flange, mixer.inlet1);
  
  connect(source_temp.flange, LTR.waterIn);
  
  connect(LTR.waterOut, mixer.inlet2);
  
  connect(mixer.outlet, T_waterOut.inlet);
  
  connect(T_waterOut.outlet, sink.flange) annotation(
    Line(points = {{1.83697e-015, -70}, {1.83697e-015, -56}, {-8.88178e-016, -56}}, thickness = 0.5, color = {0, 0, 255}));
  
  // gas/hot side
  connect(source.flange, LTR.gasIn);
  
  connect(LTR.gasOut, T_gasOut.inlet);
  connect(T_gasOut.outlet, sink_temp.flange) annotation(
    Line(points = {{46, 0}, {46, 0}, {60, 0}}, color = {159, 159, 223}, thickness = 0.5));
*/

/*
  //HTR + mixer + LTR 
  // water/cold side    
  connect(source_temp.flange, r_LTR_cin.inlet);  
  connect(r_LTR_cin.outlet, LTR.waterIn);  
  connect(source_mixer_in.flange, mixer.inlet1);  
  connect(LTR.waterOut, r_LTR_cout.inlet);
  connect(r_LTR_cout.outlet, mixer.inlet2);  
  connect(mixer.outlet, r_HTR_cin.inlet);
  connect(r_HTR_cin.outlet, HTR.waterIn);
  connect(HTR.waterOut, r_HTR_cout.inlet);  
  connect(r_HTR_cout.outlet, sink.flange) annotation(
    Line(points = {{1.83697e-015, -70}, {1.83697e-015, -56}, {-8.88178e-016, -56}}, thickness = 0.5, color = {0, 0, 255}));
  
  // gas/hot side
  connect(source.flange, r_HTR_hin.inlet);
  connect(r_HTR_hin.outlet, HTR.gasIn) annotation(
   Line(points = {{-50, 0}, {-20, 0}}, color = {159, 159, 223}, thickness = 0.5, smooth = Smooth.None));
  connect(HTR.gasOut, r_HTR_hout.inlet);
  connect(r_HTR_hout.outlet, r_LTR_hin.inlet);
  connect(r_LTR_hin.outlet, LTR.gasIn);  
  connect(LTR.gasOut, r_LTR_hout.inlet);  
  connect(r_LTR_hout.outlet, sink_temp.flange) annotation(
    Line(points = {{46, 0}, {46, 0}, {60, 0}}, color = {159, 159, 223}, thickness = 0.5));
*/

/*
  //HTR + mixer + LTR + Heater
  // main stream, water/cold side  
  connect(source_mixer_in.flange, mixer.inlet1);  
  connect(source_temp.flange, r_LTR_cin.inlet);
  connect(r_LTR_cin.outlet, LTR.waterIn);  
  connect(LTR.waterOut, r_LTR_cout.inlet);
  connect(r_LTR_cout.outlet, mixer.inlet2);  
  connect(mixer.outlet, r_HTR_cin.inlet);
  connect(r_HTR_cin.outlet, HTR.waterIn);
  connect(HTR.waterOut, r_HTR_cout.inlet);
  connect(r_HTR_cout.outlet, r_heater_cin.inlet);
  connect(r_heater_cin.outlet, heater.gasIn);  
  connect(heater.gasOut, r_heater_cout.inlet);
  connect(r_heater_cout.outlet, sink.flange) annotation(
    Line(points = {{1.83697e-015, -70}, {1.83697e-015, -56}, {-8.88178e-016, -56}}, thickness = 0.5, color = {0, 0, 255}));
  
  // main stream, gas/hot side
  connect(source.flange, r_HTR_hin.inlet);
  connect(r_HTR_hin.outlet, HTR.gasIn) annotation(
   Line(points = {{-50, 0}, {-20, 0}}, color = {159, 159, 223}, thickness = 0.5, smooth = Smooth.None));
  connect(HTR.gasOut, r_HTR_hout.inlet);
  connect(r_HTR_hout.outlet, r_LTR_hin.inlet);
  connect(r_LTR_hin.outlet,LTR.gasIn);  
  connect(LTR.gasOut, r_LTR_hout.inlet);
  connect(r_LTR_hout.outlet, sink_temp.flange) annotation(
    Line(points = {{46, 0}, {46, 0}, {60, 0}}, color = {159, 159, 223}, thickness = 0.5));
   
  // hot stream for heater
  connect(source_heater_hot.flange, r_heater_hin.inlet);
  connect(r_heater_hin.outlet, heater.waterIn);
  connect(heater.waterOut, r_heater_hout.inlet);
  connect(r_heater_hout.outlet, sink_heater_hot.flange);
*/


  // HTR + mixer + LTR + Heater + Turbine - OPEN LOOP version 2, different with one more pair of source and sink, same layout as RCBC's left part. 
  // main stream, water/cold side  
  
  // main stream, gas/hot side
  connect(source.flange, turbine.inlet);
  
    connect(turbine.shaft_b, const_speed_turb.flange) annotation(
      Line(points = {{30, 0}, {74, 0}, {74, 0}, {74, 0}}));
      
  connect(turbine.outlet, sens_turbine.inlet) 
  annotation(
   Line(points = {{-50, 0}, {-20, 0}}, color = {159, 159, 223}, thickness = 0.5, smooth = Smooth.None));   
  connect(sens_turbine.outlet, r_HTR_hin.inlet);
  connect(r_HTR_hin.outlet, HTR.gasIn) annotation(
   Line(points = {{-50, 0}, {-20, 0}}, color = {159, 159, 223}, thickness = 0.5, smooth = Smooth.None));
  connect(HTR.gasOut, r_HTR_hout.inlet);
  connect(r_HTR_hout.outlet, r_LTR_hin.inlet);
  connect(r_LTR_hin.outlet,LTR.gasIn);  
  connect(LTR.gasOut, r_LTR_hout.inlet);
  connect(r_LTR_hout.outlet, splitter.inlet); 
  
  connect(splitter.outlet1, r_cooler_hin.inlet) annotation(
    Line(points = {{46, 0}, {46, 0}, {60, 0}}, color = {159, 159, 223}, thickness = 0.5));
    
    connect(splitter.outlet2, recompressor.inlet);
    
      connect(recompressor.shaft_a, const_speed_rc.flange);
    
    connect(recompressor.outlet, mixer.inlet2);    
    
  connect(r_cooler_hin.outlet, cooler.gasIn);  
  connect(cooler.gasOut, r_cooler_hout.inlet);  
  connect(r_cooler_hout.outlet, compressor.inlet);
  
    connect(compressor.shaft_a, const_speed_mc.flange);     
  
  connect(compressor.outlet, r_LTR_cin.inlet);
  connect(r_LTR_cin.outlet, LTR.waterIn);  
  connect(LTR.waterOut, r_LTR_cout.inlet);
  connect(r_LTR_cout.outlet, mixer.inlet1);  
  connect(mixer.outlet, r_HTR_cin.inlet);
  connect(r_HTR_cin.outlet, HTR.waterIn);
  connect(HTR.waterOut, r_HTR_cout.inlet);
  connect(r_HTR_cout.outlet, r_heater_cin.inlet);
  connect(r_heater_cin.outlet, heater.gasIn);  
  connect(heater.gasOut, r_heater_cout.inlet);
  connect(r_heater_cout.outlet, sink.flange) annotation(
    Line(points = {{1.83697e-015, -70}, {1.83697e-015, -56}, {-8.88178e-016, -56}}, thickness = 0.5, color = {0, 0, 255}));   
    
  // connect(source_temp.flange, mixer.inlet2);
   
  // hot stream for heater
  connect(source_heater_hot.flange, r_heater_hin.inlet);
  connect(r_heater_hin.outlet, heater.waterIn);
  connect(heater.waterOut, r_heater_hout.inlet);
  connect(r_heater_hout.outlet, sink_heater_hot.flange);
  
  // cold side
  connect(source_cooler.flange, r_cooler_cin.inlet);
  connect(r_cooler_cin.outlet, cooler.waterIn);
  connect(cooler.waterOut, r_cooler_cout.inlet);
  connect(r_cooler_cout.outlet, sink_cooler.flange);
    


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
  connect(I2R_T_h.y, source.in_T) annotation(
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
  connect(I2R_T_c.y, source.in_T) annotation(
    Line(points = {{-37, 124}, {6, 124}, {6, 64}}, color = {0, 0, 127}));


  annotation(
    Diagram(coordinateSystem(extent = {{-100, -100}, {120, 100}})),
    experiment(StartTime = 0, StopTime = 20, Tolerance = 1e-3, Interval = 2),
    __OpenModelica_commandLineOptions = "--matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection -d=initialization,NLSanalyticJacobian",
    __OpenModelica_simulationFlags(lv = "LOG_DEBUG,LOG_NLS,LOG_NLS_V,LOG_STATS,LOG_INIT,LOG_STDOUT, -w", outputFormat = "mat", s = "dassl", nls = "homotopy"));
end TP_RCBCycle;
