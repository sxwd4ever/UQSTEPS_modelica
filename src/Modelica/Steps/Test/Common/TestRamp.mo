within Steps.Test.Common;

model TestRamp

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

  package medium_main = Steps.Media.SCO2(
    // inputChoice = ExternalMedia.Common.InputChoice.pT,
    // substanceNames = {"CO2|debug=40"}
  );
  
  // package medium_main = ExternalMedia.Examples.CO2CoolProp;
  // package medium_heater = Steps.Media.ThermiaOilD; // out of working range of this 10Mw high T loop
  package medium_heater = SolarTherm.Media.Sodium.Sodium_pT;
  package medium_cooler = ThermoPower.Water.StandardWater;//Modelica.Media.IdealGases.SingleGases.CO2; 
  
  parameter Real table_k_metalwall[:,:] = [293.15, 12.1; 373.15, 16.3; 773.15, 21.5]; 

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
    
  constant Integer SEC_PER_HOUR = integer(60 * 60 / time_scaling_factor); 
  // constant Real time_scaling_factor = 12; // 5 min = 1 hour
  constant Real time_scaling_factor = 1; // 1 hour = 1 hour


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
equation
  
  
  annotation(
    experiment(StartTime = 0, StopTime = 14400, Interval = 10, Tolerance = 1e-6),
    __OpenModelica_commandLineOptions = "--matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection -d=initialization,NLSanalyticJacobian,aliasConflicts");
end TestRamp;
