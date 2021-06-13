within Steps.Test.TPComponents;

model TestTP_TurbineFixedP_ramp_mdot
  "Test for Turbine maintaining fiexd inlet p in ThermoPower"  
  import Modelica.SIunits.Conversions.{from_degC, from_deg};
  import Modelica.SIunits.{Temperature, Pressure, SpecificEnthalpy};
  import Util = Utilities.Util;
  import Steps.Utilities.CoolProp.PropsSI;  
  import Steps.Components.{PCHEGeoParam};
  import Steps.Model.PBConfiguration;
  import Steps.Model.{PBConfiguration, SimParam, EntityConfig, EntityGeoParam, EntityThermoParam, ThermoState, HEBoundaryCondition} ;
  import ThermoPower.Choices.Init.Options;
  import ThermoPower.System;
  import ThermoPower.Gas;
  
  // package Medium = Steps.Media.CO2; 
  // package Medium = ExternalMedia.Examples.CO2CoolProp;
  package Medium = Steps.Media.SCO2; 
  // package Medium = Modelica.Media.IdealGases.SingleGases.CO2;
  
  parameter Model.RCBCycleConfig cfg(
    // mdot_heater      = 40,
    // T_heater_hot_in  = from_degC(800),
    // T_heater_hot_out = from_degC(600),
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
  /*
  // for main compressor test
  parameter Model.TurbomachineryConfig cfg_comp   = cfg.cfg_comp;
  parameter Real tableEta[:,:]  = tableEta_mc;
  parameter Real tablePhic[:,:] = tablePhic_mc;
  parameter Real tablePR[:, :]  = tablePR_mc;  
  */
  
  // for re-compressor test
  parameter Model.TurbomachineryConfig cfg_comp   = cfg.cfg_recomp;
  parameter Model.TurbomachineryConfig cfg_turb   = cfg.cfg_turb;
  parameter Real tableEta[:,:]  = tableEta_rc;
  parameter Real tablePhic[:,:] = tablePhic_rc;
  parameter Real tablePR[:, :]  = tablePR_rc;    
  
  parameter Model.ThermoState st_source           = cfg_turb.st_in;
  parameter Model.ThermoState st_sink             = cfg_comp.st_out;
  
  //configuration for main compressor          
  parameter Modelica.SIunits.Temperature T_comp_des = from_degC(45) "Design temperature for main compressor";
  
  /* 
  ThermoPower.Gas.SourcePressure SourceP1(
    redeclare package Medium = Medium, 
    T        = st_source.T,
    p0       = st_source.p,
    use_in_T = false,
    use_in_p0 = true,
    gas(
      p(start = st_source.p, nominal = st_source.p), 
      T(start = st_source.T, nominal = st_source.T))) 
  annotation(
    Placement(visible = true, transformation(origin = {-105, -39}, extent = {{-5, -5}, {5, 5}}, rotation = 0)));
  */


  ThermoPower.Gas.SourceMassFlow SourceP1(
    redeclare package Medium = Medium,
    p0        = st_source.p,
    T         = st_source.T,
    w0        = st_source.mdot,
    use_in_w0 = true,
    gas(
      p(nominal = st_source.p), 
      T(nominal=st_source.T))) 
      annotation (Placement(transformation(extent={{-80,6},{-60,26}},rotation=0)));

  
 
  ThermoPower.Gas.SinkPressure SinkP1(
    redeclare package Medium = Medium,
    p0=st_sink.p,
    T=st_sink.T,
    use_in_p0 = true,
    gas(
      p(nominal = st_sink.p), 
      T(nominal = st_sink.T))) 
      annotation (Placement(transformation(extent={{40,6},{60,26}}, rotation=0)));

/*
  ThermoPower.Gas.SinkMassFlow SinkP1(
    redeclare package Medium = Medium, 
    T        = st_sink.T,
    p0       = st_sink.p,
    use_in_T = false,
    use_in_w0 = false,
    w0       = st_sink.mdot
  ) 
  annotation(
    Placement(visible = true, transformation(origin = {-105, 33}, extent = {{-5, -5}, {5, 5}}, rotation = 180)));
*/
  
  ThermoPower.Gas.Turbine turbine(
    redeclare package Medium = Medium, 
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
 
  Steps.TPComponents.GasStateReader r_turb_in(redeclare package Medium = Medium) annotation(
    Placement(visible = true, transformation(origin = {-92, -38}, extent = {{-4, -4}, {4, 4}}, rotation = 0)));

  Steps.TPComponents.TurbineFixedPController FPC_turb(
    redeclare package Medium   = Medium,
    fileName  = Modelica.Utilities.Files.loadResource("modelica://Steps/Resources/Data/turbine_map_2.txt"),
    tablePhic = fill(0.0, 14, 12), 
    Table     = ThermoPower.Choices.TurboMachinery.TableTypes.file,
    Ndesign   = cfg_turb.N,
    Tdes_in   = cfg_turb.st_in.T,
    use_in_p0 = true,    
    T_inlet(start = cfg_turb.st_in.T, nominal = cfg_turb.st_in.T),
    p_inlet(start = cfg_turb.st_in.p, nominal = cfg_turb.st_in.p),
    mdot_inlet(start = cfg_turb.st_in.mdot, nominal = cfg_turb.st_in.mdot)
  ) "Fixed pressure controller for main compressor";
  
  Modelica.Mechanics.Rotational.Sources.ConstantSpeed const_speed_turb(
      w_fixed=cfg_turb.N, useSupport=false) annotation(
    Placement(visible = true, transformation(origin = {-49, -47}, extent = {{5, -5}, {-5, 5}}, rotation = 0)));
   
  
  Modelica.Mechanics.Rotational.Sources.Speed speed_turb(
    exact = false,
    w_ref(start = cfg_turb.N, nominal= cfg_turb.N),
    w(start = cfg_turb.N, nominal = cfg_turb.N)
  ); 
  
  /*
  ThermoPower.Gas.Compressor Compressor(
    redeclare package Medium = Medium,
    pstart_in                  = cfg_comp.st_in.p,
    pstart_out                 = cfg_comp.st_out.p,    
    Tstart_in                  = cfg_comp.st_in.T,
    Tstart_out                 = cfg_comp.st_out.T,
    tablePhic                  = tablePhic,
    tableEta                   = tableEta,
    tablePR                    = tablePR,
    Table                      = ThermoPower.Choices.TurboMachinery.TableTypes.matrix,
    Ndesign                    = cfg_comp.N,
    Tdes_in                    = cfg_comp.st_in.T,
    explicitIsentropicEnthalpy = false,    
    gas_in(
      p(nominal = cfg_comp.st_in.p), 
      T(nominal = cfg_comp.st_in.T)),
    gas_iso(
      p(nominal = cfg_comp.st_out.p), 
      T(nominal = cfg_comp.st_out.T),
      h(start   = cfg_comp.st_out.h))) annotation (Placement(transformation(extent={{-20,-20},{
              20,20}}, rotation=0)));  


  Modelica.Mechanics.Rotational.Sources.ConstantSpeed ConstantSpeed1(
      w_fixed=cfg_comp.N, useSupport=false) annotation (Placement(transformation(
          extent={{-50,-10},{-30,10}}, rotation=0)));
  */
  
  inner ThermoPower.System system
    annotation (Placement(transformation(extent={{80,80},{100,100}})));

  // value to accelerate the simulation
  constant Integer SEC_PER_HOUR = integer(60 * 60 / time_scaling_factor); 
  constant Real time_scaling_factor = 12; // 5 min = 1 hour
  // constant Real time_scaling_factor = 1; // 1 hour = 1 hour  

  // ramp input to simulate the ramp change in ASPEN+ simulation
  // ramp change for case I in Guan's report
  Steps.TPComponents.RampSeq_W ramp_mdot(    
    time_start = 3 * SEC_PER_HOUR,
    interval   = 1 * SEC_PER_HOUR,
    duration   = 0.15 * SEC_PER_HOUR,
    offset     = st_source.mdot
  );
  
  Modelica.Blocks.Sources.Constant p_constant_source(k(start = st_source.p, nominal = st_source.p));
  
  Modelica.Blocks.Sources.Constant p_constant_sink(k(start = st_sink.p, nominal = st_sink.p));
  
protected
  // performance map for main compressor
  parameter Real tableEta_mc[5, 4]  = [0, 95, 100, 105; 1, 0.85310219, 0.837591241, 0.832420925; 2, 0.868613139, 0.857238443, 0.851034063; 3, 0.860340633, 0.85, 0.842761557; 4, 0.85310219, 0.839659367, 0.816909976];
  parameter Real tablePhic_mc[5, 4] = [0, 95, 100, 105; 1, 0.000134346, 0.000150832, 0.000164161; 2, 0.000137854, 0.000153638, 0.00016802; 3, 0.000142414, 0.000158549, 0.000169774; 4, 0.000145921, 0.000161706, 0.000171528];
  parameter Real tablePR_mc[5, 4]   = [0, 95, 100, 105; 1, 1.967529638, 2.350588505, 2.785882673; 2, 1.915294338, 2.315764972, 2.681412073; 3, 1.810823737, 2.220000255, 2.524706172; 4, 1.654117837, 2.115529655, 2.359294389];

  // performance map for re compressor      
  parameter Real tableEta_rc[5, 4]  = [0, 95, 100, 105; 1, 0.85310219, 0.837591241, 0.832420925; 2, 0.868613139, 0.857238443, 0.851034063; 3, 0.860340633, 0.85, 0.842761557; 4, 0.85310219, 0.839659367, 0.816909976];
  parameter Real tablePhic_rc[5, 4] = [0, 95, 100, 105; 1, 7.17663E-05, 8.05731E-05, 8.76935E-05; 2, 7.36401E-05, 8.20721E-05, 8.97547E-05; 3, 7.6076E-05, 8.46954E-05, 9.06916E-05; 4, 7.79498E-05, 8.63819E-05, 9.16285E-05];
  parameter Real tablePR_rc[5, 4]   = [0, 95, 100, 105; 1, 1.967529638, 2.350588505, 2.785882673; 2, 1.915294338, 2.315764972, 2.681412073; 3, 1.810823737, 2.220000255, 2.524706172; 4, 1.654117837, 2.115529655, 2.359294389];
  
equation
  connect(SourceP1.flange, r_turb_in.inlet);
  connect(r_turb_in.outlet, turbine.inlet);
  
    connect(r_turb_in.T, FPC_turb.T_inlet);
    connect(r_turb_in.w, FPC_turb.mdot_inlet);
    connect(p_constant_source.y, FPC_turb.p_inlet);
    
    connect(FPC_turb.omega, speed_turb.w_ref);
    // connect(const_speed_turb.flange, turbine.shaft_a);
    connect(speed_turb.flange, turbine.shaft_a);
    
  connect(turbine.outlet, SinkP1.flange) annotation (Line(
      points={{16,16},{40,16}},
      color={159,159,223},
      thickness=0.5));
    
  connect(p_constant_sink.y, SinkP1.in_p0);
  
  connect(p_constant_sink.y, FPC_turb.in_p0);    
  
  connect(ramp_mdot.y, SourceP1.in_w0);
 
  /*
  connect(ConstantSpeed1.flange, Compressor.shaft_a) annotation (Line(
      points={{-30,0},{-30,0},{-26,-0.2},{-12,0}},
      color={0,0,0},
      thickness=0.5));
  */
  
  annotation (
    Diagram(graphics),
    experiment(StartTime = 0, StopTime = 2100, Tolerance = 1e-2, Interval = 100),
    __OpenModelica_commandLineOptions = "--matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection -d=initialization,NLSanalyticJacobian,aliasConflicts",    
    __OpenModelica_simulationFlags(lv = "LOG_DEBUG,LOG_NLS,LOG_NLS_V,LOG_STATS,LOG_INIT,LOG_STDOUT, -w", outputFormat = "mat", s = "dassl", nls = "homotopy"));
end TestTP_TurbineFixedP_ramp_mdot;
