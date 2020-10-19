within Steps.Cycle;

model OffDesignRCBCycle_v2
  "OffDesign Brayton Cycle with Recuperator + Recompression"  
// Adjustable parameters
  import Util = Utilities.Util;
  import Modelica.SIunits.Conversions.{from_degC, from_deg};
  import Steps.Utilities.CoolProp.PropsSI;
  import Steps.Utilities.Util.NewThermoState_pT;
  import Steps.Components.PCHEBoundaryCondition;
  import Steps.Components.PCHEGeoParam;
  import Modelica.SIunits.{Temperature, Pressure, SpecificEnthalpy};
  
  replaceable package PBMedia = Steps.Media.SCO2; 
 
  // efficiency of main compressor, bypass_compressor and turbine
  parameter Real eta_main_compressor = 0.89;
  
  parameter Real eta_bypass_compressor = 0.89;
  
  parameter Real eta_turbine = 0.9;
  
  // mass split ratio of splitter
  parameter Real splitter_split_ratio = mdot_bypass/mdot_main; 

  parameter Modelica.SIunits.Pressure p_pump_in = 8e6;
  parameter Modelica.SIunits.Pressure p_pump_out = 20e6;
  
  parameter Modelica.SIunits.MassFlowRate mdot_main = 51.91;
  parameter Modelica.SIunits.MassFlowRate mdot_pump = 31.31;
  parameter Modelica.SIunits.MassFlowRate mdot_bypass = mdot_main - mdot_pump;  
  
  parameter Modelica.SIunits.Temperature T_HTR_hot_out = from_degC(156.45);
  parameter Modelica.SIunits.Temperature T_LTR_cold_out = from_degC(151.45);
  
  
  parameter PCHEGeoParam geo_HTR(
    // pitch length, m
    pitch = 12e-3,
    // pitch angle
    phi = from_deg((180 - 108) /2),
    // length of pche, m
    length = 2860e-3,
    // Diameter of semi_circular, m
    d_c = 2e-3,
    // number of channels
    N_ch = integer(94e3),
    // number of segments
    N_seg = 50);
    
  parameter PCHEGeoParam geo_LTR(
    // pitch length, m
    pitch = 12e-3,
    // pitch angle
    phi = from_deg((180 - 108) /2),
    // length of pche, m
    length = 3270e-3,
    // Diameter of semi_circular, m
    d_c = 2e-3,
    // number of channels
    N_ch = integer(125e3),
    // number of segments
    N_seg = 50);
    
  // **** Boundary Conditions as Start values for all components - start ****  
  parameter PCHEBoundaryCondition bc_HTR(
    st_hot_in(p = p_pump_in, T = from_degC(578.22), h = PropsSI("H", "P",  bc_HTR.st_hot_in.p, "T", bc_HTR.st_hot_in.T, PBMedia.mediumName), mdot = mdot_main),    
    st_cold_in(p = p_pump_out, T = T_LTR_cold_out, h = PropsSI("H", "P", bc_HTR.st_cold_in.p, "T", bc_HTR.st_cold_in.T, PBMedia.mediumName), mdot = mdot_main),
    st_hot_out(p = p_pump_in, T = T_HTR_hot_out, h = PropsSI("H", "P", bc_HTR.st_hot_out.p, "T", bc_HTR.st_hot_out.T, PBMedia.mediumName), mdot = mdot_main),
    st_cold_out(p = p_pump_out, T = from_degC(533.5), h = PropsSI("H", "P", bc_HTR.st_cold_out.p, "T", bc_HTR.st_cold_out.T, PBMedia.mediumName), mdot = mdot_main));      

  // boundary condition for LTR test @ diff mdot
  parameter PCHEBoundaryCondition bc_LTR(
    st_hot_in(p = p_pump_in, T = T_HTR_hot_out, h = PropsSI("H", "P", bc_LTR.st_hot_in.p, "T", bc_LTR.st_hot_in.T, PBMedia.mediumName), mdot = mdot_main),    
    st_cold_in(p = p_pump_out, T = from_degC(62.229), h = PropsSI("H", "P", bc_LTR.st_cold_in.p, "T", bc_LTR.st_cold_in.T, PBMedia.mediumName), mdot = mdot_pump),
    st_hot_out(p = p_pump_in, T = from_degC(67.229), h = PropsSI("H", "P", bc_LTR.st_hot_out.p, "T", bc_LTR.st_hot_out.T, PBMedia.mediumName), mdot = mdot_main),
    st_cold_out(p = p_pump_out, T = T_LTR_cold_out, h = PropsSI("H", "P", bc_LTR.st_cold_out.p, "T", bc_LTR.st_cold_out.T, PBMedia.mediumName), mdot = mdot_pump)); 
     
  parameter Steps.Components.ThermoState bc_cooler_out(p = bc_LTR.st_hot_out.p, T = from_degC(33), h = PropsSI("H", "P", bc_cooler_out.p, "T", bc_cooler_out.T, PBMedia.mediumName), mdot = mdot_pump);
  
  parameter Steps.Components.ThermoState bc_heater_out(p = bc_HTR.st_cold_out.p, T = from_degC(700), h = PropsSI("H", "P", bc_heater_out.p, "T", bc_heater_out.T, PBMedia.mediumName), mdot = mdot_main);
  
  // **** Boundary Conditions as Start values for all components - end ****    
    
  /*
  Steps.Components.Regulator regulator_HTR_hot(
    T_init = from_degC(578),
    p_init = P_TURBINE_E,   
    m_flow_init = M_CO2
  );
  
  Steps.Components.Regulator regulator_HTR_cold(
    T_init = T_AMB,
    p_init = P_PUMP_E,   
    m_flow_init = M_CO2
  );
  
  Steps.Components.Regulator regulator_LTR_hot(
    T_init = from_degC(578),
    p_init = P_TURBINE_E,   
    m_flow_init = M_CO2
  );
  
  Steps.Components.Regulator regulator_LTR_cold(
    T_init = T_AMB,
    p_init = P_PUMP_E,   
    m_flow_init = M_CO2
  );*/
  
  Steps.Components.FanCooler fan_cooler(
    T_output = bc_cooler_out.T,      
    inlet.p(start = bc_LTR.st_hot_out.p),
    inlet.h_outflow(start = bc_LTR.st_hot_out.h),
    outlet.p(start = bc_cooler_out.p),
    outlet.h_outflow(start = bc_cooler_out.h)
  );
  
  // use pump as compressor
  Steps.Components.Pump pump(
    p_outlet = bc_LTR.st_cold_in.p, 
    eta = eta_main_compressor,
    outlet.p(start = bc_LTR.st_cold_in.p),
    outlet.h_outflow(start = bc_LTR.st_cold_in.h),
    inlet.p(start = bc_cooler_out.p),
    inlet.h_outflow(start = bc_cooler_out.h)
  );   
  
  Steps.Components.Pump recom_pump(
    p_outlet = bc_HTR.st_cold_in.p, 
    eta = eta_bypass_compressor,
    inlet.p(start = bc_LTR.st_hot_out.p),
    inlet.h_outflow(start = bc_LTR.st_hot_out.h),
    outlet.p(start = bc_HTR.st_cold_in.p),
    outlet.h_outflow(start = bc_HTR.st_cold_in.h)
  );
  
  Components.PCHECImpl HTR(
    name = "HTR", 
    geo = geo_HTR,   
    bc = bc_HTR,   
    // enum log_level {DEBUG = 0, INFO = 1, ERR = 2, SERVE = 3, OFF = 4};
    sim_param(log_level = 4)
    /*
    
    bc_hot_in(p = P_TURBINE_E, T = from_degC(578), h = PropsSI("H", "T", from_degC(578), "P", P_TURBINE_E, PBMedia.mediumName), mdot = M_CO2),    
    bc_cold_in(p = P_PUMP_E, T = from_degC(151), h = PropsSI("H", "T", from_degC(151), "P", P_PUMP_E, PBMedia.mediumName), mdot = M_CO2)*/   
        
  );
  
  Components.PCHECImpl LTR(
    name = "LTR", 
    geo = geo_LTR,   
    bc = bc_LTR,   
    // enum log_level {DEBUG = 0, INFO = 1, ERR = 2, SERVE = 3, OFF = 4};
    sim_param(log_level = 4, step_rel = 0.13)
    /* 
    bc_hot_in(p = P_TURBINE_E, T = from_degC(156), h = PropsSI("H", "T", from_degC(156), "P", P_TURBINE_E, PBMedia.mediumName), mdot = M_CO2),    
    bc_cold_in(p = P_PUMP_E, T = from_degC(62), h = PropsSI("H", "T", from_degC(62), "P", P_PUMP_E, PBMedia.mediumName), mdot = M_CO2)*/
  );

  Steps.Components.TemperatureOutput temp_out(
    T_start = bc_heater_out.T,
    T_stop = bc_heater_out.T,
    dT_step = 50,
    t_sim_duration = 1
  );
  
  Steps.Components.PCMHeater pcm_heater(
    mdot_init = bc_heater_out.mdot,
    T_output_set = bc_heater_out.T,
    /*
    init_outlet = true,
    p_init_outlet = P_PUMP_E,
    T_init_outlet = T_AMB
    */
    
    inlet.p(start = bc_HTR.st_cold_out.p),
    inlet.h_outflow(start = bc_HTR.st_cold_out.h),
    outlet.p(start = bc_heater_out.p),
    outlet.h_outflow(start = bc_heater_out.h)    
    
  );  
   
  Steps.Components.Turbine turbine(    
    p_out = bc_HTR.st_hot_in.p,
    eta = eta_turbine, 
    inlet.p(start = bc_heater_out.p),
    inlet.h_outflow(start = bc_heater_out.h),
    outlet.p(start = bc_HTR.st_hot_in.p),
    outlet.h_outflow(start = bc_HTR.st_hot_in.h)    
  ); 
  
  Steps.Components.Splitter splitter(
    split_ratio = splitter_split_ratio,
    outlet.p(start = bc_LTR.st_hot_out.p),
    outlet.h_outflow(start = bc_LTR.st_hot_out.h),
    inlet.p(start = bc_LTR.st_hot_out.p),
    inlet.h_outflow(start = bc_LTR.st_hot_out.h),
    outlet_split.p(start = bc_LTR.st_hot_out.p),
    outlet_split.h_outflow(start = bc_LTR.st_hot_out.h)
  );
  
  Steps.Components.Merger merger(
    //T_init = T_AMB,
    //p_init = P_ATM,
    outlet.p(start = bc_HTR.st_cold_in.p),
    outlet.h_outflow(start = bc_HTR.st_cold_in.h),
    inlet.p(start = bc_LTR.st_cold_out.p),
    inlet.h_outflow(start = bc_LTR.st_cold_out.h),    
    inlet_merge.p(start = bc_LTR.st_cold_out.p),
    inlet_merge.h_outflow(start = bc_LTR.st_cold_out.h)
  );
  
  //total efficiency
  Real eta_total;


equation
  
  // connect(regulator.outlet, turbine.inlet);
  
  connect(turbine.outlet, HTR.inlet_hot);
  
  // connect(turbine.outlet, regulator_HTR_hot.inlet);  
  // connect(regulator_HTR_hot.outlet, HTR.inlet_hot);
  
  connect(HTR.outlet_hot,LTR.inlet_hot);
  // connect(HTR.outlet_hot,regulator_LTR_hot.inlet);
  // connect(regulator_LTR_hot.outlet, LTR.inlet_hot); 
 
  connect(LTR.outlet_hot, splitter.inlet);
  //recompression loop
  connect(splitter.outlet_split, recom_pump.inlet);
  connect(recom_pump.outlet, merger.inlet_merge);  
  
  connect(splitter.outlet, fan_cooler.inlet);
  
  connect(fan_cooler.outlet, pump.inlet);
  
  connect(pump.outlet, LTR.inlet_cold);  
  
  // connect(pump.outlet, regulator_LTR_cold.inlet);  
  
  // connect(regulator_LTR_cold.outlet, LTR.inlet_cold);  
   
  connect(LTR.outlet_cold, merger.inlet);
  
  connect(merger.outlet, HTR.inlet_cold);
  
  // connect(merger.outlet, regulator_HTR_cold.inlet);
  
  // connect(regulator_HTR_cold.outlet, HTR.inlet_cold);
  
  connect(HTR.outlet_cold, pcm_heater.inlet);    
  
  //connect(temp_out.y, pcm_heater.T_input);
  // connect(pcm_heater.outlet, regulator.inlet);
  connect(pcm_heater.outlet, turbine.inlet);

// algorithm
  eta_total = if initial() then 0 else (turbine.W_turbine - pump.W_comp - recom_pump.W_comp) / pcm_heater.Q * 100;
annotation(
    experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-6, Interval = 1),
    __OpenModelica_commandLineOptions = "--matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection -d=initialization,NLSanalyticJacobian,aliasConflicts");
end OffDesignRCBCycle_v2;
