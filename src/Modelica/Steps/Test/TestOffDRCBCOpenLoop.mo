within Steps.Test;

model TestOffDRCBCOpenLoop
  "Open loop off-Design test for  Brayton Cycle with Recuperator + Recompression"  
// Adjustable parameters
  import Modelica.SIunits.Conversions.{from_degC, from_deg};
  import Modelica.SIunits.{Temperature, Pressure, SpecificEnthalpy};
  import Util = Utilities.Util;
  import Steps.Utilities.CoolProp.PropsSI;  
  import Steps.Components.{PCHEBoundaryCondition, ThermoState, PCHEGeoParam};
  import Steps.Cycle.OffDPBParamSet;   
  
  parameter OffDPBParamSet param; 
  
  parameter PCHEBoundaryCondition bc_LTR = param.bc_LTR;
  parameter PCHEBoundaryCondition bc_HTR = param.bc_HTR;
  parameter ThermoState bc_cooler_out = param.bc_cooler_out;
  parameter ThermoState bc_heater_out = param.bc_heater_out;
  parameter PCHEGeoParam geo_HTR = param.geo_HTR;
  parameter PCHEGeoParam geo_LTR = param.geo_LTR;
    
  Components.Source source(
    p_outlet = bc_HTR.st_hot_in.p,
    T_outlet = bc_HTR.st_hot_in.T,
    mdot_init = bc_HTR.st_hot_in.mdot,
    fix_state = true,
    outlet.p(start = bc_HTR.st_hot_in.p),
    outlet.h_outflow(start = bc_HTR.st_hot_in.h)      
  );

  Components.Sink sink(
    p_inlet = bc_HTR.st_hot_out.p,
    T_inlet = bc_HTR.st_hot_out.T,
    mdot_init = bc_HTR.st_hot_out.mdot,    
    fix_state = false,
    inlet.p(start = bc_HTR.st_hot_out.p),
    inlet.h_outflow(start = bc_HTR.st_hot_out.h)    
  );   
  
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
    eta = param.eta_main_compressor,
    outlet.p(start = bc_LTR.st_cold_in.p),
    outlet.h_outflow(start = bc_LTR.st_cold_in.h),
    inlet.p(start = bc_cooler_out.p),
    inlet.h_outflow(start = bc_cooler_out.h)
  );   
  
  Steps.Components.Pump recom_pump(
    p_outlet = bc_HTR.st_cold_in.p, 
    eta = param.eta_bypass_compressor,
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

  Steps.Components.Turbine turbine(    
    p_out = bc_HTR.st_hot_in.p,
    eta = param.eta_turbine, 
    inlet.p(start = bc_heater_out.p),
    inlet.h_outflow(start = bc_heater_out.h),
    outlet.p(start = bc_HTR.st_hot_in.p),
    outlet.h_outflow(start = bc_HTR.st_hot_in.h)    
  ); 
  
  Steps.Components.Splitter splitter(
    split_ratio = param.splitter_split_ratio,
    outlet.p(start = bc_LTR.st_hot_out.p),
    outlet.h_outflow(start = bc_LTR.st_hot_out.h),
    inlet.p(start = bc_LTR.st_hot_out.p),
    inlet.h_outflow(start = bc_LTR.st_hot_out.h),
    outlet_split.p(start = bc_LTR.st_hot_out.p),
    outlet_split.h_outflow(start = bc_LTR.st_hot_out.h)
  );
  
  Steps.Components.Mixer mixer(
    //T_init = T_AMB,
    //p_init = P_ATM,
    outlet.p(start = bc_HTR.st_cold_in.p),
    outlet.h_outflow(start = bc_HTR.st_cold_in.h),
    inlet.p(start = bc_LTR.st_cold_out.p),
    inlet.h_outflow(start = bc_LTR.st_cold_out.h),    
    inlet_mix.p(start = bc_LTR.st_cold_out.p),
    inlet_mix.h_outflow(start = bc_LTR.st_cold_out.h)
  );
  
  //total efficiency
  Real eta_total;


equation
  
  // connect(regulator.outlet, turbine.inlet);
  //connect(temp_out.y, pcm_heater.T_input);
  // connect(pcm_heater.outlet, regulator.inlet);
  connect(source.outlet, turbine.inlet);
    
  connect(turbine.outlet, HTR.inlet_hot);
  
  // connect(turbine.outlet, regulator_HTR_hot.inlet);  
  // connect(regulator_HTR_hot.outlet, HTR.inlet_hot);
  
  connect(HTR.outlet_hot,LTR.inlet_hot);
  // connect(HTR.outlet_hot,regulator_LTR_hot.inlet);
  // connect(regulator_LTR_hot.outlet, LTR.inlet_hot); 
 
  connect(LTR.outlet_hot, splitter.inlet);
  //recompression loop
  connect(splitter.outlet_split, recom_pump.inlet);
  connect(recom_pump.outlet, mixer.inlet_mix);  
  
  connect(splitter.outlet, fan_cooler.inlet);
  
  connect(fan_cooler.outlet, pump.inlet);
  
  connect(pump.outlet, LTR.inlet_cold);  
  
  // connect(pump.outlet, regulator_LTR_cold.inlet);  
  
  // connect(regulator_LTR_cold.outlet, LTR.inlet_cold);  
   
  connect(LTR.outlet_cold, mixer.inlet);
  
  connect(mixer.outlet, HTR.inlet_cold);
  
  // connect(mixer.outlet, regulator_HTR_cold.inlet);
  
  // connect(regulator_HTR_cold.outlet, HTR.inlet_cold);
  
  connect(HTR.outlet_cold, sink.inlet);     


// algorithm
  eta_total = 1; //if initial() then 0 else (turbine.W_turbine - pump.W_comp - recom_pump.W_comp) / pcm_heater.Q * 100;
annotation(
    experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-2, Interval = 1),
    __OpenModelica_commandLineOptions = "--matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection -d=initialization,NLSanalyticJacobian,aliasConflicts");
end TestOffDRCBCOpenLoop;
