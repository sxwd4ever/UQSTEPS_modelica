within Steps.Test;

model TestComponentSeries
  "Test a series of Components for  Brayton Cycle with Recuperator + Recompression"  
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
  parameter ThermoState bc_bypass = param.bc_bypass;
  parameter PCHEGeoParam geo_HTR = param.geo_HTR;
  parameter PCHEGeoParam geo_LTR = param.geo_LTR;
    
  Components.Source hot_source(
    p_outlet = bc_HTR.st_hot_in.p,
    T_outlet = bc_HTR.st_hot_in.T,
    mdot_init = bc_HTR.st_hot_in.mdot,
    fix_state = true/*,
    outlet.p(start = bc_HTR.st_hot_in.p),
    outlet.h_outflow(start = bc_HTR.st_hot_in.h)*/      
  );

  Components.Sink cold_sink(
    p_inlet = bc_HTR.st_cold_out.p,
    T_inlet = bc_HTR.st_cold_out.T,
    mdot_init = bc_HTR.st_cold_out.mdot,    
    fix_state = false/*,
    inlet.p(start = bc_HTR.st_cold_out.p),
    inlet.h_outflow(start = bc_HTR.st_cold_out.h)*/    
  ); 

  Components.Source cold_source(
    p_outlet = bc_LTR.st_cold_out.p,
    T_outlet = bc_LTR.st_cold_out.T,
    mdot_init = bc_LTR.st_cold_out.mdot,
    fix_state = true/*,
    outlet.p(start = bc_LTR.st_cold_in.p),
    outlet.h_outflow(start = bc_LTR.st_cold_in.h)*/     
  );

  Components.Sink hot_sink(
    p_inlet = bc_LTR.st_hot_in.p,
    T_inlet = bc_LTR.st_hot_in.T,
    mdot_init = bc_LTR.st_hot_in.mdot,    
    fix_state = false/*,
    inlet.p(start = bc_LTR.st_hot_out.p),
    inlet.h_outflow(start = bc_LTR.st_hot_out.h)*/    
  );  

  Components.Source merger_source(
    p_outlet = bc_bypass.p,
    T_outlet = bc_bypass.T,
    mdot_init = bc_bypass.mdot,
    fix_state = true/*,
    outlet.p(start = bc_bypass.p),
    outlet.h_outflow(start = bc_bypass.h)*/      
  );

  Components.PCHECImpl HTR(
    name = "HTR", 
    geo = geo_HTR,   
    bc = bc_HTR,   
    // enum log_level {DEBUG = 0, INFO = 1, ERR = 2, SERVE = 3, OFF = 4};
    sim_param(err = 5e-3, log_level = 1, step_rel = 0.13)
    /*    
    bc_hot_in(p = P_TURBINE_E, T = from_degC(578), h = PropsSI("H", "T", from_degC(578), "P", P_TURBINE_E, PBMedia.mediumName), mdot = M_CO2),    
    bc_cold_in(p = P_PUMP_E, T = from_degC(151), h = PropsSI("H", "T", from_degC(151), "P", P_PUMP_E, PBMedia.mediumName), mdot = M_CO2)*/           
  );
  /* 
  Components.PCHECImpl LTR(
    name = "LTR", 
    geo = geo_LTR,   
    bc = bc_LTR,   
    // enum log_level {DEBUG = 0, INFO = 1, ERR = 2, SERVE = 3, OFF = 4};
    sim_param(log_level = 1, step_rel = 0.13)
    
    bc_hot_in(p = P_TURBINE_E, T = from_degC(156), h = PropsSI("H", "T", from_degC(156), "P", P_TURBINE_E, PBMedia.mediumName), mdot = M_CO2),    
    bc_cold_in(p = P_PUMP_E, T = from_degC(62), h = PropsSI("H", "T", from_degC(62), "P", P_PUMP_E, PBMedia.mediumName), mdot = M_CO2)
  ); */ 

  Steps.Components.Merger merger(
      //T_init = T_AMB,
      //p_init = P_ATM,
      /*
      outlet.p(start = bc_HTR.st_cold_in.p),
      outlet.h_outflow(start = bc_HTR.st_cold_in.h),
      inlet.p(start = bc_LTR.st_cold_out.p),
      inlet.h_outflow(start = bc_LTR.st_cold_out.h),    
      inlet_merge.p(start = bc_bypass.p),
      inlet_merge.h_outflow(start = bc_bypass.h)*/
    );
equation
  /*
  connect(hot_source.outlet, HTR.inlet_hot);
  connect(HTR.outlet_hot, LTR.inlet_hot); 
  connect(LTR.outlet_hot, hot_sink.inlet);

  //recompression loop  
  connect(merger_source.outlet, merger.inlet_merge);  
  connect(cold_source.outlet, LTR.inlet_cold); 
  connect(LTR.outlet_cold, merger.inlet);
  connect(merger.outlet, HTR.inlet_cold);
  connect(HTR.outlet_cold, cold_sink.inlet);*/  
  
  connect(hot_source.outlet, HTR.inlet_hot);
  connect(HTR.outlet_hot, hot_sink.inlet);

  //recompression loop  
  connect(merger_source.outlet, merger.inlet_merge);  
  connect(cold_source.outlet, merger.inlet);
  connect(merger.outlet, HTR.inlet_cold);
  connect(HTR.outlet_cold, cold_sink.inlet);   

annotation(
    experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-2, Interval = 1),
    __OpenModelica_commandLineOptions = "--matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection -d=initialization,NLSanalyticJacobian,aliasConflicts");
end TestComponentSeries;
