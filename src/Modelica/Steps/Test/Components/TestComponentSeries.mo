within Steps.Test.Components;

model TestComponentSeries
  "Test a series of Components for  Brayton Cycle with Recuperator + Recompression"  
// Adjustable parameters
  import Modelica.SIunits.Conversions.{from_degC, from_deg};
  import Modelica.SIunits.{Temperature, Pressure, SpecificEnthalpy};
  import Util = Utilities.Util;
  import Steps.Utilities.CoolProp.PropsSI;  
  import Steps.Components.{PCHEGeoParam, PCHEBoundaryCondition};
  import Steps.Model.{PBConfiguration, SimParam, EntityGeoParam, ThermoState} ;
  
  // parameter OffDPBParamSet param(bc_HTR.st_cold_in.T = from_degC(151.993)); 
  // parameter for modelica implementation
  parameter PBConfiguration param_MImpl(geo_HTR.N_seg = 50, geo_LTR.N_seg = 50); 
  
  // parameter for C++ implementation of PCHE - based on Modelica impl's result
  parameter Model.PBConfiguration param_CImpl(
  bc_HTR.st_cold_out.T = from_degC(495.302),
  bc_HTR.st_cold_in.T = from_degC(141.3),
  bc_HTR.st_hot_out.T = from_degC(141.041),
  bc_LTR.st_cold_out.T = from_degC(130.631),
  bc_LTR.st_hot_in.T = from_degC(141.007),
  bc_LTR.st_hot_out.T = from_degC(63.6726)); 
  
  parameter PBConfiguration param_def; 
  
  // parameter I get after a success simulation with a CImpl model
  parameter PBConfiguration param_CImpl_v2(
  bc_HTR.st_cold_out.T = from_degC(532.62),
  bc_HTR.st_cold_in.T = from_degC(143.574),
  bc_HTR.st_hot_out.T = from_degC(144.604),
  bc_HTR.st_hot_out.p = from_bar(77.3621),
  bc_LTR.st_cold_out.T = from_degC(138.571),
  bc_LTR.st_hot_in.T = from_degC(144.604),
  bc_LTR.st_hot_in.p = from_bar(77.3621),
  bc_LTR.st_hot_out.T = from_degC(65.0493),
  bc_LTR.st_hot_out.p = from_bar(76.3718)); 
  
  parameter PBConfiguration param = param_CImpl_v2; 
  
  parameter PCHEBoundaryCondition bc_LTR = param.bc_LTR;
  parameter PCHEBoundaryCondition bc_HTR = param.bc_HTR;
  parameter ThermoState bc_bypass = param.bc_bypass;
  parameter PCHEGeoParam geo_HTR = param.geo_HTR;
  parameter PCHEGeoParam geo_LTR = param.geo_LTR;
  parameter SimParam sim_param = param.sim_param_def;
    
  parameter Modelica.SIunits.ReynoldsNumber Re_hot_start = 14.5e3 "Hot stream's start value of Reynolds Number, used to increase convergence";  
  parameter Modelica.SIunits.ReynoldsNumber Re_cold_start = 14.5e3 "Cold stream's start value of Reynolds Number, used to increase convergence";  
  parameter Boolean SourceFixed_hot = true;  
  parameter Boolean SourceFixed_cold = true;    
    
  Components.Source hot_source(
    p_outlet = bc_HTR.st_hot_in.p,
    T_outlet = bc_HTR.st_hot_in.T,
    mdot_init = bc_HTR.st_hot_in.mdot,
    fix_state = true,
    outlet.p(start = bc_HTR.st_hot_in.p),
    outlet.h_outflow(start = bc_HTR.st_hot_in.h)     
  );

  Components.Sink cold_sink(
    p_inlet = bc_HTR.st_cold_out.p,
    T_inlet = bc_HTR.st_cold_out.T,
    mdot_init = bc_HTR.st_cold_out.mdot,    
    fix_state = false,
    inlet.p(start = bc_HTR.st_cold_out.p),
    inlet.h_outflow(start = bc_HTR.st_cold_out.h)  
  ); 

  Components.Source cold_source(
    p_outlet = bc_LTR.st_cold_in.p,
    T_outlet = bc_LTR.st_cold_in.T,
    mdot_init = bc_LTR.st_cold_in.mdot,
    fix_state = true,
    outlet.p(start = bc_LTR.st_cold_in.p),
    outlet.h_outflow(start = bc_LTR.st_cold_in.h)     
  );

  Components.Sink hot_sink(
    p_inlet = bc_LTR.st_hot_out.p,
    T_inlet = bc_LTR.st_hot_out.T,
    mdot_init = bc_LTR.st_hot_out.mdot,    
    fix_state = false,
    inlet.p(start = bc_LTR.st_hot_out.p),
    inlet.h_outflow(start = bc_LTR.st_hot_out.h)    
  );  

  Components.Source mixer_source(
    p_outlet = bc_bypass.p,
    T_outlet = bc_bypass.T,
    mdot_init = bc_bypass.mdot,
    fix_state = true 
  );  

  Components.PCHECImpl HTR(
    name = "HTR", 
    geo = geo_HTR,   
    bc = bc_HTR,   
    init_cold_in = true,
    init_hot_out = false,
    // enum log_level {DEBUG = 0, INFO = 1, ERR = 2, SERVE = 3, OFF = 4};
    sim_param = sim_param        
  );

  /* 
  Components.PCHeatExchanger HTR(
    geo = geo_HTR,
    bc = bc_HTR,
    Re_cold_start = Re_cold_start,
    Re_hot_start = Re_hot_start,     
    ByInlet_hot = SourceFixed_hot,
    ByInlet_cold = SourceFixed_cold 
  );
  */
  
  Components.PCHECImpl LTR(
    name = "LTR", 
    geo = geo_LTR,   
    bc = bc_LTR, 
    init_hot_in = false,
    // enum log_level {DEBUG = 0, INFO = 1, ERR = 2, SERVE = 3, OFF = 4};
    sim_param = sim_param
  );
   
  /*
  Components.PCHeatExchanger LTR(
    geo = geo_LTR,
    bc = bc_LTR,
    Re_cold_start = Re_cold_start,
    Re_hot_start = Re_hot_start,     
    ByInlet_hot = SourceFixed_hot,
    ByInlet_cold = SourceFixed_cold 
  );
  */
  Steps.Components.Mixer mixer(
      //T_init = T_AMB,
      //p_init = P_ATM,
      /*
      outlet.p(start = bc_HTR.st_cold_in.p),
      outlet.h_outflow(start = bc_HTR.st_cold_in.h),
      inlet.p(start = bc_LTR.st_cold_out.p),
      inlet.h_outflow(start = bc_LTR.st_cold_out.h),    
      inlet_mix.p(start = bc_bypass.p),
      inlet_mix.h_outflow(start = bc_bypass.h)*/
  );
  
   
equation
  
  // HTR + LTR + Mixer
  connect(hot_source.outlet, HTR.inlet_hot);
  connect(HTR.outlet_hot, LTR.inlet_hot); 
  connect(LTR.outlet_hot, hot_sink.inlet);

  //recompression loop  
  connect(mixer_source.outlet, mixer.inlet_mix);  
  connect(cold_source.outlet, LTR.inlet_cold); 
  connect(LTR.outlet_cold, mixer.inlet);
  connect(mixer.outlet, HTR.inlet_cold);
  connect(HTR.outlet_cold, cold_sink.inlet);  
  
  /*
  // LTR + Mixer   
  connect(hot_source.outlet, LTR.inlet_hot);
  connect(LTR.outlet_hot, hot_sink.inlet);

  //recompression loop  
  //connect(mixer_source.outlet, mixer.inlet_mix);  
  //recompression loop  
  connect(mixer_source.outlet, mixer.inlet_mix);  
  connect(cold_source.outlet, LTR.inlet_cold); 
  connect(LTR.outlet_cold, mixer.inlet);
  connect(mixer.outlet, cold_sink.inlet);    
  */
  /*
  // HTR + Mixer 
  connect(hot_source.outlet, HTR.inlet_hot);
  connect(HTR.outlet_hot, hot_sink.inlet);

  //recompression loop  
  connect(mixer_source.outlet, mixer.inlet_mix);  
  connect(cold_source.outlet, mixer.inlet);
  connect(mixer.outlet, HTR.inlet_cold);
  connect(HTR.outlet_cold, cold_sink.inlet);   
  */
annotation(
    experiment(StartTime = 0, StopTime = 1, Tolerance = 1, Interval = 1),
    __OpenModelica_commandLineOptions = "--matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection -d=initialization,NLSanalyticJacobian,aliasConflicts,bltdump");/*, 
    __OpenModelica_simulationFlags(lv = "LOG_DEBUG,LOG_EVENTS_V,LOG_INIT,LOG_INIT_V,LOG_LS,LOG_LS_V,LOG_NLS,LOG_NLS_V,LOG_NLS_JAC,LOG_NLS_RES,LOG_RES_INIT,LOG_STATS,LOG_SOTI,LOG_SIMULATION", outputFormat = "mat", s = "dassl", nls = "homotopy"));*/
end TestComponentSeries;
