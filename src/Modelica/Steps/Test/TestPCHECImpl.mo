within Steps.Test;

model TestPCHECImpl
  "PCHE c-based implementation"      
  
  import Modelica.SIunits.Conversions.{from_bar, from_deg, from_degC};   
  import Steps.Utilities.CoolProp.PropsSI;
  import Steps.Components.{PCHEBoundaryCondition, PCHEGeoParam, ThermoState};
  
  replaceable package PBMedia = Steps.Media.SCO2; 
  
  parameter Modelica.SIunits.Pressure p_pump_in = 8e6;
  parameter Modelica.SIunits.Pressure p_pump_out = 20e6;
  
  parameter Modelica.SIunits.MassFlowRate mdot_main = 51.91;  
  parameter Modelica.SIunits.MassFlowRate mdot_pump = 31.31;  
  
  parameter Modelica.SIunits.Temperature T_HTR_hot_out = from_degC(156.45);  
  parameter Modelica.SIunits.Temperature T_LTR_cold_out = from_degC(151.45);  
   
  parameter Boolean SourceFixed_hot = true;  
  parameter Boolean SourceFixed_cold = true;
  
  // **** candinate parameter set - start ****
  parameter PCHEBoundaryCondition bc_HTR(
    st_hot_in(p = p_pump_in, T = from_degC(578.22), h = PropsSI("H", "P", bc_HTR.st_hot_in.p, "T", bc_HTR.st_hot_in.T, PBMedia.mediumName), mdot = mdot_main),    
    st_cold_in(p = p_pump_out, T = T_LTR_cold_out, h = PropsSI("H", "P", bc_HTR.st_cold_in.p, "T", bc_HTR.st_cold_in.T, PBMedia.mediumName), mdot = mdot_main),
    st_hot_out(p = p_pump_in, T = T_HTR_hot_out, h = PropsSI("H", "P", bc_HTR.st_hot_out.p, "T", bc_HTR.st_hot_out.T, PBMedia.mediumName), mdot = mdot_main),
    st_cold_out(p = p_pump_out, T = from_degC(533.5), h = PropsSI("H", "P", bc_HTR.st_cold_out.p, "T", bc_HTR.st_cold_out.T, PBMedia.mediumName), mdot = mdot_main));
  
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

  // boundary condition for LTR test @ diff mdot
  parameter PCHEBoundaryCondition bc_LTR(
    st_hot_in(p = p_pump_in, T = T_HTR_hot_out, h = PropsSI("H", "P", bc_LTR.st_hot_in.p, "T", bc_LTR.st_hot_in.T, PBMedia.mediumName), mdot = mdot_main),    
    st_cold_in(p = p_pump_out, T = from_degC(62.229), h = PropsSI("H", "P", bc_LTR.st_cold_in.p, "T", bc_LTR.st_cold_in.T, PBMedia.mediumName), mdot = mdot_pump),
    st_hot_out(p = p_pump_in, T = from_degC(67.229), h = PropsSI("H", "P", bc_LTR.st_hot_out.p, "T", bc_LTR.st_hot_out.T, PBMedia.mediumName), mdot = mdot_main),
    st_cold_out(p = p_pump_out, T = T_LTR_cold_out, h = PropsSI("H", "P", bc_LTR.st_cold_out.p, "T", bc_LTR.st_cold_out.T, PBMedia.mediumName), mdot = mdot_pump));
    
  parameter PCHEGeoParam geo_LTR(
    // pitch length, m
    pitch = 12e-3,
    // pitch angle
    phi = from_deg((180 - 108) /2),
    // length of pche, m
    length = 3270e-3, // 1100e-3,
    // Diameter of semi_circular, m
    d_c = 2e-3,
    // number of channels
    N_ch = integer(125e3),
    // number of segments
    N_seg = 50);
  // **** candinate parameter set - end ****  
  
  // select the group of parameter set as input here
  parameter PCHEBoundaryCondition bc = bc_LTR;
  parameter PCHEGeoParam geo = geo_LTR;  

  Components.Source source_hot(
    p_outlet = bc.st_hot_in.p,
    T_outlet = bc.st_hot_in.T,
    mdot_init = bc.st_hot_in.mdot,
    fix_state = SourceFixed_hot
  );

  Components.Source source_cold(
    p_outlet = bc.st_cold_in.p,
    T_outlet = bc.st_cold_in.T,
    mdot_init = bc.st_cold_in.mdot,
    fix_state = SourceFixed_cold // True if set its state as boundary condition
  );

  Components.Sink sink_hot(
    p_inlet = bc.st_hot_out.p,
    T_inlet = bc.st_hot_out.T,
    mdot_init = bc.st_hot_out.mdot,
    fix_state = not SourceFixed_hot
  );

  Components.Sink sink_cold(
    p_inlet = bc.st_cold_out.p,
    T_inlet = bc.st_cold_out.T,
    mdot_init = bc.st_cold_out.mdot,
    fix_state = not SourceFixed_cold
  );

  Components.PCHECImpl pche(
    geo = geo,   
    bc = bc,   
    sim_param(log_level = 1, step_rel = 0.13) // step_rel will affect result's error and simulation speed
  );

equation  
  
  connect(source_hot.outlet, pche.inlet_hot);
  connect(pche.outlet_hot, sink_hot.inlet);
  connect(source_cold.outlet, pche.inlet_cold);
  connect(pche.outlet_cold, sink_cold.inlet); 
   
annotation(
  experiment(StartTime = 0, StopTime = 1, Interval = 1, Tolerance = 1e-6),
    __OpenModelica_commandLineOptions = "--matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection -d=initialization,NLSanalyticJacobian,aliasConflicts");  
    
end TestPCHECImpl;
