within Steps.Test;

model TestPCHEMImpl
  "Modelica implemented PCHE Test against Meshram [2016]"      
  
  import Modelica.SIunits.Conversions.{from_bar, from_degC}; 
  import Steps.Utilities.CoolProp.PropsSI;
  import Steps.Components.{BoundaryCondition, PCHEGeoParam, ThermoState};
  
  replaceable package PBMedia = Steps.Media.SCO2;  

  parameter Modelica.SIunits.MassFlowRate mdot_hot = 51.91 "mass flow rate for hot stream";
  
  parameter Modelica.SIunits.MassFlowRate mdot_cold = 51.91 "mass flow rate for cold stream";
  
  parameter BoundaryCondition bc(
    st_hot_in(p = from_bar(80), T = from_degC(578.22), h = PropsSI("H", "P", bc.st_hot_in.p, "T", bc.st_hot_in.T, PBMedia.mediumName), mdot = mdot_hot),    
    st_cold_in(p = from_bar(200), T = from_degC(151.45), h = PropsSI("H", "P", bc.st_cold_in.p, "T", bc.st_cold_in.T, PBMedia.mediumName), mdot = mdot_hot),
    st_hot_out(p = from_bar(80), T = from_degC(156.5), h = PropsSI("H", "P", bc.st_hot_out.p, "T", bc.st_hot_out.T, PBMedia.mediumName), mdot = mdot_hot),
    st_cold_out(p = from_bar(200), T = from_degC(533.5), h = PropsSI("H", "P", bc.st_cold_out.p, "T", bc.st_cold_out.T, PBMedia.mediumName), mdot = mdot_cold));
  
  parameter PCHEGeoParam geo(
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

  parameter Modelica.SIunits.ReynoldsNumber Re_hot_start = 14.5e3 "Hot stream's start value of Reynolds Number, used to increase convergence";
  
  parameter Modelica.SIunits.ReynoldsNumber Re_cold_start = 14.5e3 "Cold stream's start value of Reynolds Number, used to increase convergence";  
  
  parameter Boolean SourceFixed_hot = true;
  
  parameter Boolean SourceFixed_cold = true;
  
  Components.Source source_hot(
    p_outlet = bc.st_hot_in.p,
    T_outlet = bc.st_hot_in.T,
    mdot_init = mdot_hot,
    fix_state = SourceFixed_hot
  );

  Components.Source source_cold(
    p_outlet = bc.st_cold_in.p,
    T_outlet = bc.st_cold_in.T,
    mdot_init = mdot_cold,
    fix_state = SourceFixed_cold // True if set its state as boundary condition
  );

  Components.Sink sink_hot(
    p_inlet = bc.st_hot_out.p,
    T_inlet = bc.st_hot_out.T,
    mdot_init = mdot_hot,
    fix_state = not SourceFixed_hot
  );

  Components.Sink sink_cold(
    p_inlet = bc.st_cold_out.p,
    T_inlet = bc.st_cold_out.T,
    mdot_init = mdot_cold,
    fix_state = not SourceFixed_cold
  );
 
  Components.PCHeatExchanger pche(
    geo = geo,   
    bc = bc,
    Re_cold_start = Re_cold_start,
    Re_hot_start = Re_hot_start,     
    ByInlet_hot = SourceFixed_hot,
    ByInlet_cold = SourceFixed_cold
  );
  
equation
  
  connect(source_hot.outlet, pche.inlet_hot);
  connect(pche.outlet_hot, sink_hot.inlet);
  connect(source_cold.outlet, pche.inlet_cold);
  connect(pche.outlet_cold, sink_cold.inlet);  

annotation(
    experiment(StartTime = 0, StopTime = 1, Interval = 1, Tolerance = 1e-6),
    __OpenModelica_commandLineOptions = "--matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection -d=initialization,NLSanalyticJacobian,aliasConflicts");
  
end TestPCHEMImpl;
