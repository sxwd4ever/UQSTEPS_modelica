within Steps.Test;

model TestPCHEMeshram
  "PCHE Test against Meshram [2016]"      
  
  import Modelica.SIunits.Conversions.{from_bar, from_degC}; 
  import Steps.Utilities.CoolProp.PropsSI;
  import Steps.Cycle.OffDPBParamSet;
  import Steps.Components.ThermoState;
  import Steps.Components.PCHEGeoParam;
  import Steps.Components.PCHEBoundaryCondition;
  
  replaceable package PBMedia = Steps.Media.SCO2;
  
  OffDPBParamSet param;
  
  parameter PCHEGeoParam geo = param.geo_HTR;
  parameter PCHEBoundaryCondition bc = param.bc_HTR;
  
  parameter Modelica.SIunits.ReynoldsNumber Re_hot_start = 14.5e3 "Hot stream's start value of Reynolds Number, used to increase convergence";  
  parameter Modelica.SIunits.ReynoldsNumber Re_cold_start = 14.5e3 "Cold stream's start value of Reynolds Number, used to increase convergence";  
  parameter Boolean SourceFixed_hot = true;  
  parameter Boolean SourceFixed_cold = true;
  
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
  
end TestPCHEMeshram;
