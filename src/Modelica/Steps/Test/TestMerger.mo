within Steps.Test;

model TestMerger
  "StandAlone Component Test For FanCooler"      
  
  import Modelica.SIunits.Conversions.{from_bar, from_degC};
  import Steps.Components.PCHEBoundaryCondition;
  import Steps.Components.ThermoState;

  parameter Steps.Cycle.OffDPBParamSet param;  
  
  parameter PCHEBoundaryCondition bc_LTR = param.bc_LTR;
  parameter PCHEBoundaryCondition bc_HTR = param.bc_HTR;  
    
  // **** Arbitary inputs ****
  parameter Modelica.SIunits.Pressure P_ATM = 101325; // Pa  
  parameter Modelica.SIunits.Temperature T_AMB = from_degC(15) "Ambinent temperature";
  parameter Modelica.SIunits.AbsolutePressure p_sys = from_bar(90) "p of hot inlet";   
  parameter Modelica.SIunits.Temp_K T_in = 730 "T of hot inlet";      
  parameter Modelica.SIunits.Temp_K T_out = 576.69 "T of hot outlet";    
  parameter Modelica.SIunits.MassFlowRate mdot = 10 "mass flow rate for hot stream";  
  parameter Modelica.SIunits.TemperatureDifference DT_COOLER = 18.0;
  // **** Arbitary inputs - end ****
  
  Components.Source source(
    p_outlet = bc_LTR.st_cold_out.p,
    T_outlet = bc_LTR.st_cold_out.T,
    mdot_init = bc_LTR.st_cold_out.mdot,
    fix_state = true
  );

  Components.Source source_2(
    p_outlet = bc_HTR.st_cold_in.p,
    T_outlet = bc_HTR.st_cold_in.T,
    mdot_init = bc_HTR.st_cold_in.mdot - bc_LTR.st_cold_out.mdot,
    fix_state = true
  );

  Components.Sink sink(
    p_inlet = bc_HTR.st_cold_in.p,
    T_inlet = bc_HTR.st_cold_in.T,
    mdot_init = bc_HTR.st_cold_in.mdot,
    fix_state = false
  );

  Components.Merger merger(
    // T_init = T_AMB,
    // p_init = P_ATM
  );

equation
  
  connect(source.outlet, merger.inlet); 
   
  connect(source_2.outlet, merger.inlet_merge);  
  
  connect(merger.outlet, sink.inlet);

end TestMerger;
