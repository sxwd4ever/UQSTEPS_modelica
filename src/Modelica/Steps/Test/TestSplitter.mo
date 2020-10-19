within Steps.Test;

model TestSplitter
  "StandAlone Component Test For FanCooler"      
  
  import Modelica.SIunits.Conversions.{from_bar, from_degC}; 
  import Steps.Components.PCHEBoundaryCondition;
  import Steps.Components.ThermoState;
  
  // **** parameter from on-design simulation - START ****
  // UQMECH05_99_CL02_B 1 MW Power Block (RCBC)_2.pdf
  parameter Steps.Cycle.OffDPBParamSet param;    
  
  parameter ThermoState bc_in = param.bc_LTR.st_hot_out; 
  parameter ThermoState bc_out = param.bc_LTR.st_cold_out;
  parameter ThermoState bc_out_split = param.bc_LTR.st_hot_out;  
  // **** parameter from on-design simulation - END ****
    
  // **** Arbitary inputs - START ****      
  parameter Modelica.SIunits.Pressure P_ATM = 101325; // Pa  
  parameter Modelica.SIunits.Temperature T_AMB = from_degC(15) "Ambinent temperature";
  parameter Modelica.SIunits.AbsolutePressure p_sys = from_bar(90) "p of hot inlet";  
  parameter Modelica.SIunits.Temp_K T_in = 730 "T of hot inlet";     
  parameter Modelica.SIunits.Temp_K T_out = 576.69 "T of hot outlet";  
  parameter Modelica.SIunits.MassFlowRate mdot = 10 "mass flow rate for hot stream"; 
  parameter Modelica.SIunits.TemperatureDifference DT_COOLER = 18.0;    
  // **** Arbitary inputs - END ****      
  
  Components.Source source(
    p_outlet = bc_in.p,
    T_outlet = bc_in.T,
    mdot_init = param.mdot_main,
    fix_state = true
  );

  Components.Sink sink(
    p_inlet = bc_out_split.p,
    T_inlet = bc_out_split.T,
    mdot_init = param.mdot_pump,
    fix_state = false
  );

  Components.Sink sink_2(
    p_inlet = bc_out.p,
    T_inlet = bc_out.T,
    mdot_init = param.mdot_bypass,
    fix_state = false
  );

  Components.Splitter splitter(
    split_ratio = 0.4
  );

equation
  
  connect(source.outlet, splitter.inlet);    
  connect(splitter.outlet, sink.inlet);   
  connect(splitter.outlet_split, sink_2.inlet);  

end TestSplitter;
