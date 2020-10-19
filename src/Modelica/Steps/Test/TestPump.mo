within Steps.Test;

model TestPump
  "StandAlone Component Test For Pump"      
  
  import Modelica.SIunits.Conversions.{from_bar, from_degC};
  import Steps.Components.PCHEBoundaryCondition;
  import Steps.Components.ThermoState;
  
  // **** parameter from on-design simulation - START ****
  // UQMECH05_99_CL02_B 1 MW Power Block (RCBC)_2.pdf
  parameter Steps.Cycle.OffDPBParamSet param;    
  
  parameter ThermoState bc_in = param.bc_LTR.st_hot_out; 
  parameter ThermoState bc_out = param.bc_LTR.st_cold_in;
  // **** parameter from on-design simulation - END ****
    
  // **** Arbitary inputs - START **** 
  parameter Modelica.SIunits.AbsolutePressure p_sys = from_bar(90) "p of hot inlet";   
  parameter Modelica.SIunits.Temp_K T_in = 730 "T of hot inlet";      
  parameter Modelica.SIunits.Temp_K T_out = 576.69 "T of hot outlet";    
  parameter Modelica.SIunits.MassFlowRate mdot = 10 "mass flow rate for hot stream";  
  parameter Modelica.SIunits.TemperatureDifference DT_COOLER = 18.0;  
  parameter Modelica.SIunits.Temperature T_AMB = from_degC(15) "Ambinent temperature";  
  parameter Real eta_compressor = param.eta_main_compressor;
  // **** Arbitary inputs - END **** 
  
  Components.Source source(
    p_outlet = bc_in.p,
    T_outlet = bc_in.T,
    mdot_init = bc_in.mdot,
    fix_state = true
  );

  Components.Sink sink(
    p_inlet = bc_out.p,
    T_inlet = bc_out.T,
    mdot_init = bc_out.mdot,
    fix_state = false
  );
  
  // use pump as compressor
  Steps.Components.Pump pump(
    p_outlet = bc_out.p, 
    eta = eta_compressor
  ); 

equation
  
  connect(source.outlet, pump.inlet);  
  connect(pump.outlet, sink.inlet);

end TestPump;
