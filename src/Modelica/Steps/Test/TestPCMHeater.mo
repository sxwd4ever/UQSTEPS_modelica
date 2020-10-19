within Steps.Test;

model TestPCMHeater
  "StandAlone Component Test For PCMHeater"      
  
  import Modelica.SIunits.Conversions.{from_bar, from_degC};  
  import Steps.Components.PCHEBoundaryCondition;
  import Steps.Components.ThermoState;
  
  // **** parameter from on-design simulation - START ****
  // UQMECH05_99_CL02_B 1 MW Power Block (RCBC)_2.pdf
  parameter Steps.Cycle.OffDPBParamSet param;    
  
  parameter ThermoState bc_heater_in = param.bc_HTR.st_cold_out;  
  parameter ThermoState bc_heater_out = param.bc_heater_out;
  // **** parameter from on-design simulation - END ****
    
  // **** Arbitary inputs ****   
  parameter Modelica.SIunits.AbsolutePressure p_sys = from_bar(90) "p of hot inlet";   
  parameter Modelica.SIunits.Temp_K T_in = 730 "T of hot inlet";      
  parameter Modelica.SIunits.Temp_K T_out = 576.69 "T of hot outlet";    
  parameter Modelica.SIunits.MassFlowRate mdot = 10 "mass flow rate for hot stream";  
  parameter Modelica.SIunits.TemperatureDifference DT_COOLER = 18.0;  
  parameter Modelica.SIunits.Temperature T_AMB = from_degC(15) "Ambinent temperature";      
  parameter Modelica.SIunits.Time stop_time = 1.0 "time length of the experiment";
  // **** Arbitary inputs - end ****   
  
  Components.Source source(
    p_outlet = bc_heater_in.p,
    T_outlet = bc_heater_in.T,
    mdot_init = bc_heater_in.mdot,
    fix_state = true
  );

  Components.Sink sink(
    p_inlet = bc_heater_out.p,
    T_inlet = bc_heater_out.T,
    mdot_init = bc_heater_out.mdot,
    fix_state = false
  );
  /*
  Steps.Components.TemperatureOutput temp_out(
    T_start = bc_heater_out.T - 10,
    T_stop = bc_heater_out.T + 10,
    dT_step = 50,
    t_sim_duration = stop_time
  );
  */
  // use pump as compressor
  Steps.Components.PCMHeater pcm_heater(
    T_output_set = bc_heater_out.T,
    SetMdot = false
  ); 


equation
  
  //connect(temp_out.y, pcm_heater.T_output_set);
  connect(source.outlet, pcm_heater.inlet);  
  connect(pcm_heater.outlet, sink.inlet);

end TestPCMHeater;
