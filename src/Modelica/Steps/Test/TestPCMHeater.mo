within Steps.Test;

model TestPCMHeater
  "StandAlone Component Test For PCMHeater"      
  
  import Modelica.SIunits.Conversions.{from_bar, from_degC};    

  parameter Modelica.SIunits.AbsolutePressure p_sys = from_bar(90) "p of hot inlet"; 
  
  parameter Modelica.SIunits.Temp_K T_in = 730 "T of hot inlet";  
    
  parameter Modelica.SIunits.Temp_K T_out = 576.69 "T of hot outlet";  
  
  parameter Modelica.SIunits.MassFlowRate mdot = 10 "mass flow rate for hot stream";
  
  parameter Modelica.SIunits.TemperatureDifference DT_COOLER = 18.0;
  
  parameter Modelica.SIunits.Temperature T_AMB = from_degC(15) "Ambinent temperature";  
    
  parameter Modelica.SIunits.Time stop_time = 1.0 "time length of the experiment";
  
  Components.Source source(
    p_outlet = p_sys,
    T_outlet = T_in,
    mdot_init = mdot,
    fix_state = true
  );


  Components.Sink sink(
    p_inlet = p_sys,
    T_inlet = T_out,
    mdot_init = mdot,
    fix_state = false
  );

  Steps.Components.TemperatureOutput temp_out(
    T_start = from_degC(700),
    T_stop = from_degC(700),
    dT_step = 50,
    t_sim_duration = stop_time
  );
  
  // use pump as compressor
  Steps.Components.PCMHeater pcm_heater(
  ); 


equation
  
  connect(temp_out.y, pcm_heater.T_input);
  connect(source.outlet, pcm_heater.inlet);  
  connect(pcm_heater.outlet, sink.inlet);

end TestPCMHeater;
