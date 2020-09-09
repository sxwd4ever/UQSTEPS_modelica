within Steps.Test;

model TestFanCooler
  "StandAlone Component Test For FanCooler"      
  
  import Modelica.SIunits.Conversions.{from_bar, from_degC};    

  parameter Modelica.SIunits.AbsolutePressure p_sys = from_bar(90) "p of hot inlet"; 
  
  parameter Modelica.SIunits.Temp_K T_in = 730 "T of hot inlet";  
    
  parameter Modelica.SIunits.Temp_K T_out = 576.69 "T of hot outlet";  
  
  parameter Modelica.SIunits.MassFlowRate mdot = 10 "mass flow rate for hot stream";
  
  parameter Modelica.SIunits.TemperatureDifference DT_COOLER = 18.0;
  
  parameter Modelica.SIunits.Temperature T_AMB = from_degC(15) "Ambinent temperature";
  
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

  Components.FanCooler cooler(
    T_amb = T_AMB,
    delta_T = DT_COOLER
  );

equation
  
  connect(source.outlet, cooler.inlet);  
  connect(cooler.outlet, sink.inlet);

end TestFanCooler;
