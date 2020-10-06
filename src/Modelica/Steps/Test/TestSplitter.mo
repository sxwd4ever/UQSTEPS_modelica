within Steps.Test;

model TestSplitter
  "StandAlone Component Test For FanCooler"      
  
  import Modelica.SIunits.Conversions.{from_bar, from_degC};    

  parameter Modelica.SIunits.Pressure P_ATM = 101325; // Pa
  
  parameter Modelica.SIunits.Temperature T_AMB = from_degC(15) "Ambinent temperature";

  parameter Modelica.SIunits.AbsolutePressure p_sys = from_bar(90) "p of hot inlet"; 
  
  parameter Modelica.SIunits.Temp_K T_in = 730 "T of hot inlet";  
    
  parameter Modelica.SIunits.Temp_K T_out = 576.69 "T of hot outlet";  
  
  parameter Modelica.SIunits.MassFlowRate mdot = 10 "mass flow rate for hot stream";
  
  parameter Modelica.SIunits.TemperatureDifference DT_COOLER = 18.0;
  
  parameter Modelica.SIunits.Temperature T_AMB = from_degC(15) "Ambinent temperature";
  
  Components.Source source(
    p_outlet = p_sys,
    T_outlet = T_in + 10,
    mdot_init = mdot / 2,
    fix_state = true
  );

  Components.Sink sink(
    p_inlet = p_sys,
    T_inlet = T_out,
    mdot_init = mdot,
    fix_state = false
  );

  Components.Sink sink_2(
    p_inlet = p_sys,
    T_inlet = T_out,
    mdot_init = mdot,
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
