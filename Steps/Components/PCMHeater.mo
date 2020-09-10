within Steps.Components;

model PCMHeater
  "Phase change material Heater"
  extends Steps.Components.TwoPorts;
  
  replaceable package PBMedia = Steps.Media.SCO2;
  
  // Modelica.SIunits.TemperatureDifference delta_T "The exit temperature difference between PCM and fluid";

  Modelica.Blocks.Interfaces.RealInput T_input;
  
  Modelica.SIunits.Temperature T_output;
  
  Modelica.SIunits.HeatFlowRate Q;
 
  import Steps.Utilities.CoolProp.PropsSI;

protected
  
  Modelica.SIunits.Temperature T_inlet;

equation
  
  T_inlet = PropsSI("T", "P", inlet.p, "H", inlet.h_outflow, PBMedia.mediumName);
 
  outlet.m_flow + inlet.m_flow = 0;
  T_output = max(T_input, T_inlet);
  
  outlet.h_outflow =  PropsSI("H", "P", inlet.p, "T", T_output, PBMedia.mediumName);
  //outlet.T = medium_out.T;  
  outlet.p = inlet.p;  
  inlet.h_outflow = inStream(inlet.h_outflow);
  
  Q = inlet.m_flow * (outlet.h_outflow - inlet.h_outflow);

end PCMHeater;
