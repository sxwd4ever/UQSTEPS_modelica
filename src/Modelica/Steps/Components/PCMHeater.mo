within Steps.Components;

model PCMHeater
  "Phase change material Heater"
  extends Steps.Components.TwoPorts;
  
  replaceable package PBMedia = Steps.Media.SCO2;
  
  // Modelica.SIunits.TemperatureDifference delta_T "The exit temperature difference between PCM and fluid";

  parameter Modelica.Blocks.Interfaces.RealInput T_output_set;
  
  Modelica.SIunits.Temperature T_output;
  
  Modelica.SIunits.HeatFlowRate Q;
 
  import Steps.Utilities.CoolProp.PropsSI;
  
  parameter Modelica.SIunits.MassFlowRate mdot_init = 100;

equation
 
  // outlet.m_flow + inlet.m_flow = 0; 
  
  outlet.m_flow = - mdot_init;
  
  T_output = max(T_output_set, T_inlet);
  
  outlet.h_outflow =  PropsSI("H", "P", inlet.p, "T", T_output, PBMedia.mediumName);
  //outlet.T = medium_out.T;  
  outlet.p = inlet.p;  
  inlet.h_outflow = inStream(inlet.h_outflow);
  
  Q = inlet.m_flow * (outlet.h_outflow - inlet.h_outflow);

end PCMHeater;
