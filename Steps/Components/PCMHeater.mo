within Steps.Components;

model PCMHeater
  "Phase change material Heater"
  extends Steps.Components.TwoPorts;
  
  replaceable package PBMedia = Steps.Media.SCO2;
  
  Modelica.SIunits.TemperatureDifference delta_T "The exit temperature difference between PCM and fluid";

  Modelica.Blocks.Interfaces.RealInput T_input;
  Modelica.SIunits.HeatFlowRate Q;
  
  import Steps.Utilities.CoolProp.PropsSI;

protected
  
  Modelica.SIunits.Temperature T_inlet;

equation
  
  T_inlet = PropsSI("T", "P", inlet.p, "H", inStream(inlet.h_outflow), PBMedia.mediumName);
  
  medium_in.state = PBMedia.setState_phX(p = inlet.p, h = inStream(inlet.h_outflow));
  medium_out.state = PBMedia.setState_pTX(p = inlet.p, T = max(T_input, T_inlet));
  
  outlet.m_flow + inlet.m_flow = 0;
  outlet.h_outflow = medium_out.h;
  //outlet.T = medium_out.T;  
  outlet.p = inlet.p;  
  inlet.h_outflow = inStream(outlet.h_outflow);
  
  Q = inlet.m_flow * (medium_out.h - medium_in.h);

end PCMHeater;
