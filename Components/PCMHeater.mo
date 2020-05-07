within Steps.Components;

model PCMHeater
  "Phase change material Heater"
  extends Steps.Components.TwoPorts;
  replaceable package PBMedia = Steps.Media.SCO2;
  
  Modelica.SIunits.TemperatureDifference delta_T "The exit temperature difference between PCM and fluid";
  
  //Real h_i, h_e;
  Real h_e;

  Modelica.Blocks.Interfaces.RealInput T_input;
  PBMedia.ThermodynamicState state_inlet "inlet state";

equation
  
  outlet.m_flow + inlet.m_flow = 0;
  outlet.h_outflow = inStream(inlet.h_outflow);
  outlet.p = inlet.p;
  outlet.T = max(T_input, inlet.T);// - delta_T;
  
  h_e = PBMedia.specificEnthalpy_pTX(outlet.p, outlet.T);

end PCMHeater;
