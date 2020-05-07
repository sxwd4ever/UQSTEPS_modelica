within Steps.Components;

model PCMHeater
  "Phase change material Heater"
  extends Steps.Components.TwoPorts;
  
  replaceable package PBMedia = Steps.Media.SCO2;
  
  Modelica.SIunits.TemperatureDifference delta_T "The exit temperature difference between PCM and fluid";
  
  //Real h_i, h_e;
  //Real h_e;

  Modelica.Blocks.Interfaces.RealInput T_input;
  //PBMedia.ThermodynamicState state_inlet "inlet state";

equation
  
  medium_in.state = SCO2.setState_pTX(p = inlet.p, T = inlet.T);
  medium_out.state = SCO2.setState_pTX(p = inlet.p, T = max(T_input, inlet.T));
  
  outlet.m_flow + inlet.m_flow = 0;
  outlet.h_outflow = medium_out.h;
  //outlet.h_outflow = inStream(inlet.h_outflow);
  outlet.p = inlet.p;
  outlet.T = medium_out.T;// - delta_T;
  
  inlet.h_outflow = inStream(outlet.h_outflow);
  //h_e = PBMedia.specificEnthalpy_pTX(outlet.p, outlet.T);

end PCMHeater;
