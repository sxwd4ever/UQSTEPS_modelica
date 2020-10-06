within Steps.Components;

model Turbine
  extends Steps.Components.TwoPorts;  
  import CP = Steps.Utilities.CoolProp;
  
  //Initial parameters
  parameter Modelica.SIunits.AbsolutePressure p_out "fixed output pressure";
  parameter Real eta "Efficiency of this turbine";

  //Intermediate variables
  Modelica.SIunits.SpecificEntropy s; 
  
  Modelica.SIunits.SpecificEnthalpy h_isen "Isentropic enthaplpy";  
  
  Modelica.SIunits.Power W_turbine;  
  
equation
  
  s = CP.PropsSI("S", "P", inlet.p, "H", inlet.h_outflow, PBMedia.mediumName);  
  h_isen = CP.PropsSI("H", "P", p_out, "S", s, PBMedia.mediumName); 
  
  //outlet.T = medium_out.T;  
  outlet.m_flow + inlet.m_flow = 0;
  outlet.p = p_out;
  outlet.h_outflow = inlet.h_outflow - eta * (inlet.h_outflow - h_isen);
  
  inlet.h_outflow = inStream(inlet.h_outflow); 
  
  W_turbine = inlet.m_flow * (inlet.h_outflow - outlet.h_outflow);
end Turbine;
