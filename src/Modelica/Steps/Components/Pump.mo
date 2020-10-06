within Steps.Components;

model Pump
  "class for a Pump"
  extends TwoPorts;
  
  import SI = Modelica.SIunits;   
  import CP = Steps.Utilities.CoolProp;
  
  import Steps.Interfaces.PortType;  
  
  input Real eta  "efficiency of this pump";   
  input SI.AbsolutePressure p_outlet "fixed outlet pressure of pump"; 
  SI.Power W_comp "power input";   

  
  Modelica.SIunits.Entropy s "Isentropic entropy";
  
  Modelica.SIunits.SpecificEnthalpy h_isen "Isentropic enthaplpy";
  
equation
 
  s = CP.PropsSI("S", "P", inlet.p, "H", inlet.h_outflow, PBMedia.mediumName);  
  h_isen = CP.PropsSI("H", "P", p_outlet, "S", s, PBMedia.mediumName);
    
  inlet.m_flow + outlet.m_flow = 0;
  outlet.p = p_outlet;
  
  outlet.h_outflow = (inlet.h_outflow - (inlet.h_outflow - h_isen) / eta);
  
  inlet.h_outflow = inStream(inlet.h_outflow);
  
  W_comp = inlet.m_flow * (outlet.h_outflow - inlet.h_outflow);
  
end Pump;
