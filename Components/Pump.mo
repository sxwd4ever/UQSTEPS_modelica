within Steps.Components;

model Pump
  "class for a Pump"
  extends TwoPorts;
  
  import SI = Modelica.SIunits;   
  import Steps.Interfaces.PortType;
  
  
  input Real eta  "efficiency of this pump";   
  input SI.AbsolutePressure p_outlet "fixed outlet pressure of pump"; 
  SI.Power W_comp "power input";
    
  PBMedia.CO2_pT medium_isen;
  
equation
 
  medium_in.state = PBMedia.setState_phX(p = inlet.p, h = inStream(inlet.h_outflow));
  medium_isen.state = PBMedia.setState_psX(p = p_outlet, s = medium_in.s); 
  medium_out.state = PBMedia.setState_phX(p = p_outlet, h = (medium_in.h - (medium_in.h - medium_isen.h) / eta));
  
  inlet.m_flow + outlet.m_flow = 0;
  outlet.p = medium_out.p;
  outlet.T = medium_out.T;
  outlet.h_outflow = medium_out.h; 
  inlet.h_outflow = inStream(outlet.h_outflow); 
  
  W_comp = inlet.m_flow * (medium_out.h - medium_in.h);
  
end Pump;
