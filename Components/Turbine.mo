within Steps.Components;

model Turbine
  extends Steps.Components.TwoPorts;  

  //Initial parameters
  Modelica.SIunits.AbsolutePressure p_out "fixed output pressure";
  Real eta "Efficiency of this turbine";

  //Intermediate variables
  PBMedia.CO2_pT medium_isen "medium under isentropic process";  
  
  Modelica.SIunits.Power W_turbine;  
  
equation
  medium_in.state = PBMedia.setState_pTX(p = inlet.p, T = inlet.T);
  medium_isen.state = PBMedia.setState_psX(p = p_out, s = medium_in.s);
  medium_out.state = PBMedia.setState_phX(p = p_out, h = medium_in.h - eta * (medium_in.h - medium_isen.h));
  
  outlet.T = medium_out.T;  
  outlet.m_flow + inlet.m_flow = 0;
  outlet.p = medium_out.p;
  outlet.h_outflow = medium_out.h;
  inlet.h_outflow = inStream(outlet.h_outflow); 
  
  W_turbine = inlet.m_flow * (medium_in.h - medium_out.h);
end Turbine;
