within Steps.Components;

model Turbine
  extends Steps.Components.TwoPorts;
  
  replaceable package PBMedia = Steps.Media.SCO2;
  
  //Modelica.SIunits.Enthalpy h_i;
  //Modelica.SIunits.Enthalpy h_es "idealized outlet enthalpy";
  //Modelica.SIunits.Enthalpy h_ea "actual outlet enthalpy";
  
  PBMedia.CO2_pT medium_isen "medium under isentropic process";
  
  //PBMedia.ThermodynamicState state_inlet;
  //PBMedia.ThermodynamicState state_outlet;
  
  //Modelica.SIunits.Enthalpy s_i "inlet entropy";
  
  Modelica.SIunits.AbsolutePressure p_out "fixed output pressure";
  
  Real eta "Efficiency of this turbine";
  
equation
  //fluid.cp = PropsSI("C", "T", inlet.T, "P", inlet.p, fluid.name);
  medium_in.state = PBMedia.setState_pTX(p = inlet.p, T = inlet.T);
  medium_isen.state = PBMedia.setState_psX(p = p_out, s = medium_in.s);
  
  
  
  //h_i = inlet.Medium.specificEnthalpy(inlet.Medium.state);// PropsSI("H", "T", inlet.T, "P", inlet.p, fluid.name);
  //s_i = PropsSI("S", "T", inlet.T, "P", inlet.p, fluid.name);
  
  //outlet
  //h_es = PropsSI("H", "P", outlet.p, "S", s_i, fluid.name); // no entropy increase
  //h_ea = h_i - eta * (h_i - h_es);
  medium_out.state = PBMedia.setState_phX(p = p_out, h = medium_in.h - eta * (medium_in.h - medium_isen.h));
  
  
  outlet.T = medium_out.T; // PropsSI("T", "H", h_ea, "P", outlet.p, fluid.name);  
  outlet.m_flow + inlet.m_flow = 0;
  outlet.p = medium_out.p;
  outlet.h_outflow = medium_out.h;
  inlet.h_outflow = inStream(outlet.h_outflow); //PBMedia.specificEnthalpy(state_inlet);
end Turbine;
