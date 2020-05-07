within Steps.Components;

model Turbine
  extends Steps.Components.TwoPorts;
  
  replaceable package PBMedia = Steps.Media.SCO2;
  
  Modelica.SIunits.Enthalpy h_i;
  Modelica.SIunits.Enthalpy h_es "idealized outlet enthalpy";
  Modelica.SIunits.Enthalpy h_ea "actual outlet enthalpy";
  
  PBMedia.ThermodynamicState state_inlet;
  PBMedia.ThermodynamicState state_outlet;
  
  Modelica.SIunits.Enthalpy s_i "inlet entropy";
  
  Real eta "Efficiency of this turbine";
  
equation
  //fluid.cp = PropsSI("C", "T", inlet.T, "P", inlet.p, fluid.name);
  state_inlet = PBMedia.setState_pTX(inlet.p, inlet.T);
  
  inlet.h_outflow = PBMedia.specificEnthalpy(state_inlet);
  
  //h_i = inlet.Medium.specificEnthalpy(inlet.Medium.state);// PropsSI("H", "T", inlet.T, "P", inlet.p, fluid.name);
  s_i = PropsSI("S", "T", inlet.T, "P", inlet.p, fluid.name);
  
  //outlet
  h_es = PropsSI("H", "P", outlet.p, "S", s_i, fluid.name); // no entropy increase
  h_ea = h_i - eta * (h_i - h_es);
  
  outlet.T = PropsSI("T", "H", h_ea, "P", outlet.p, fluid.name);
  outlet.m_flow = inlet.m_flow;

end Turbine;
