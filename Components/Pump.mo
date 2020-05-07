within Steps.Components;

model Pump
  "class for a Pump"
  extends TwoPorts;
  
  import SI = Modelica.SIunits;  
  
  //Intermedia variables
  //SCO2.CO2_pT medium_isen "ideal medium under isentropic process";
  //SCO2.CO2_pT medium_out "actual medium at outlet"; 
  
  input Real eta  "efficiency of this pump";   
  input SI.AbsolutePressure p_outlet "fixed outlet pressure of pump"; 
  
  SCO2.ThermodynamicState state_isen;
  SI.SpecificEnthalpy h_i;
  SI.SpecificEnthalpy h_es;
     
equation
 
  state_a = SCO2.setState_phX(p = inlet.p, h = inStream(inlet.h_outflow));
  state_isen = SCO2.setState_psX(p = p_outlet, s = SCO2.specificEntropy(state_a)); 
   
  h_i = SCO2.specificEnthalpy(state_a);
  h_es = SCO2.specificEnthalpy(state_isen);
  
  state_b = SCO2.setState_phX(p = p_outlet, h = (h_i - (h_i - h_es) / eta));
  
  inlet.m_flow + outlet.m_flow = 0;
  outlet.p = state_b.p;
  outlet.h_outflow = SCO2.specificEnthalpy(state_b);
  inlet.h_outflow = inStream(outlet.h_outflow);
  //inlet.T = SCO2.temperature(state_a);
  outlet.T = SCO2.temperature(state_b);
  //fluid.cp = PropsSI("C", "T", inlet.T, "P", inlet.p, fluid.name);
 /* 
  h_i = PropsSI("H", "T", inlet.T, "P", inlet.p, fluid.name) "inlet enthalpy";
  s_i =  PropsSI("S", "T", inlet.T, "P", inlet.p, fluid.name) "inlet entropy";
  
  h_es = PropsSI("H", "P", outlet.p, "S", s_i, fluid.name) "idealized outlet enthalpy";
  h_ea = h_i - (h_i - h_es) / eta "actual outlet enthalpy";       
  
  outlet.T = PropsSI("T", "H", h_ea, "P", outlet.p, fluid.name);
  outlet.m_flow = inlet.m_flow;
  */
end Pump;
