within Steps.Components;

model Valve
  "This class is basic valve"
  extends TwoPorts;
  
  parameter Modelica.SIunits.AbsolutePressure p_outlet "fixed outlet pressure";
  //Steps.Media.SCO2.ThermodynamicState state_b ;
  //protected
    //Real h_i "specific enthalpy of inlet";  
    //Real h_e "specific enthalpy of outlet";
equation
  state_a = SCO2.setState_phX(p = inlet.p, h = inStream(inlet.h_outflow));
  state_b = SCO2.setState_phX(p = p_outlet, h = inStream(inlet.h_outflow));
  
  //h_i = PropsSI("H", "T", inlet.T, "P", inlet.p, fluid.name);
  //h_e = h_i;  
 
  //outlet.T = PropsSI("T", "H", h_e, "P", outlet.p, fluid.name);
  
  outlet.m_flow + inlet.m_flow = 0;
  outlet.h_outflow = inStream(inlet.h_outflow);// medium.h;
  outlet.p = p_outlet;
  inlet.h_outflow = inStream(outlet.h_outflow);
  //inlet.T = SCO2.temperature(state_a);
  outlet.T = SCO2.temperature(state_b);
end Valve;
