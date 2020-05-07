within Steps.Components;

model Valve  "This class is basic valve"
  extends TwoPorts;
  
  parameter Modelica.SIunits.AbsolutePressure p_outlet "fixed outlet pressure";

equation
  medium_in.state = SCO2.setState_phX(p = inlet.p, h = inStream(inlet.h_outflow));
  medium_out.state = SCO2.setState_phX(p = p_outlet, h = inStream(inlet.h_outflow));
  
  outlet.m_flow + inlet.m_flow = 0;
  outlet.h_outflow = inStream(inlet.h_outflow);// medium.h;
  outlet.p = p_outlet;
  inlet.h_outflow = inStream(outlet.h_outflow);

  outlet.T = medium_out.T;
end Valve;
