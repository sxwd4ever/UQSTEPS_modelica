within Steps.Components;

model EndPoint
  "EndPoint for component test with fixed inlet/outlet temperature and pressure"
  extends Steps.Components.TwoPorts;
  // fixed inlet temperature and pressure    
  parameter Modelica.SIunits.AbsolutePressure p_inlet;
  parameter  Modelica.SIunits.Temperature T_inlet; 
  
  // fixed outlet temperautre and pressure
  parameter Modelica.SIunits.AbsolutePressure p_outlet;
  parameter Modelica.SIunits.AbsolutePressure T_outlet;

	// fixed mass flow
  parameter Modelica.SIunits.MassFlowRate m_dot_flow;   
      
equation  
  medium_in.state = PBMedia.setState_pTX(p = p_inlet, T = T_inlet);
  medium_out.state = PBMedia.setState_pTX(p = p_outlet, T = T_outlet);
  
  outlet.p = medium_out.p;
  outlet.T = medium_out.T;
  outlet.h_outflow = medium_out.h;
  outlet.m_flow = - m_dot_flow;
  //outlet.m_flow + inlet.m_flow = 0;  

  inlet.p = medium_in.p;
  inlet.T = medium_in.T;  
  inlet.h_outflow = medium_in.h;
  inlet.m_flow = m_dot_flow;
  //inlet.h_outflow = inStream(outlet.h_outflow);

end EndPoint;
