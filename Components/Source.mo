within Steps.Components;

model Source 
  "EndPoint for component test with fixed inlet/outlet temperature and pressure"
  //extends Steps.Components.TwoPorts;
  
  import CP = Steps.Utilities.CoolProp;  
  replaceable package PBMedia = Steps.Media.SCO2;
  
  // fixed inlet temperature and pressure    
  //parameter Modelica.SIunits.AbsolutePressure p_inlet;
  //parameter  Modelica.SIunits.Temperature T_inlet; 
  
  // fixed outlet temperautre and pressure
  parameter Modelica.SIunits.AbsolutePressure p_outlet(start = 9 * 1e6);
  parameter Modelica.SIunits.AbsolutePressure T_outlet(start = 400);

  // fixed mass flow
  parameter Modelica.SIunits.MassFlowRate m_dot_flow(start = 8.3);
  //PBMedia.CO2_pT medium_out;    
  
  replaceable Steps.Interfaces.PBFluidPort_b outlet(redeclare package Medium = PBMedia) "Outlet port, next component";  
  
  //parameter Boolean IsSource = true;
  
  //replaceable package PBMedia = Steps.Media.SCO2;
 
  //replaceable Steps.Interfaces.PBFluidPort_a inlet(redeclare package Medium = PBMedia) "Inlet port, previous component";
  //replaceable Steps.Interfaces.PBFluidPort_b outlet(redeclare package Medium = PBMedia) "Outlet port, next component";

initial equation
  // use initial equation to set initial fixed value for the cross variable in connector
  outlet.p = p_outlet;
  outlet.T = T_outlet;
  outlet.m_flow = - m_dot_flow;
  /*
  if IsSource == true then
    outlet.m_flow = - m_dot_flow;
  end if;
             
algorithm
  outlet.p := p_outlet;
  outlet.T := T_outlet;
  */ 
equation  

  //medium_in.state = PBMedia.setState_phX(p = inlet.p, h = inStream(inlet.h_outflow));
  //medium_out.state = PBMedia.setState_pTX(p = p_outlet, T = T_outlet); 
  
  //outlet.m_flow = - m_dot_flow;
  //outlet.h_outflow = CP.PropsSI("H", "P", outlet.p, "T", outlet.T, "CO2");    
  
  outlet.p = p_outlet;
  outlet.T = T_outlet;
  outlet.m_flow = - m_dot_flow;
  outlet.h_outflow = - CP.PropsSI("H", "P", p_outlet, "T", T_outlet, PBMedia.mediumName);
  //inlet.h_outflow = inStream(outlet.h_outflow);
    
  //outlet.m_flow + inlet.m_flow = 0;  

  //inlet.p = medium_in.p;
  //inlet.T = medium_in.T;  
  //inlet.h_outflow = medium_in.h;
  //inlet.m_flow = m_dot_flow;
  //inlet.h_outflow = inStream(outlet.h_outflow);

end Source;
