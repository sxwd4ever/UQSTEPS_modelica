within Steps.Components;

model Sink 
  "Sink for component test with fixed inlet temperature and pressure"
  
  import CP = Steps.Utilities.CoolProp;  
  replaceable package PBMedia = Steps.Media.SCO2;
  
  // fixed inlet temperature and pressure    
  //parameter Modelica.SIunits.AbsolutePressure p_inlet;
  //parameter  Modelica.SIunits.Temperature T_inlet; 
  
  // fixed outlet temperautre and pressure
  parameter Modelica.SIunits.AbsolutePressure p_inlet(start = 9 * 1e6);
  parameter Modelica.SIunits.Temperature T_inlet(start = 400);

  // fixed mass flow
  parameter Modelica.SIunits.MassFlowRate mdot_init(start = 8.3);     
  
  
  //parameter Boolean IsSource = true;
  
  //replaceable package PBMedia = Steps.Media.SCO2;
 
  replaceable Steps.Interfaces.PBFluidPort_a inlet(redeclare package Medium = PBMedia) "Inlet port, previous component";
  
  parameter Boolean fix_state = true;
  
  //PBMedia.CO2_pT medium_in;
  //replaceable Steps.Interfaces.PBFluidPort_b outlet(redeclare package Medium = PBMedia) "Outlet port, next component";

initial equation
  // use initial equation to set initial fixed value for the cross variable in connector
  //medium_in.state = PBMedia.setState_pTX(p = inlet.p, T = inStream(inlet.h_outflow));
  //inlet.T = T_inlet;
  //inlet.p = p_inlet;
  /*
  if fix_state then
    inlet.T = T_inlet;
    inlet.p = p_inlet;
    inlet.m_flow = mdot_init;
    inlet.h_outflow = CP.PropsSI("H", "P", p_inlet, "T", T_inlet, PBMedia.mediumName);
  end if;
  */
  /*
  inlet.m_flow = 
  outlet.T = T_outlet;
  outlet.p = p_outlet;
  outlet.m_flow = - m_dot_flow;
  */
  /*
  if IsSource == true then
    outlet.m_flow = - m_dot_flow;
  end if;
             
algorithm
  outlet.p := p_outlet;
  outlet.T := T_outlet;
  */ 
equation  
  
  //medium_in.state = PBMedia.setState_pTX(p = p_inlet, T = T_inlet);
  //inlet.T = T_inlet;//T_inlet;
  //inlet.p = p_inlet;

  if fix_state then
    inlet.T = T_inlet;
    inlet.p = p_inlet;
    inlet.m_flow = mdot_init;
    inlet.h_outflow = CP.PropsSI("H", "P", p_inlet, "T", T_inlet, PBMedia.mediumName);
  else
    inlet.h_outflow = CP.PropsSI("H", "P", inlet.p, "T", inlet.T, PBMedia.mediumName);
  end if;

  
  /*
  medium_in.state = PBMedia.setState_phX(p = inlet.p, h = inStream(inlet.h_outflow));
  medium_out.state = PBMedia.setState_pTX(p = p_outlet, T = T_outlet); 
  
  //outlet.m_flow = - m_dot_flow;
  //outlet.h_outflow = CP.PropsSI("H", "P", outlet.p, "T", outlet.T, "CO2");    
  
  outlet.p = medium_out.p;
  outlet.T = medium_out.T;
  outlet.m_flow = - m_dot_flow;
  outlet.h_outflow = medium_out.h;
  inlet.h_outflow = inStream(outlet.h_outflow);
    
  //outlet.m_flow + inlet.m_flow = 0;  

  //inlet.p = medium_in.p;
  //inlet.T = medium_in.T;  
  //inlet.h_outflow = medium_in.h;
  //inlet.m_flow = m_dot_flow;
  //inlet.h_outflow = inStream(outlet.h_outflow);
  */
end Sink;
