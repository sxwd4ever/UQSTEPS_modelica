within Steps.Components;

model Regulator
  extends Steps.Components.TwoPorts;
      
  parameter Modelica.SIunits.AbsolutePressure p_init;
  //parameter Modelica.SIunits.Pressure P_E = -1.0;
  parameter  Modelica.SIunits.Temperature T_init; 
  
  parameter Modelica.SIunits.MassFlowRate m_flow_init;   
  
  //Modelica.Blocks.Interfaces.RealInput T_input;
      
equation  
  state_a = SCO2.setState_phX(p = inlet.p, h = inStream(inlet.h_outflow));
  state_b = SCO2.setState_pTX(p = p_init, T = T_init);
  //medium.p = p_init;
  //medium.T = T_init; // .state = SCO2.setState_pTX(p = p_init, T = T_init);
  //inlet.p = p_init;
  //inlet.T = T_init;
  //inlet.h_outflow = - medium.h;
  
  outlet.p = state_b.p;
  outlet.m_flow = - m_flow_init;
  //outlet.m_flow + inlet.m_flow = 0;
  outlet.h_outflow = SCO2.specificEnthalpy(state_b);
  outlet.T = SCO2.temperature(state_a);
  
  inlet.h_outflow = SCO2.specificEnthalpy(state_a);
  //inlet.T = SCO2.temperature(state_b);
end Regulator;
