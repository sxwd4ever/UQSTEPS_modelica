within Steps.Components;

model Regulator
  extends Steps.Components.TwoPorts;
      
  parameter Modelica.SIunits.AbsolutePressure p_init;
  //parameter Modelica.SIunits.Pressure P_E = -1.0;
  parameter  Modelica.SIunits.Temperature T_init; 
  
  parameter Modelica.SIunits.MassFlowRate m_flow_init;   
  
  //Modelica.Blocks.Interfaces.RealInput T_input;
      
equation  
  medium_in.state = PBMedia.setState_phX(p = inlet.p, h = inStream(inlet.h_outflow));
  medium_out.state = PBMedia.setState_pTX(p = p_init, T = T_init);
  
  outlet.p = medium_out.p;
  outlet.m_flow = - m_flow_init;
  outlet.h_outflow = medium_out.h;
  outlet.T = medium_in.T;
  
  inlet.h_outflow = medium_in.h;
end Regulator;
