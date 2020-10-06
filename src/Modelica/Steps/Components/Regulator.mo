within Steps.Components;

model Regulator
  extends Steps.Components.TwoPorts;
  
  import CP = Steps.Utilities.CoolProp;
      
  parameter Modelica.SIunits.AbsolutePressure p_init;
  //parameter Modelica.SIunits.Pressure P_E = -1.0;
  parameter  Modelica.SIunits.Temperature T_init; 
  
  parameter Modelica.SIunits.MassFlowRate m_flow_init;   
  
  //Modelica.Blocks.Interfaces.RealInput T_input;
      
equation  

  outlet.p = p_init;
  outlet.m_flow = - m_flow_init;
  outlet.h_outflow = CP.PropsSI("H", "P", p_init, "T", T_init, PBMedia.mediumName); 
  
  inlet.h_outflow = inStream(inlet.h_outflow);
end Regulator;
