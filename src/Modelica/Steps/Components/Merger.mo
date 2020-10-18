within Steps.Components;

model Merger
  extends Steps.Components.TwoPorts;
  
  import Modelica.SIunits;
  
  Steps.Interfaces.PBFluidPort_a inlet_merge(redeclare package Medium = PBMedia);
  
  // parameter SIunits.Temperature T_init "initial temperature when no fluid from merge_in";
  // parameter SIunits.Pressure p_init "initial pressure when no fluid from merge_in";

  Modelica.SIunits.SpecificEnthalpy h_mix "mixes h of the output stream";

equation  
  
  h_mix = (inlet.h_outflow * inlet.m_flow + inlet_merge.h_outflow * inlet_merge.m_flow) / (inlet.m_flow + inlet_merge.m_flow);
  
  outlet.m_flow + inlet.m_flow + inlet_merge.m_flow = 0; 
  outlet.p = inlet.p;
  //outlet.T = medium_out.T; 
  outlet.h_outflow = h_mix;  
  
  inlet.h_outflow = inStream(inlet.h_outflow);
  inlet_merge.h_outflow = inStream(inlet_merge.h_outflow);

end Merger;
