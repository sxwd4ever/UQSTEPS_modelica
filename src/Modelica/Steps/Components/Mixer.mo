within Steps.Components;

model Mixer
  extends Steps.Components.TwoPorts;
  
  import Modelica.SIunits;
  
  Steps.Interfaces.PBFluidPort_a inlet_mix(redeclare package Medium = PBMedia);
  
  // parameter SIunits.Temperature T_init "initial temperature when no fluid from merge_in";
  // parameter SIunits.Pressure p_init "initial pressure when no fluid from merge_in";

  // Modelica.SIunits.SpecificEnthalpy h_mix "mixes h of the output stream";

equation  
  
  //h_mix = ;
  
  outlet.m_flow + inlet.m_flow + inlet_mix.m_flow = 0; 
  outlet.p = inlet.p;
  // outlet.p = inlet_mix.p;  
  //outlet.T = medium_out.T; 
 
  outlet.h_outflow = (inlet.h_outflow * inlet.m_flow + inlet_mix.h_outflow * inlet_mix.m_flow) / (inlet.m_flow + inlet_mix.m_flow);      
  inlet.h_outflow = inStream(inlet.h_outflow);  
  inlet_mix.h_outflow = inStream(inlet_mix.h_outflow);
 
  /*
  outlet.h_outflow * outlet.m_flow = inlet.h_outflow * inlet.m_flow + inlet_mix.h_outflow * inlet_mix.m_flow;
  inStream(inlet.h_outflow) = inStream(outlet.h_outflow);  
  inlet_mix.h_outflow = inStream(inlet_mix.h_outflow);
  */

end Mixer;
