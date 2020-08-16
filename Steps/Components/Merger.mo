within Steps.Components;

model Merger
  extends Steps.Components.TwoPorts;
  
  import Modelica.SIunits;
  
  Steps.Interfaces.PBFluidPort_a inlet_merge(redeclare package Medium = PBMedia);
  PBMedia.CO2_pT medium_merge;
  
  parameter SIunits.Temperature T_init "initial temperature when no fluid from merge_in";
  parameter SIunits.Pressure p_init "initial pressure when no fluid from merge_in";

equation  
  
  medium_in.state = PBMedia.setState_pTX(p = inlet.p, T = inlet.T);
  medium_merge.state = PBMedia.setState_pTX(p = inlet_merge.p, T = inlet_merge.T);
  
  medium_out.state = PBMedia.setState_phX( 
    p = inlet.p, 
    h = (medium_in.h * inlet.m_flow +  medium_merge.h * inlet_merge.m_flow) / (inlet.m_flow + inlet_merge.m_flow)
  );  
 
  outlet.m_flow + inlet.m_flow + inlet_merge.m_flow = 0; 
  outlet.p = medium_out.p;
  outlet.T = medium_out.T; 
  outlet.h_outflow = medium_out.h;
  
  inlet.h_outflow = inStream(outlet.h_outflow);
  inlet_merge.h_outflow = inStream(outlet.h_outflow);

end Merger;
