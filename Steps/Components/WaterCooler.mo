within Steps.Components;

model WaterCooler
  extends TwoPorts;
  
  import SI = Modelica.SIunits;

  parameter SI.Temperature T_cold_in;
  parameter Real pinch = 0.0; 
    
equation
  
  medium_in.state = PBMedia.setState_phX(p = inlet.p, h = inStream(inlet.h_outflow));
  
  outlet.p = inlet.p;
  medium_out.state = PBMedia.setState_pTX(p = outlet.p, T = T_cold_in + pinch);  

  outlet.m_flow + inlet.m_flow = 0;  
  outlet.h_outflow = medium_out.h; 
  inlet.h_outflow = inStream(outlet.h_outflow);
  outlet.T = medium_out.T; 
  
end WaterCooler;
