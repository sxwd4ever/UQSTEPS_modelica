within Steps.Components;

model FanCooler
  extends TwoPorts;
  
  import Steps.Interfaces.PortType;
  
  Modelica.SIunits.TemperatureDifference delta_T "Incoming temperature difference";  
  
  Modelica.SIunits.Temperature T_amb "ambinent temperature";
  
algorithm

  //assert(outlet.PT == PortType.free, "error port type");
  
equation
  // Note that, no matter what the CO2 inlet conditions are, the fan cooler CO2 
  // exit temperature is always the same for fixed ITD and ambient air temperature.  
  medium_in.state = PBMedia.setState_pTX(p = inlet.p, T = inlet.T);  
  medium_out.state = PBMedia.setState_pTX(p = inlet.p, T = T_amb + delta_T);
  
  outlet.T =  medium_out.T;
  outlet.p = inlet.p;
  outlet.m_flow + inlet.m_flow = 0;
  outlet.h_outflow = medium_out.h;
  
  inlet.h_outflow = inStream(outlet.h_outflow);

end FanCooler;
