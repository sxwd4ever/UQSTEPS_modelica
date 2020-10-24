within Steps.Components;

model Peeker
  "Path peeker - to output the variable values at initialization stage"
  extends TwoPorts;
  
  import Steps.Interfaces.PortType;
  import Steps.Components.ThermoState;
  import Steps.Utilities.CoolProp.PropsSI;
  import Steps.Utilities.CoolProp.PrintPathState;
  
  //parameter Modelica.SIunits.TemperatureDifference delta_T = 18.0 "Incoming temperature difference";  
  
  ThermoState st_in(T = 0);
  ThermoState st_out(T = 0);
  
  parameter String path_name = "path";
  // enum log_level {DEBUG = 0, INFO = 1, ERR = 2, SERVE = 3, OFF = 4};
  parameter Integer log_level = 0;
  
equation

  st_in.p = inlet.p;
  st_in.h = inlet.h_outflow;
  st_in.mdot = inlet.m_flow;
  
  st_out.p = outlet.p;
  st_out.h = outlet.h_outflow;
  st_out.mdot = - outlet.m_flow;  
  
  PrintPathState(path_name, PBMedia.mediumName, st_in);
    
  outlet.p = inlet.p;
  outlet.m_flow + inlet.m_flow = 0;
  inlet.h_outflow = inStream(inlet.h_outflow); 
  inlet.h_outflow = outlet.h_outflow;

end Peeker;
