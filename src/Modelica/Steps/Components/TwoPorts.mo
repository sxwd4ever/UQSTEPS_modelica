within Steps.Components;

model TwoPorts
  "Abstracted classes for the components with two ports"
  replaceable package PBMedia = Steps.Media.SCO2;
  
  import Steps.Interfaces.PortType;
  import Steps.Utilities.CoolProp.PropsSI;
 
  replaceable Steps.Interfaces.PBFluidPort_a inlet(redeclare package Medium = PBMedia) "Inlet port, previous component";
  replaceable Steps.Interfaces.PBFluidPort_b outlet(redeclare package Medium = PBMedia) "Outlet port, next component";
   
  parameter Boolean debug_mode = false;
  
  // runtime T of inlet/outlet for debug purpose
  Modelica.SIunits.Temperature T_inlet_rt;  
  Modelica.SIunits.Temperature T_outlet_rt;
  
equation

  T_inlet_rt = PropsSI("T","P", inlet.p, "H", inlet.h_outflow, PBMedia.mediumName); 
  
  T_outlet_rt = PropsSI("T", "P", outlet.p, "H", outlet.h_outflow, PBMedia.mediumName);

end TwoPorts;
