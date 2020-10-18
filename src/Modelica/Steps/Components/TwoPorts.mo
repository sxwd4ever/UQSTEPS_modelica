within Steps.Components;

model TwoPorts
  "Abstracted classes for the components with two ports"
  replaceable package PBMedia = Steps.Media.SCO2;
  
  import Steps.Interfaces.PortType;
  import Steps.Utilities.CoolProp.PropsSI;
 
  replaceable Steps.Interfaces.PBFluidPort_a inlet(redeclare package Medium = PBMedia) "Inlet port, previous component";
  replaceable Steps.Interfaces.PBFluidPort_b outlet(redeclare package Medium = PBMedia) "Outlet port, next component";
   
  parameter Boolean debug_mode = false;
  
  Modelica.SIunits.Temperature T_inlet;
  
  Modelica.SIunits.Temperature T_outlet;
  
equation

  T_inlet = PropsSI("T","P", inlet.p, "H", inlet.h_outflow, PBMedia.mediumName); 
  
  T_outlet = PropsSI("T", "P", outlet.p, "H", outlet.h_outflow, PBMedia.mediumName);

end TwoPorts;
