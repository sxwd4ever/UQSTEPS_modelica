within Steps.Components;

model TwoPorts
  "Abstracted classes for the components with two ports"
  replaceable package PBMedia = Steps.Media.SCO2;
  
  import Steps.Interfaces.PortType;
 
  replaceable Steps.Interfaces.PBFluidPort_a inlet(redeclare package Medium = PBMedia) "Inlet port, previous component";
  replaceable Steps.Interfaces.PBFluidPort_b outlet(redeclare package Medium = PBMedia) "Outlet port, next component";
   
  parameter Boolean debug_mode = false;
  
end TwoPorts;
