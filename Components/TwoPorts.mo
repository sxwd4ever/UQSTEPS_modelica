within Steps.Components;

model TwoPorts
  "Abstracted classes for the components with two ports"
  replaceable package PBMedia = Steps.Media.SCO2;
 
  replaceable Steps.Interfaces.PBFluidPort_a inlet(redeclare package Medium = PBMedia) "Inlet port, previous component";
  replaceable Steps.Interfaces.PBFluidPort_b outlet(redeclare package Medium = PBMedia) "Outlet port, next component";
  
  // Common intermediate variables for states of inlet and outlet
  PBMedia.CO2_pT medium_in;
  PBMedia.CO2_pT medium_out;
  
end TwoPorts;
