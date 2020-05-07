within Steps.Components;

model TwoPorts
  "Abstracted classes for the components with two ports"
  replaceable package SCO2 = Steps.Media.SCO2;
 
  replaceable Steps.Interfaces.PBFluidPort_a inlet(redeclare package Medium = SCO2) "Inlet port, previous component";
  replaceable Steps.Interfaces.PBFluidPort_b outlet(redeclare package Medium = SCO2) "Outlet port, next component";
  
  //SCO2.CO2_pT medium "medim in this component";  
  SCO2.ThermodynamicState state_a;
  SCO2.ThermodynamicState state_b;

end TwoPorts;
