within Steps.Components;

model ExchangerContainer "Base class for container components for exchangers,recuperators or heaters"
  replaceable package PBMedia = Steps.Media.SCO2;
    
  // connectors inside this container for counter port 
  replaceable Steps.Interfaces.PBFluidPort_a inlet_hot_in(redeclare package Medium = PBMedia) "hot Inlet port, previous component";
  replaceable Steps.Interfaces.PBFluidPort_b outlet_hot_in(redeclare package Medium = PBMedia) "hot Outlet port, next component";
  replaceable Steps.Interfaces.PBFluidPort_a inlet_cool_in(redeclare package Medium = PBMedia) "cool inlet port";
  replaceable Steps.Interfaces.PBFluidPort_b outlet_cool_in(redeclare package Medium = PBMedia) "Recuperator outlet";   

end ExchangerContainer;
