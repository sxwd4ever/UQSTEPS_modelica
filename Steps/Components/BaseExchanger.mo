within Steps.Components;

model BaseExchanger "Base class for components such as heat exchangers,recuperators or heaters"
  replaceable package PBMedia = Steps.Media.SCO2;
  
  replaceable Steps.Interfaces.PBFluidPort_a inlet_hot(redeclare package Medium = PBMedia) "Inlet port, previous component";
  replaceable Steps.Interfaces.PBFluidPort_b outlet_hot(redeclare package Medium = PBMedia) "Outlet port, next component";
  replaceable Steps.Interfaces.PBFluidPort_a inlet_cool(redeclare package Medium = PBMedia) "Recuperator inlet";
  replaceable Steps.Interfaces.PBFluidPort_b outlet_cool(redeclare package Medium = PBMedia) "Recuperator outlet";
  
  // Common intermediate variables for states of inlet and outlet
  PBMedia.CO2_pT medium_cool_in;
  PBMedia.CO2_pT medium_cool_out;
  PBMedia.CO2_pT medium_hot_in;
  PBMedia.CO2_pT medium_hot_out;   
    
  parameter Boolean debug_mode = false;
  
  /*
  // Port type of the ports
  // To specify if the value(pressure, Temperature) of a port should be fixed by preset value or
  // determined by conncted port.
  parameter PortType PT_inlet_cool = PortType.free;
  parameter PortType PT_outlet_cool = PortType.free;
  parameter PortType PT_inlet_hot = PortType.free;
  parameter PortType PT_outlet_hot = PortType.free;
  */
end BaseExchanger;
