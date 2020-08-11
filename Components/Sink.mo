within Steps.Components;

model Sink 
  "Sink for component test with fixed inlet temperature and pressure"
  
  import CP = Steps.Utilities.CoolProp;  
  replaceable package PBMedia = Steps.Media.SCO2;  

  // fixed outlet temperautre and pressure
  parameter Modelica.SIunits.AbsolutePressure p_inlet(start = 9 * 1e6);
  parameter Modelica.SIunits.Temperature T_inlet(start = 400);

  // fixed mass flow
  parameter Modelica.SIunits.MassFlowRate mdot_init(start = 8.3);     
  
  Modelica.SIunits.Temperature T;
 
  replaceable Steps.Interfaces.PBFluidPort_a inlet(redeclare package Medium = PBMedia) "Inlet port, previous component";
  
  parameter Boolean fix_state = true;
    
equation   

  if fix_state then
    inlet.p = p_inlet;
    inlet.m_flow = mdot_init;
    inlet.h_outflow = CP.PropsSI("H", "P", p_inlet, "T", T_inlet, PBMedia.mediumName);
  end if;
  
  T = CP.PropsSI("T", "P", inlet.p, "H", inStream(inlet.h_outflow), PBMedia.mediumName);  

end Sink;
