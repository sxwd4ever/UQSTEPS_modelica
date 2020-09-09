within Steps.Components;

model Sink 
  "Sink for component test with fixed inlet temperature and pressure"
  
  import CP = Steps.Utilities.CoolProp;  
  import Modelica.SIunits.Conversions.{from_bar, from_degC};
  
  replaceable package PBMedia = Steps.Media.SCO2;  

  // fixed outlet temperautre and pressure
  parameter Modelica.SIunits.AbsolutePressure p_inlet = 9 * 1e6;
  parameter Modelica.SIunits.Temperature T_inlet = from_degC(400);

  // fixed mass flow
  parameter Modelica.SIunits.MassFlowRate mdot_init = 8.3;     
  
  Modelica.SIunits.Temperature T;
 
  replaceable Steps.Interfaces.PBFluidPort_a inlet(redeclare package Medium = PBMedia) "Inlet port, previous component";
  
  parameter Boolean fix_state = true;
    
equation   

  if fix_state then
    inlet.p = p_inlet;
    inlet.m_flow = mdot_init;
    inlet.h_outflow = CP.PropsSI("H", "P", p_inlet, "T", T_inlet, PBMedia.mediumName);
    T = T_inlet;
  else    
    inlet.h_outflow = inStream(inlet.h_outflow); 
  end if; 
   
  T = CP.PropsSI("T", "P", inlet.p, "H", inStream(inlet.h_outflow), PBMedia.mediumName); 
end Sink;
