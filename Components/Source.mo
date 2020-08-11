within Steps.Components;

model Source 
  "EndPoint for component test with fixed inlet/outlet temperature and pressure"
  //extends Steps.Components.TwoPorts;
  
  import CP = Steps.Utilities.CoolProp;  
  replaceable package PBMedia = Steps.Media.SCO2;
  
  // fixed outlet temperautre and pressure
  parameter Modelica.SIunits.AbsolutePressure p_outlet(start = 9 * 1e6);
  parameter Modelica.SIunits.Temperature T_outlet(start = 400);

  // fixed mass flow
  parameter Modelica.SIunits.MassFlowRate mdot_init(start = 8.3);
  
  parameter Boolean fix_state = true;
  
  Modelica.SIunits.Temperature T;
  
  replaceable Steps.Interfaces.PBFluidPort_b outlet(redeclare package Medium = PBMedia) "Outlet port, next component";  
 
equation  
  //determine the boundary condition
  if fix_state then 
    outlet.p = p_outlet;

    outlet.m_flow = - mdot_init;  
    outlet.h_outflow = CP.PropsSI("H", "P", p_outlet, "T", T_outlet, PBMedia.mediumName);
  end if;  
  
  T = CP.PropsSI("T", "P", outlet.p, "H", inStream(outlet.h_outflow), PBMedia.mediumName);

end Source;
