within Steps.Components;

model Source 
  "EndPoint for component test with fixed inlet/outlet temperature and pressure"
  //extends Steps.Components.TwoPorts;
  
  import CP = Steps.Utilities.CoolProp;  
  import Modelica.SIunits.Conversions.{from_bar, from_degC};
  import Modelica.Blocks.Interfaces.RealInput;
  
  replaceable package PBMedia = Steps.Media.SCO2;
  
  // fixed outlet temperautre and pressure
  parameter Modelica.SIunits.AbsolutePressure p_outlet = 9 * 1e6;
  parameter Modelica.SIunits.Temperature T_outlet = from_degC(400);
  
  // fixed mass flow
  parameter Modelica.SIunits.MassFlowRate mdot_init = 8.3;
  
  parameter Boolean fix_state = true "if use this component as boundary condition";
  
  Modelica.SIunits.Temperature T;
  
  replaceable Steps.Interfaces.PBFluidPort_b outlet(redeclare package Medium = PBMedia, p.start = p_outlet) "Outlet port, next component";  
 
equation  
  //determine the boundary condition
  if fix_state then 
    outlet.p = p_outlet;
    outlet.m_flow = - mdot_init;  
    // outlet.h_outflow = max(CP.PropsSI("H", "P", p_outlet, "T", T_outlet, PBMedia.mediumName), h_out_update);
    outlet.h_outflow = CP.PropsSI("H", "P", p_outlet, "T", T_outlet, PBMedia.mediumName);
    T = T_outlet;
  else     
    //outlet.h_outflow = inStream(outlet.h_outflow); 
    T = CP.PropsSI("T", "P", p_outlet, "H", outlet.h_outflow, PBMedia.mediumName);
  end if;     
  
end Source;
