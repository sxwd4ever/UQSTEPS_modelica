within Steps.Components;

model Merger
  extends Steps.Components.TwoPorts;
  
  import Modelica.SIunits;
  //import Steps.Constants;
  
  //Steps.InPort merge_in;
  
  Steps.Interfaces.PBFluidPort_a inlet_merge(redeclare package Medium = SCO2);
  SCO2.CO2_pT medium_merge;
  
  //SIunits.Enthalpy h1(start=100.0), h2(start=100.0), h_out;  
  
  //SIunits.Temperature T1, T2;
  //SIunits.Pressure p1, p2;
  //Real m1(start=1.0), m2(start=1.0), t_e;
  
  parameter SIunits.Temperature T_init "initial temperature when no fluid from merge_in";
  parameter SIunits.Pressure p_init "initial pressure when no fluid from merge_in";
  //parameter Real m_init "initial mass flow of the merge_in";

equation  
  
//algorithm
  //T1 = RestrictValue(inlet.T, Constants.T_min, Constants.T_max);
  //T1 = inlet.T;
  //p1 = inlet.p;
  //m1 = inlet.m_flow;  
  medium_in.state = SCO2.setState_pTX(p = inlet.p, T = inlet.T);
  
  //T2 = RestrictValue(merge_in.T, Constants.T_min, Constants.T_max);
  medium_merge.state = SCO2.setState_pTX(p = inlet_merge.p, T = inlet_merge.T);
  
  
  //T2 = merge_in.T;
  //p2 = merge_in.p;
  //m2 = merge_in.m_flow; 

  //h1 = PropsSI("H", "T", T1, "P", p1, fluid.name);  
  //h2 = PropsSI("H", "T", T2, "P", p2, fluid.name); 
  
  medium_out.state = 
  SCO2.setState_phX( 
    p = inlet.p, 
    h = (medium_in.h * inlet.m_flow +  medium_merge.h * inlet_merge.m_flow) / (inlet.m_flow + inlet_merge.m_flow)
  );
  
  
  //t_e = PropsSI("T", "H", h_out, "P", outlet.p, fluid.name);
 
  outlet.m_flow + inlet.m_flow + inlet_merge.m_flow = 0; //m1 + m2;  
  outlet.p = medium_out.p;
  outlet.T = medium_out.T; 
  outlet.h_outflow = medium_out.h;
  
  inlet.h_outflow = inStream(outlet.h_outflow);
  inlet_merge.h_outflow = inStream(outlet.h_outflow);
  //outlet.T = t_e;  

end Merger;
