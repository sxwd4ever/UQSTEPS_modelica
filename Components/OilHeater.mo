within Steps.Components;

model OilHeater
  extends TwoPorts;
  import CP = Steps.Utilities.CoolProp;
  import SI = Modelica.SIunits;
  import HO = Steps.Media.Oil;  
   
  import Modelica.Block.Math.Min;
  
  parameter Real eta "efficiency of the Oil heat exchanger";
  parameter SI.Temperature T_hot_in = 0 "input temperature of the hot inlet";
  parameter SI.AbsolutePressure p_hot_in = 10e5;
  parameter SI.MassFlowRate m_dot_hot "mass flow rate of the hot oil";
  parameter SI.SpecificHeatCapacity cp_hot = 3.0 * 1e3 "cp of hot oil";
  
  SI.SpecificHeatCapacity cp_in;
  
 //protected    
    Real C_min "c_min = min(m_flow_oil * cp_oil, m_flow_fluid * cp_fluid)";
    
    SI.HeatFlowRate Q "heat exchanged";
    
equation
 
  medium_in.state = SCO2.setState_phX(p = inlet.p, h = inStream(inlet.h_outflow));    
  cp_in = SCO2.specificHeatCapacityCp(medium_in.state);

  C_min = min(inlet.m_flow * cp_in , m_dot_hot * cp_hot);
  // maximum possible exchanged heat
  Q = eta * C_min * (T_hot_in - medium_in.T);
  
  outlet.p = inlet.p;
  medium_out.state = SCO2.setState_phX(p = outlet.p, h = inStream(inlet.h_outflow) + Q / inlet.m_flow);
 
  outlet.m_flow + inlet.m_flow = 0;
  outlet.h_outflow = medium_out.h;
  inlet.h_outflow = inStream(outlet.h_outflow);
  
  outlet.T = medium_out.T;
end OilHeater;
