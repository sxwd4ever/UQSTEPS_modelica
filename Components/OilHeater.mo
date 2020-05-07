within Steps.Components;

model OilHeater
  extends TwoPorts;
  import CP = Steps.Utilities.CoolProp;
  import SI = Modelica.SIunits;
  import HO = Steps.Media.Oil;
  
  //Modelica.Fluid.Interfaces.FluidPort_a hot_in(redeclare package Medium = HO);
  //Modelica.Fluid.Interfaces.FluidPort_b hot_out(redeclare package Medium = HO);
  
  //SCO2.CO2_pT medium_out "medium at outlet";
  //parameter Steps.Fluid hot_fluid;
    
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
    //Real h_fi "fluid inlet enthalpy";
    //Real h_fe "fluid outlet enthalpy";
    
equation
 
  state_a = SCO2.setState_phX(p = inlet.p, h = inStream(inlet.h_outflow));    
  cp_in = SCO2.specificHeatCapacityCp(state_a);
  C_min = min(inlet.m_flow * cp_in , m_dot_hot * cp_hot);
  
  // maximum possible exchanged heat
  Q = eta * C_min * (T_hot_in - SCO2.temperature(state_a));
  
  //medium_oil_out.state =
  
  //hot_out.T = T_hot_in - Q / (m_dot_hot * cp_hot);
  //hot_out.p = p_hot_in;
  
  //h_fi = PropsSI("H", "T", inlet.T, "P", inlet.p, fluid.name);
  outlet.p = inlet.p;
  state_b = SCO2.setState_phX(p = outlet.p, h = inStream(inlet.h_outflow) + Q / inlet.m_flow);
   
  //outlet.T = medium_out.T;
  outlet.m_flow + inlet.m_flow = 0;
  outlet.h_outflow = SCO2.specificEnthalpy(state_b);
  inlet.h_outflow = inStream(outlet.h_outflow);
  //inlet.T = SCO2.temperature(state_a);
  outlet.T = SCO2.temperature(state_b);
end OilHeater;
