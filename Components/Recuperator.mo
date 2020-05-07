within Steps.Components;

model Recuperator
  "Recuperator to regenerate heat"
  extends Steps.Components.TwoPorts;
  
  Steps.Interfaces.PBFluidPort_a inlet_cool(redeclare package Medium = SCO2) "Recuperator inlet";
  //Steps.Interfaces.PBFluidPort_a inlet_hot(redeclare package Medium = SCO2) = inlet;
  Steps.Interfaces.PBFluidPort_b outlet_cool(redeclare package Medium = SCO2) "Recuperator outlet";
  //Steps.Interfaces.PBFluidPort_b outlet_hot(redeclare package Medium = SCO2) = outlet;
  //Steps.Fluid crec_fluid "Recuperator fluid";  
 
  Media.SCO2.CO2_pT medium_cool_in;
  Media.SCO2.CO2_pT medium_cool_out;
  //Media.SCO2.CO2_pT medium_hot_in = medium_in;
  //Media.SCO2.CO2_pT medium_hot_out = medium_out;
  
  //Intermedian variables for calculation of the max possible transfered heat
  Media.SCO2.CO2_pT medium_cool_max;
  Media.SCO2.CO2_pT medium_hot_max;
 
  //Modelica.SIunits.Enthalpy h_hot_i "enthalpy of hot inlet";
  //Modelica.SIunits.Enthalpy h_hot_e "enthalpy of hot exit";  
  //Modelica.SIunits.Enthalpy h_cool_i "enthalpy of cool inlet";
  //Modelica.SIunits.Enthalpy h_cool_e "enthalpy of cool outlet";
  
  //Modelica.SIunits.Temperature T_hot_i "Temperature of hot inlet";
  //Modelica.SIunits.Temperature T_cool_i "Temperature of cool inlet";
  
  //Modelica.SIunits.Pressure p_hot_i "Pressure of hot inlet";
  //Modelica.SIunits.Pressure p_cool_i "Pressure of cool inlet";  
  
  Real Q_actual;
  
  Real eta "heat exchange efficiency";
  
  //Real m2;
  
  //Real m_min = 1, m_max = 200;
  
  Real Q_max_hot "The maximum possible heat exchanged by hot pass";
  
  Real Q_max_cool "The maximum possible heat exchanged by cool pass";
  
equation
 
  //T_hot_i = RestrictValue(inlet.T, Constants.T_min, Constants.T_max); 
  medium_in.state = SCO2.setState_pTX(p = inlet.p, T = inlet.T);
  //T_hot_i = inlet.T;
  
  //T_cool_i = RestrictValue(inlet.T, Constants.T_min, Constants.T_max);
  medium_cool_in.state = SCO2.setState_pTX(p = inlet_cool.p, T = inlet_cool.T);
  //T_cool_i = crec_in.T;
    
  //p_hot_i = RestrictValue(inlet.p, p_min, p_max);
  //p_hot_i  = inlet.p;
  
  //p_cool_i = RestrictValue(crec_in.p, p_min, p_max);  
  //p_cool_i = crec_in.p;
  
  //h_hot_i = PropsSI("H", "T", T_hot_i, "P", p_hot_i, fluid.name);
  //h_cool_i = PropsSI("H", "T", T_cool_i, "P", p_cool_i, crec_fluid.name);
  
  //m2 =  RestrictValue(crec_in.m_flow, m_min, m_max);
    
  // Q_max_hot
  medium_cool_max.state = SCO2.setState_pTX(p = medium_in.p, T = medium_cool_in.T);
  Q_max_hot = inlet.m_flow * (medium_in.h - medium_cool_max.h);
  
  //Q_max_cool
  medium_hot_max.state = SCO2.setState_pTX(p = medium_cool_in.p, T = medium_in.T);
  Q_max_cool = inlet_cool.m_flow * (medium_hot_max.h - medium_cool_in.h); 
  
  Q_actual = eta * min(Q_max_hot,Q_max_cool); // maximum possible heat exchanged 
  //Q_actual = eta * Q1; // maximum possible heat exchanged 
  //h_hot_e = h_hot_i - Q_actual / inlet.m_flow;
  //h_cool_e = h_cool_i + Q_actual / m2;
  
  medium_cool_out.state = SCO2.setState_phX(p = inlet_cool.p, h = medium_cool_in.h + Q_actual / inlet_cool.m_flow);
  outlet_cool.T = medium_cool_out.T;// PropsSI("T", "P", crec_out.p, "H", h_cool_e, crec_fluid.name);
  outlet_cool.p = medium_cool_out.p;
  outlet_cool.m_flow + inlet_cool.m_flow = 0;
  outlet_cool.h_outflow = medium_cool_out.h;
  
  inlet_cool.h_outflow = inStream(outlet_cool.h_outflow);
  
  //crec_out.p = crec_in.p;
  //crec_out.m_flow = m2;
  
  medium_out.state = SCO2.setState_phX(p = inlet.p, h = medium_in.h - Q_actual / inlet.m_flow); 
  outlet.T = medium_out.T;
  outlet.p = medium_out.p;
  outlet.m_flow + inlet.m_flow = 0;
  outlet.h_outflow = medium_out.h;
  
  inlet.h_outflow = inStream(outlet.h_outflow);
  //outlet.T = PropsSI("T", "H", h_hot_e, "P", outlet.p, fluid.name);
  //outlet.p = inlet.p;
  //outlet.m_flow = inlet.m_flow;   

end Recuperator;
