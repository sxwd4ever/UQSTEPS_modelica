within Steps.Components;

model WaterCooler
  extends TwoPorts;
  
  import SI = Modelica.SIunits;

  parameter SI.Temperature T_cool_in;
  parameter Real pinch = 0.0;
  
  //SCO2.CO2_pT medium_out;
  
  //protected
    //Real h_i "enthalpy of inlet";
    //Real h_e "enthalpy of outlet";
    //Real Q "Heat Exchanged";
    
equation
  
  state_a = SCO2.setState_phX(p = inlet.p, h = inStream(inlet.h_outflow));
  
  //outlet.T = cool_in.T + pinch;
  outlet.p = inlet.p;
  state_b = SCO2.setState_pTX(p = outlet.p, T = T_cool_in + pinch);
  
  //h_i = PropsSI("H", "T", inlet.T, "P", inlet.p, fluid.name);
  //h_e = PropsSI("H", "T", outlet.T, "P", outlet.p, fluid.name);
  
  //Q = inlet.m_flow * (inStream(inlet.h_outflow) - outlet.h_outflow);
  //cool_out.T = cool_in.T + Q / (cool_in.m_flow * cool_fluid.cp);  
  //cool_out.p = cool_in.p;
  
  outlet.m_flow + inlet.m_flow = 0;  
  outlet.h_outflow = SCO2.specificEnthalpy(state_b);
  inlet.h_outflow = inStream(outlet.h_outflow);
  //inlet.T = SCO2.temperature(state_a);
  outlet.T = SCO2.temperature(state_b);  
end WaterCooler;
