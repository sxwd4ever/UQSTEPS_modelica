within Steps.Components;

model Recuperator
  "Recuperator to regenerate heat"
  extends Steps.Components.BaseExchanger;

  Real eta "heat exchange efficiency";
  
  //Intermedian variables for calculation of the max possible transfered heat
  Media.SCO2.CO2_pT medium_cool_max;
  Media.SCO2.CO2_pT medium_hot_max;
 
  Real Q_actual;  
  Real Q_max_hot "The maximum possible heat exchanged by hot pass";  
  Real Q_max_cool "The maximum possible heat exchanged by cool pass";
  
equation
 
  medium_hot_in.state = PBMedia.setState_pTX(p = inlet_hot.p, T = inlet_hot.T);
  medium_cool_in.state = PBMedia.setState_pTX(p = inlet_cool.p, T = inlet_cool.T);
    
  // Q_max_hot
  medium_cool_max.state = PBMedia.setState_pTX(p = medium_hot_in.p, T = medium_cool_in.T);
  Q_max_hot = inlet_hot.m_flow * (medium_hot_in.h - medium_cool_max.h);
  
  //Q_max_cool
  medium_hot_max.state = PBMedia.setState_pTX(p = medium_cool_in.p, T = medium_hot_in.T);
  Q_max_cool = inlet_cool.m_flow * (medium_hot_max.h - medium_cool_in.h); 
  
  Q_actual = eta * min(Q_max_hot,Q_max_cool); // maximum possible heat exchanged 
  
  // for out_cool
  medium_cool_out.state = PBMedia.setState_phX(p = inlet_cool.p, h = medium_cool_in.h + Q_actual / inlet_cool.m_flow);
  outlet_cool.T = medium_cool_out.T;
  outlet_cool.p = medium_cool_out.p;
  outlet_cool.m_flow + inlet_cool.m_flow = 0;
  outlet_cool.h_outflow = medium_cool_out.h;  
  inlet_cool.h_outflow = inStream(outlet_cool.h_outflow);
  
  // for outlet_hot = outlet_main
  medium_hot_out.state = PBMedia.setState_phX(p = inlet_hot.p, h = medium_hot_in.h - Q_actual / inlet_hot.m_flow); 
  outlet_hot.T = medium_hot_out.T;
  outlet_hot.p = medium_hot_out.p;
  outlet_hot.m_flow + inlet_hot.m_flow = 0;
  outlet_hot.h_outflow = medium_hot_out.h;
  
  inlet_hot.h_outflow = inStream(outlet_hot.h_outflow);

end Recuperator;
