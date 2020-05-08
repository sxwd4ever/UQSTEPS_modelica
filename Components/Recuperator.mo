within Steps.Components;

model Recuperator
  "Recuperator to regenerate heat"
  extends Steps.Components.TwoPorts;

  Real eta "heat exchange efficiency";

  Steps.Interfaces.PBFluidPort_a inlet_cool(redeclare package Medium = PBMedia) "Recuperator inlet";
  Steps.Interfaces.PBFluidPort_b outlet_cool(redeclare package Medium = PBMedia) "Recuperator outlet";
 
  Media.SCO2.CO2_pT medium_cool_in;
  Media.SCO2.CO2_pT medium_cool_out;
  
  //Intermedian variables for calculation of the max possible transfered heat
  Media.SCO2.CO2_pT medium_cool_max;
  Media.SCO2.CO2_pT medium_hot_max;
 
  Real Q_actual;  
  Real Q_max_hot "The maximum possible heat exchanged by hot pass";  
  Real Q_max_cool "The maximum possible heat exchanged by cool pass";
  
equation
 
  medium_in.state = PBMedia.setState_pTX(p = inlet.p, T = inlet.T);
  medium_cool_in.state = PBMedia.setState_pTX(p = inlet_cool.p, T = inlet_cool.T);
    
  // Q_max_hot
  medium_cool_max.state = PBMedia.setState_pTX(p = medium_in.p, T = medium_cool_in.T);
  Q_max_hot = inlet.m_flow * (medium_in.h - medium_cool_max.h);
  
  //Q_max_cool
  medium_hot_max.state = PBMedia.setState_pTX(p = medium_cool_in.p, T = medium_in.T);
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
  medium_out.state = PBMedia.setState_phX(p = inlet.p, h = medium_in.h - Q_actual / inlet.m_flow); 
  outlet.T = medium_out.T;
  outlet.p = medium_out.p;
  outlet.m_flow + inlet.m_flow = 0;
  outlet.h_outflow = medium_out.h;
  
  inlet.h_outflow = inStream(outlet.h_outflow);

end Recuperator;
