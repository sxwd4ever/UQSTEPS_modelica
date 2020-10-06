within Steps.Components;

model Recuperator
  "Recuperator to regenerate heat"
  extends Steps.Components.BaseExchanger;
  
  import Steps.Utilities.Util.myAssert;

  Real eta "heat exchange efficiency";
  
  //Intermedian variables for calculation of the max possible transfered heat
  Steps.Media.SCO2.BaseProperties medium_cold_max;
  Steps.Media.SCO2.BaseProperties medium_hot_max;
 
  Real Q_actual;  
  Real Q_max_hot "The maximum possible heat exchanged by hot pass";  
  Real Q_max_cold "The maximum possible heat exchanged by cold pass";
  
equation
 
  medium_hot_in.state = PBMedia.setState_phX(p = inlet_hot.p, h = inStream(inlet_hot.h_outflow)); // T = inlet_hot.T);
  medium_cold_in.state = PBMedia.setState_phX(p = inlet_cold.p, h = inStream(inlet_cold.h_outflow)); // T = inlet_cold.T);
    
  // Q_max_hot
  medium_cold_max.state = PBMedia.setState_pTX(p = medium_hot_in.p, T = medium_cold_in.T);
  Q_max_hot = inlet_hot.m_flow * (medium_hot_in.h - medium_cold_max.h);
  
  //Q_max_cold
  medium_hot_max.state = PBMedia.setState_pTX(p = medium_cold_in.p, T = medium_hot_in.T);
  Q_max_cold = inlet_cold.m_flow * (medium_hot_max.h - medium_cold_in.h); 
  
  Q_actual = eta * min(Q_max_hot,Q_max_cold); // maximum possible heat exchanged 
  
  // for out_cold
  medium_cold_out.state = PBMedia.setState_phX(p = inlet_cold.p, h = medium_cold_in.h + Q_actual / inlet_cold.m_flow);
  
  myAssert(debug = true, val_test = medium_cold_out.h, min = 0, max = 1e5, name_val = "h_cold_out", val_ref = {medium_cold_in.h, inlet_cold.p, inlet_cold.m_flow, Q_actual}, name_val_ref = {"medium_cold_in.h", "inlet_cold.p", "inlet_cold.m_flow", "Q_actual"});  
  
  //outlet_cold.T = medium_cold_out.T;
  outlet_cold.p = medium_cold_out.p;
  outlet_cold.m_flow + inlet_cold.m_flow = 0;
  outlet_cold.h_outflow = medium_cold_out.h;  
  inlet_cold.h_outflow = inStream(outlet_cold.h_outflow);
  
  // for outlet_hot = outlet_main
  medium_hot_out.state = PBMedia.setState_phX(p = inlet_hot.p, h = medium_hot_in.h - Q_actual / inlet_hot.m_flow); 
  //outlet_hot.T = medium_hot_out.T;
  outlet_hot.p = medium_hot_out.p;
  outlet_hot.m_flow + inlet_hot.m_flow = 0;
  outlet_hot.h_outflow = medium_hot_out.h;
  
  inlet_hot.h_outflow = inStream(outlet_hot.h_outflow);

end Recuperator;
