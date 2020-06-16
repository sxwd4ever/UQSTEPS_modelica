within Steps.Components;

model Splitter
  extends Steps.Components.TwoPorts;
  
  //Steps.OutPort split_out "the splitted output port";
  Steps.Interfaces.PBFluidPort_b outlet_split(redeclare package Medium = PBMedia);
  PBMedia.CO2_pT medium_split "splitted medium";
  
  input Real split_ratio(start = 0.6) "mass split ratio for the input work fluid";
  
equation
  
  medium_in.state = PBMedia.setState_pTX(p = inlet.p, T = inlet.T);
  
  // split outlet
  medium_out.state = medium_in.state;
  
  outlet.m_flow = - inlet.m_flow * (1 - split_ratio); // first output
  outlet.T = medium_out.T;
  outlet.p = medium_out.p;
  outlet.h_outflow = medium_out.h;
  
  inlet.h_outflow = inStream(outlet.h_outflow);
  
  medium_split.state = medium_in.state;
  outlet_split.m_flow = - inlet.m_flow * split_ratio; // second output
  outlet_split.T = medium_split.T;
  outlet_split.p = medium_split.p;
  outlet_split.h_outflow = medium_out.h;      
    
end Splitter;
