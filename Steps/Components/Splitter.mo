within Steps.Components;

model Splitter
  extends Steps.Components.TwoPorts;
  
  //Steps.OutPort split_out "the splitted output port";
  Steps.Interfaces.PBFluidPort_b outlet_split(redeclare package Medium = PBMedia);
  
  parameter Real split_ratio(start = 0.6) "mass split ratio for the input work fluid";
  
equation  

  inlet.h_outflow = inStream(inlet.h_outflow);  

  outlet.p = inlet.p;
  outlet.h_outflow = inlet.h_outflow;
  outlet.m_flow = - inlet.m_flow * (1 - split_ratio); // first output
  
  outlet_split.p = inlet.p;
  outlet_split.h_outflow = inlet.h_outflow;      
  outlet_split.m_flow =  - inlet.m_flow * split_ratio; // second output 
    
end Splitter;
