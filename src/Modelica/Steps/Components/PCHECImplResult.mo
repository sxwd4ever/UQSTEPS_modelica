within Steps.Components;

record PCHECImplResult
  "Record for Result of PCHE's CImpl, align with the struct definition in C"
  
  Real p_hot "hot side calculated pressuare, could be in/out"; 
  Real h_hot "hot side calculated specific enthalphy";
  
  Real p_cold "cold side calculated pressuare, could be in/out"; 
  Real h_cold "cold side calculated specific enthalphy, could be in/out"; 
  
end PCHECImplResult;
