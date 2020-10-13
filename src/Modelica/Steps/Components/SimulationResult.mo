within Steps.Components;

record SimulationResult
  "Record for the ThermoState, align with the struct definition in C"
  /*
    ThermoState st_hot_in;
    ThermoState st_cold_in;
    ThermoState st_hot_out;
    ThermoState st_cold_out;
  */
  
  Real T_hot[2];
  Real p_hot[2];
  Real T_cold[2];
  Real p_cold[2];
  
end SimulationResult;
