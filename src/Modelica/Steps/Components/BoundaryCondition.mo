within Steps.Components;

record BoundaryCondition
  "Record for Kim correlation coefficients, aligns with the struct definition in C"
  
    ThermoState st_hot_in;
    ThermoState st_cold_in;
    ThermoState st_hot_out;
    ThermoState st_cold_out;
    
    String media_hot;
    String media_cold;
  
end BoundaryCondition;
