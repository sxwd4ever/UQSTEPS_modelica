within Steps.Components;

record PCHEBoundaryCondition
  "Record for PCHE's PCHEBoundaryCondition, aligns with the struct definition in C"
    // error in transferring a string together with array of records within a record
    // so I move these two variables to the interface. 
    // String media_hot;
    // String media_cold;
  
    ThermoState st_hot_in;
    ThermoState st_cold_in;
    ThermoState st_hot_out;
    ThermoState st_cold_out;
  
end PCHEBoundaryCondition;
