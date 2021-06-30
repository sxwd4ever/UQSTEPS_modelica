within Steps.Model;

record HEBoundaryCondition
  "Record for Heat Exchangers's PCHEBoundaryCondition, aligns with the struct definition in C"
    import Model.ThermoState;
    
    // error in transferring a string together with array of records within a record
    // so I move these two variables to the interface. 
    // String media_hot;
    // String media_cold;
  
    Model.ThermoState st_hot_in;
    Model.ThermoState st_cold_in;
    Model.ThermoState st_hot_out;
    Model.ThermoState st_cold_out;
  
end HEBoundaryCondition;
