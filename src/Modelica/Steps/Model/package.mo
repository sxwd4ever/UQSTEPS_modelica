within Steps;

package Model
  "Data model definition for consistent parameter configuration"
  extends Modelica.Icons.Package;

  record EntityConfig
    "Record for the ThermoState, align with the struct definition in C"
    
    ThermoState st_start "thermo state for start point";
    ThermoState st_nom "nominal state values";
    EntityGeoParam geo "Geometry values";   
  
  end EntityConfig;

  record EntityGeoParam
    "Record for PCHE's PCHEBoundaryCondition, aligns with the struct definition in C"      
      
      // error in transferring a string together with array of records within a record
      // so I move these two variables to the interface. 
      // String media_hot;
      // String media_cold;
    
      ThermoState st_hot_in;
      ThermoState st_cold_in;
      ThermoState st_hot_out;
      ThermoState st_cold_out;
    
  end EntityGeoParam;

  record ThermoState
    "Record for the ThermoState, align with the struct definition in C"
    
    Real T "Temperature, K";
    Real p "Pressure, pa";
    Real h "specific enthalpy, J/kg";
    Real mdot "mass flow rate kg/s";
    
    // optional id, to specify which point this state is assigned to 
    // because of the mapping constraintns of Modelica and external functions
    // integer id is used instead of string(char *) name
    // more detailed naming specificaiton should be provided
    // hot_in=0, cold_in=1, hot_out=2, cold_out=3
    Integer id = -1; 
  
  end ThermoState;

end Model;
