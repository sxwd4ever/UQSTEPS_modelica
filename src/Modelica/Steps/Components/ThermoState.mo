within Steps.Components;

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
