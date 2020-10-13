within Steps.Components;

record ThermoState
  "Record for the ThermoState, align with the struct definition in C"
  
  Real T "Temperature, K";
  Real p "Pressure, pa";
  Real h "specific enthalpy, J/kg";
  Real mdot "mass flow rate kg/s";
  
end ThermoState;
