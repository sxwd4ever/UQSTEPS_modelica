within Steps.Model;

record HEConfig
  "Record for the simulation param, align with the struct definition in C"
  
    Real rho_mcm "Metal heat capacity per unit volume [J/m^3.K]";
    Real lambda "Thermal conductivity of the metal (density by specific heat capacity)";
    Boolean SSInit "Steady-state initialization";

    EntityConfig cfg_hot "hot side configuration";
    EntityConfig cfg_cold "cold side configuration";
    EntityConfig cfg_tube "tube configuration";
  
end HEConfig;
