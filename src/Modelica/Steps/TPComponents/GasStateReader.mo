within Steps.TPComponents;

model GasStateReader   
  "Gas state reader which extends ThermoPower's StateReader_gas"
  extends ThermoPower.PowerPlants.HRSG.Components.StateReader_gas;
  
  Medium.Density rho "Temperature";
  Medium.SpecificEntropy s "Pressure";
  String fluid = Medium.mediumName "Fluid name";
    
equation   

  rho = gas.d;
  s = 0; // when used for thermal oil, use this trial value. 
  // s = gas.s;  
  // s = Medium.specificEntropy(gas);
  
end GasStateReader;
