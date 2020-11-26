within Steps.Components;

model WaterStateReader 
  "Water state reader which extends ThermoPower's StateReader_water"
  extends ThermoPower.PowerPlants.HRSG.Components.StateReader_water;
  
  Medium.Density rho "Temperature";
  Medium.SpecificEntropy s "Pressure";
  String fluid = Medium.mediumName "Fluid name";
    
equation 
  rho = Medium.density(fluidState);  
  s = Medium.specificEntropy(fluidState);
  
end WaterStateReader;
