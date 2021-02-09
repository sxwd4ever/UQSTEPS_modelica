within Steps.Test;

model TestThermiaOilD "Test ThermiaOilD Medium model"
  extends Modelica.Icons.Example;
  package Medium = Media.ThermiaOilD "Medium model (ThermiaOilD)";
  Medium.ThermodynamicState state;

  Medium.DynamicViscosity eta=Medium.dynamicViscosity(state);
  Medium.ThermalConductivity lambda=Medium.thermalConductivity(state);
  Medium.SpecificEntropy s=Medium.specificEntropy(state);
  Medium.SpecificHeatCapacity cv=Medium.specificHeatCapacityCv(state);
  Medium.SpecificInternalEnergy u=Medium.specificInternalEnergy(state);
  Medium.SpecificInternalEnergy h=Medium.specificEnthalpy(state);
  Medium.SpecificInternalEnergy d=Medium.density(state);
  protected
  constant Modelica.SIunits.Time timeUnit = 1;
  constant Modelica.SIunits.Temperature Ta = 400;
equation
  state = Medium.setState_pT(p=1.013e5, T=Medium.T_min + time/timeUnit*Ta);
  // medium.p = 1.013e5;
  // medium.T = ;
    annotation (experiment(StopTime=1.01));
end TestThermiaOilD;
