within Steps.Test;

model TestThermiaOilD "Test ThermiaOilD Medium model"
  extends Modelica.Icons.Example;
  package Medium = Media.ThermiaOilD "Medium model (ThermiaOilD)";
  Medium.BaseProperties medium;

  Medium.DynamicViscosity eta=Medium.dynamicViscosity(medium.state);
  Medium.ThermalConductivity lambda=Medium.thermalConductivity(medium.state);
  Medium.SpecificEntropy s=Medium.specificEntropy(medium.state);
  Medium.SpecificHeatCapacity cv=Medium.specificHeatCapacityCv(medium.state);
  Medium.SpecificInternalEnergy u=Medium.specificInternalEnergy(medium.state);
  Medium.SpecificInternalEnergy h=Medium.specificEnthalpy(medium.state);
  Medium.SpecificInternalEnergy d=Medium.density(medium.state);
  protected
  constant Modelica.SIunits.Time timeUnit = 1;
  constant Modelica.SIunits.Temperature Ta = 1;
equation
  medium.p = 1.013e5;
  medium.T = Medium.T_min + time/timeUnit*Ta;
    annotation (experiment(StopTime=1.01));
end TestThermiaOilD;
