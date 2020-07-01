within Steps.Interfaces;

connector PBFluidPort_b
  extends Modelica.Fluid.Interfaces.FluidPort_b;
  import SI = Modelica.SIunits;
  SI.ThermodynamicTemperature T "absolute temperature of the port";
equation

end PBFluidPort_b;