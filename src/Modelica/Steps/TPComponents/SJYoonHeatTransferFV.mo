within Steps.TPComponents;

model SJYoonHeatTransferFV 
  "SJ Yoon heat transfer Correlation for zigzag PCHE with low Reynolds Number "
  extends ThermoPower.Thermal.BaseClasses.DistributedHeatTransferFV;
  
  // Turned out SJYoon Nusselt number correlation, although works for low Reynolds number case,
  // May not suitable for the Thermo Oil with high viscosity (PH HPL test case). 
  // So this correlation has been abondaned. 
  
  import ThermoPower.Thermal.BaseClasses;
  import ThermoPower.Thermal;
  import ThermoPower.Functions;
  parameter Real phi = 35 "pitch angle, alpha in SJYoon-2017, degree °, [5°,45°]";  
  parameter Real Re_min = 100 "Minimum Reynolds number";
  final parameter Real Pr_min = 0.7 "Minimum Prandtl number";
  final parameter Real Pr_max = 160 "Maximum Prandtl number";
  SI.CoefficientOfHeatTransfer gamma[Nf] "Heat transfer coefficients at the nodes";
  SI.CoefficientOfHeatTransfer gamma_vol[Nw] "Heat transfer coefficients at the volumes";
  Medium.Temperature Tvol[Nw] "Fluid temperature in the volumes";
  Medium.DynamicViscosity mu[Nf] "Dynamic viscosity";
  Medium.ThermalConductivity k[Nf] "Thermal conductivity";
  Medium.SpecificHeatCapacity cp[Nf] "Heat capacity at constant pressure";
  SI.PerUnit Re[Nf] "Reynolds number";
  SI.PerUnit Re_vol[Nw] "Reynolds number per volume, for comparison (with Other HT models) purpose only";
  SI.PerUnit Pr[Nf] "Prandtl numbers";
  SI.PerUnit Nu[Nf] "Nusselt numbers";
  SI.PerUnit Nu_vol[Nw] "Nusselt numbers per volume, for comparison (with Other HT models) purpose only";
  SI.PerUnit Re_l[Nf] "Reynolds number limited to validity range";
  SI.PerUnit Pr_l[Nf] "Prandtl number limited to validity range";
equation
  assert(Nw == Nf - 1, "Number of volumes Nw on wall side should be equal to number of volumes fluid side Nf - 1");
// Fluid properties at the nodes
  for j in 1:Nf loop
    mu[j] = Medium.dynamicViscosity(fluidState[j]);
    k[j] = Medium.thermalConductivity(fluidState[j]);
    cp[j] = Medium.heatCapacity_cp(fluidState[j]);
    Re[j] = abs(w[j] * Dhyd / (A * mu[j]));
    Pr[j] = cp[j] * mu[j] / k[j];
    Re_l[j] = Functions.smoothSat(Re[j], Re_min, 1e9, Re_min / 2);
    Pr_l[j] = Functions.smoothSat(Pr[j], Pr_min, Pr_max, Pr_min / 2, Pr_max / 10);
    // calculate the Nusselt Number explicitly    
    Nu[j] = 5.05 + (0.02 * phi + 0.003) * Re_l[j] * (Pr_l[j] ^ 0.6);

    gamma[j] = Nu[j] * k[j] /Dhyd;
  end for;
  for j in 1:Nw loop
    Tvol[j] = if useAverageTemperature then (T[j] + T[j + 1]) / 2 else T[j + 1];
    gamma_vol[j] = if useAverageTemperature then (gamma[j] + gamma[j + 1]) / 2 else gamma[j + 1];
    Nu_vol[j] = if useAverageTemperature then (Nu[j] + Nu[j + 1]) / 2 else Nu[j + 1];
    Re_vol[j] = if useAverageTemperature then (Re[j] + Re[j + 1]) / 2 else Re[j + 1];
    Qw[j] = (Tw[j] - Tvol[j]) * kc * omega * l * gamma_vol[j] * Nt;
  end for;
  annotation(
    Diagram(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}}), graphics),
    Icon(graphics = {Text(extent = {{-100, -52}, {100, -80}}, textString = "%name")}));
end SJYoonHeatTransferFV;
