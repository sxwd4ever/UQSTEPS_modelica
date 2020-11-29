within Steps.Components;


model KimPCHEHeatTransferFV "Kim [2012] heat transfer Correlation"
  extends BaseClasses.DistributedHeatTransferFV(redeclare replaceable package Medium = ExternalMedia.Media.BaseClasses.ExternalTwoPhaseMedium);
  import SI = Modelica.SIunits;
  import MyUtil = Steps.Utilities.Util;
  import ThermoPower.Thermal.BaseClasses;
  import ThermoPower.Thermal;
  import ThermoPower.Functions;
  import Steps.Components.KimCorrelations;
  import Steps.Components.MaterialConductivity;
  parameter Modelica.SIunits.Length pitch = 24.6 * 1e-3;
  parameter Modelica.SIunits.Angle phi = 0.0 "unit rad";
  parameter Real Re_min = 10000 "Minimum Reynolds number";
  parameter String name_material = "inconel 750";
  SI.Length d_c = (8 * A / Modelica.Constants.pi) ^ 0.5 "Diameter of semi_circular";
  /*
                                                      inner Modelica.SIunits.Diameter d_h = 4 * A_c / peri_c "Hydraulic Diameter";  
                                                      inner Modelica.SIunits.Area A_c = Modelica.Constants.pi * d_c * d_c / 8 "Area of semi-circular tube"; 
                                                      inner Modelica.SIunits.Length peri_c = d_c * Modelica.Constants.pi / 2 + d_c "perimeter of semi-circular";   
                                                      */
  KimCorrelations kim_cor(phi = phi, pitch = pitch, d_h = Dhyd);
  MaterialConductivity mc(name_material = name_material);
  SI.CoefficientOfHeatTransfer gamma[Nw] "Heat transfer coefficients at the volumes";
  Medium.Temperature Tvol[Nw] "Fluid temperature in the volumes";
  Medium.DynamicViscosity mu[Nf] "Dynamic viscosity";
  Medium.ThermalConductivity k[Nf] "Thermal conductivity";
  Medium.Density rho[Nf] "Density of fluid";
  SI.Velocity u[Nf] "Local velocity of fluid";
  SI.PerUnit hc[Nf] "Local heat transfer coefficient";
  SI.PerUnit Re[Nf] "Reynolds number";
  SI.PerUnit Nu[Nf] "Nussult numbers";
  SI.PerUnit Re_l[Nf] "Reynolds number limited to validity range";
  Real G[Nf] "mass flow flux";
  Real f[Nf] "Fanning Friction Factor - used to calculate pressure drop";
  SI.PressureDifference dp[Nf];
  SI.ThermalConductivity k_wall[Nw] "wall thermal conductivity - determined by material of wall and local temperature";
  Modelica.SIunits.Length t_wall = (2 - Modelica.Constants.pi / 4) * (d_c / 2) "thickness of wall between two neighboring hot and cold";
equation
  assert(Nw == Nf - 1, "Number of volumes Nw on wall side should be equal to number of volumes fluid side Nf - 1");
// Fluid properties at the nodes
  for j in 1:Nf loop
    mu[j] = noEvent(Medium.dynamicViscosity(fluidState[j]));
    k[j] = noEvent(Medium.thermalConductivity(fluidState[j]));
    rho[j] = noEvent(Medium.density(fluidState[j]));
    G[j] = noEvent(abs(w[j] / A));
    Re[j] =  noEvent(G[j]/ A * Dhyd / mu[j]);
    Re_l[j] = noEvent(Functions.smoothSat(Re[j], Re_min, 1e9, Re_min / 2));
    Nu[j] = noEvent(4.089 + kim_cor.c * Re_l[j] ^ kim_cor.d);
    u[j] = noEvent(abs(w[j]) / A / rho[j]);
    hc[j] = noEvent(Nu[j] * k[j] / Dhyd);
    f[j] = noEvent((15.78 + kim_cor.a * Re_l[j] ^ kim_cor.b) / Re_l[j]);
//pressure drop, unit Pa
    dp[j] = noEvent(2 * f[j] * l * rho[j] * u[j] ^ 2 / Dhyd);
  end for;
  for j in 1:Nw loop
    Tvol[j] = if useAverageTemperature then (T[j] + T[j + 1]) / 2 else T[j + 1];
    k_wall[j] = MyUtil.thermal_conductivity(tableID = mc.table_th_inconel_750, name = name_material, temperature = (Tw[j] + Tvol[j]) / 2);
    gamma[j] = 1 / (1 / hc[j] + t_wall / k_wall[j]);
    Qw[j] = (Tw[j] - Tvol[j]) * kc * omega * l * gamma[j] * Nt;
  end for;
end KimPCHEHeatTransferFV;
