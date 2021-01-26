within Steps.TPComponents;


model MarchionniPCHEHeatTransferFV "Marchionni [2019] heat transfer Correlation"
  extends BaseClasses.DistributedHeatTransferFV;
  import SI = Modelica.SIunits;
  import MyUtil = Steps.Utilities.Util;
  import ThermoPower.Thermal.BaseClasses;
  import ThermoPower.Thermal;
  import ThermoPower.Functions;

  //parameter Modelica.SIunits.Length pitch = 24.6 * 1e-3;
  //parameter Modelica.SIunits.Angle phi = 0.0 "unit rad";
  parameter Real Re_min = 2300 "Minimum Reynolds number";
  
  SI.Length d_c = (8 * A / Modelica.Constants.pi) ^ 0.5 "Diameter of semi_circular";
  /*
  inner Modelica.SIunits.Diameter d_h = 4 * A_c / peri_c "Hydraulic Diameter";  
  inner Modelica.SIunits.Area A_c = Modelica.Constants.pi * d_c * d_c / 8 "Area of semi-circular tube"; 
  inner Modelica.SIunits.Length peri_c = d_c * Modelica.Constants.pi / 2 + d_c "perimeter of semi-circular";   
  */
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
  Real co_A[Nf] "coefficients in Eq. 2 of [Marchionni 2019]";
  Real co_B[Nf] "coefficients in Eq. 3 of [Marchionni 2019]";  
  Real f[Nf] "Fanning Friction Factor - used to calculate pressure drop";
  Real dp[Nf];
  Real Pr[Nf];
  SI.ThermalConductivity k_wall[Nw] "wall thermal conductivity - determined by material of wall and local temperature";
  Modelica.SIunits.Length t_wall = (2 - Modelica.Constants.pi / 4) * (d_c / 2) "thickness of wall between two neighboring hot and cold";
  parameter SI.Length pitch "pitch length";
  parameter Real phi"pitch angle, degree";
  parameter String name_material = "inconel 750";
  parameter Real kc_dp = 2 "Pressure drop correction coeffecient";
  parameter Real kc_cf =if abs(phi - 0) < Modelica.Constants.eps then 1 else 12; 
  // DO NOT set values for these two parameters since use_rho_bar will be used as flags in if-statement. 
  parameter Real use_rho_bar "> 0, use rho_bar for dp calculation. error in passing a boolean from OMPython so Real type variable is used here";
  parameter Real rho_bar "Averaged rho, >0: valid rho and will be used for dp calculation";   
  
  Real C1 = kc_cf / kc_dp "calibration coefficient in Eq. 1 in [Marchionni 2019], = 1 / kc_dp for straight Channel and = 12 / kc_dp for zigzag Channel, which is conducted by numerical simulation against CFD result in [Meshram 2012]";
  Real C2 = 1/C1 "calibration coefficient in Eq. 4 in [Marchionni 2019], = 1.0 / C1 ,which is conducted by numerical simulation against CFD result in [Meshram 2012]";
  parameter Real table_th_conductivity[:, :] = [149, 16.9; 316, 20.5; 538, 26.5; 649, 28.7; 760, 31.4; 871, 35.3];
  
  Real kc_T[Nf];
  Modelica.Blocks.Tables.CombiTable1D th_conductivity(tableOnFile = false, table = table_th_conductivity, tableName = "conductivity", fileName = "conductivity", smoothness = Modelica.Blocks.Types.Smoothness.ContinuousDerivative) annotation(
    Placement(transformation(extent = {{60, -80}, {80, -60}}, rotation = 0)));    
   
equation
  assert(Nw == Nf - 1, "Number of volumes Nw on wall side should be equal to number of volumes fluid side Nf - 1");
   
// Fluid properties at the nodes    
  for j in 1:Nf loop
    mu[j] = noEvent(Medium.dynamicViscosity(fluidState[j]));
    k[j] = noEvent(Medium.thermalConductivity(fluidState[j]));
    rho[j] = noEvent(Medium.density(fluidState[j]));
    
    kc_T[j] = if Medium.temperature(fluidState[j]) < 600 then 2.0 else 1.0;
    
    G[j] = noEvent(abs(w[j] / A));
    Re[j] =  noEvent(G[j] * Dhyd / mu[j]);
    Re_l[j] = noEvent(Functions.smoothSat(Re[j], Re_min, 1e9, Re_min / 2));
    u[j] = noEvent(abs(w[j]) / A / rho[j]);
    // pressure drop calculation in Marchionni 2019
    co_A[j] = -2.0 * Modelica.Math.log10( 12 / Re_l[j]) "neglect roughness";
    co_B[j] = -2.0 * Modelica.Math.log10( 2.51 * co_A[j] / Re_l[j]);  

    f[j] = C1 * (0.25 * (( 4.781 - (co_A[j] - 4.781) ^ 2 / (co_B[j] - 2 * co_A[j] + 4.781)) ^(-2)));
    // dp[j] = noEvent(kc_T[j] * kc_dp * f[j] * l * rho[j] * u[j] ^ 2 / Dhyd);
    // error in passing a boolean parameter by OMPython, so I use negetive value to indicate a invalid rho, which will not be used for dp
    //dp[j] = noEvent(kc_dp * f[j] * l * rho_bar * u[j] ^ 2 / Dhyd);
    
    if use_rho_bar > 0 then    
      dp[j] = noEvent(kc_dp * f[j] * l * rho_bar * u[j] ^ 2 / Dhyd);
    else
      dp[j] = noEvent(kc_dp * f[j] * l * rho[j] * u[j] ^ 2 / Dhyd);
    end if;
    
    Pr[j] = noEvent(Medium.specificHeatCapacityCp(fluidState[j]) * mu[j] / k[j]);    
    Nu[j] = noEvent(C2 * ((f[j]/2 * (Re_l[j] - 1000) * Pr[j]) / (1 + 12.7 * ( Pr[j] ^(2/3) - 1) * (f[j]/2)^0.5)));

    hc[j] = noEvent(Nu[j] * k[j] / Dhyd);    
    
  end for;
  th_conductivity.u[1] = 0.0;
  for j in 1:Nw loop
    Tvol[j] = if useAverageTemperature then (T[j] + T[j + 1]) / 2 else T[j + 1];
    //th_conductivity.u[1] = (Tw[j] + Tvol[j]) / 2;
    //k_wall[j] =  MyUtil.metal_conductivity(th_conductivity.tableID, (Tw[j] + Tvol[j]) / 2);
    k_wall[j] =  MyUtil.metal_conductivity(th_conductivity.tableID, Tw[j]);
    gamma[j] = 1 / (1 / hc[j] + t_wall / k_wall[j]);
    Qw[j] = (Tw[j] - Tvol[j]) * kc * omega * l * gamma[j] * Nt;
  end for;
end MarchionniPCHEHeatTransferFV;
