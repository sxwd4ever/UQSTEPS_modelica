within Steps.Components;


model KimPCHEHeatTransferFV "Kim [2012] heat transfer Correlation"
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
  Real f[Nf] "Fanning Friction Factor - used to calculate pressure drop";
  SI.PressureDifference dp[Nf];
  SI.ThermalConductivity k_wall[Nw] "wall thermal conductivity - determined by material of wall and local temperature";
  Modelica.SIunits.Length t_wall = (2 - Modelica.Constants.pi / 4) * (d_c / 2) "thickness of wall between two neighboring hot and cold";
  parameter SI.Length pitch = 12.3e-3 "pitch length";
  parameter Real phi = 35 "pitch angle °";
  parameter String name_material = "inconel 750";

  // default value for d_c = 2mm(Dhyd=0.922mm) in Kim [2012]
  parameter Real table_kim_cor_a[:, :] = [0, 12.3e-3, 24.6e-3; 5, 0.00366, 0.0019; 10, 0.02536, 0.01775; 15, 0.0696 , 0.06455; 20, 0.12817, 0.08918; 25, 0.19392, 0.1442; 30, 0.27701 , 0.21115; 35, 0.37934, 0.29342 ; 40, 0.49995 , 0.39158; 45, 0.65898, 0.51539]  "Table for kim_correlation_a";

  parameter Real table_kim_cor_b[:, :] = [0, 12.3e-3, 24.6e-3; 5, 1.07165, 1.1187; 10, 0.91728, 0.90795 ; 15, 0.85362, 0.81021 ; 20, 0.83085, 0.8136 ; 25, 0.82838, 0.79529 ; 30, 0.82549, 0.78615 ; 35, 0.82413, 0.78118 ; 40, 0.82555, 0.7799 ; 45, 0.82508, 0.77944]  "Table for kim_correlation_a";

  parameter Real table_kim_cor_c[:, :] = [0, 12.3e-3, 24.6e-3; 5, 0.00061,   0.00036; 10, 0.00251 ,  0.0022; 15, 0.00607,  0.00544; 20, 0.0118,  0.00894; 25, 0.02165,  0.01321; 30, 0.03006,  0.01877; 35, 0.03845,  0.02578; 40, 0.04829,  0.03489; 45, 0.06368,  0.04885]  "Table for kim_correlation_a";

  parameter Real table_kim_cor_d[:, :] = [0, 12.3e-3, 24.6e-3; 5, 1.13917,   1.18213; 10, 1.01844,  0.99841; 15, 0.93106,  0.91361; 20, 0.85964,  0.86708; 25, 0.79051,  0.82738; 30, 0.76007,  0.78945; 35, 0.73793,  0.75403; 40, 0.7159,  0.71951; 45, 0.68516,  0.67896]  "Table for kim_correlation_a";
  
  Modelica.Blocks.Tables.CombiTable2D kim_cor_a(tableOnFile = false, table = table_kim_cor_a, tableName = "kim_cor_a", fileName = "kim_cor_a", smoothness = Modelica.Blocks.Types.Smoothness.ContinuousDerivative) annotation(
    Placement(transformation(extent = {{60, -80}, {80, -60}}, rotation = 0)));
    
  Modelica.Blocks.Tables.CombiTable2D kim_cor_b(tableOnFile = false, table = table_kim_cor_b, tableName = "kim_cor_b", fileName = "kim_cor_b", smoothness = Modelica.Blocks.Types.Smoothness.ContinuousDerivative) annotation(
    Placement(transformation(extent = {{60, -80}, {80, -60}}, rotation = 0)));
    
  Modelica.Blocks.Tables.CombiTable2D kim_cor_c(tableOnFile = false, table = table_kim_cor_c, tableName = "kim_cor_c", fileName = "kim_cor_c", smoothness = Modelica.Blocks.Types.Smoothness.ContinuousDerivative) annotation(
    Placement(transformation(extent = {{60, -80}, {80, -60}}, rotation = 0)));
    
  Modelica.Blocks.Tables.CombiTable2D kim_cor_d(tableOnFile = false, table = table_kim_cor_d, tableName = "kim_cor_d", fileName = "kim_cor_d", smoothness = Modelica.Blocks.Types.Smoothness.ContinuousDerivative) annotation(
    Placement(transformation(extent = {{60, -80}, {80, -60}}, rotation = 0)));            
  
  parameter Real table_th_conductivity[:, :] = [149, 16.9; 316, 20.5; 538, 26.5; 649, 28.7; 760, 31.4; 871, 35.3];
  
  Modelica.Blocks.Tables.CombiTable1D th_conductivity(tableOnFile = false, table = table_th_conductivity, tableName = "conductivity", fileName = "conductivity", smoothness = Modelica.Blocks.Types.Smoothness.ContinuousDerivative) annotation(
    Placement(transformation(extent = {{60, -80}, {80, -60}}, rotation = 0)));    
        
  //Real kim_cor[4] = {kim.a, kim.b, kim.c, kim.d};
equation
  assert(Nw == Nf - 1, "Number of volumes Nw on wall side should be equal to number of volumes fluid side Nf - 1");
  kim_cor_a.u1 = phi;
  kim_cor_a.u2 = pitch;  
  kim_cor_b.u1 = phi;
  kim_cor_b.u2 = pitch;  
  kim_cor_c.u1 = phi;
  kim_cor_c.u2 = pitch;
  kim_cor_d.u1 = phi;
  kim_cor_d.u2 = pitch;
    
// Fluid properties at the nodes    
  for j in 1:Nf loop
    mu[j] = noEvent(Medium.dynamicViscosity(fluidState[j]));
    k[j] = noEvent(Medium.thermalConductivity(fluidState[j]));
    rho[j] = noEvent(Medium.density(fluidState[j]));
    G[j] = noEvent(abs(w[j] / A));
    Re[j] =  noEvent(G[j] * Dhyd / mu[j]);
    Re_l[j] = noEvent(Functions.smoothSat(Re[j], Re_min, 1e9, Re_min / 2));
    Nu[j] = noEvent(4.089 + kim_cor_c.y * Re_l[j] ^ kim_cor_d.y);
    //Nu[j] = noEvent(4.089 + 0.04247 * Re_l[j] ^ 0.70055);
    u[j] = noEvent(abs(w[j]) / A / rho[j]);
    hc[j] = noEvent(Nu[j] * k[j] / Dhyd);
    f[j] = noEvent((15.78 + kim_cor_a.y * Re_l[j] ^ kim_cor_b.y) / Re_l[j]);
    //f[j] = noEvent((15.78 + 0.35159 * Re_l[j] ^ 0.78015) / Re_l[j]);
//pressure drop, unit Pa
    dp[j] = noEvent(2 * f[j] * l * rho[j] * u[j] ^ 2 / Dhyd);
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
end KimPCHEHeatTransferFV;
