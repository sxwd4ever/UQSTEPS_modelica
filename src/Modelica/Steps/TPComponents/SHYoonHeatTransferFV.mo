within Steps.TPComponents;

model SHYoonHeatTransferFV "Marchionni [2019] heat transfer Correlation"
  extends BaseClasses.DistributedHeatTransferFV;
  import SI = Modelica.SIunits;
  import MyUtil = Steps.Utilities.Util;
  import ThermoPower.Thermal.BaseClasses;
  import ThermoPower.Thermal;
  import ThermoPower.Functions;
  import Modelica.Math.log;

  //parameter Modelica.SIunits.Length pitch = 24.6 * 1e-3;
  //parameter Modelica.SIunits.Angle phi = 0.0 "unit rad";
  parameter Real Re_min = 2300 "Minimum Reynolds number";
  
  SI.Length d_c = (8 * A / Modelica.Constants.pi) ^ 0.5 "Diameter of semi_circular";
  /*
  inner Modelica.SIunits.Diameter d_h = 4 * A_c / peri_c "Hydraulic Diameter";  
  inner Modelica.SIunits.Area A_c = Modelica.Constants.pi * d_c * d_c / 8 "Area of semi-circular tube"; 
  inner Modelica.SIunits.Length peri_c = d_c * Modelica.Constants.pi / 2 + d_c "perimeter of semi-circular";   
  */
  
  Real g = Modelica.Constants.g_n;
  
  // volume properties
  // averaged with upstream/downstream node property
  Medium.Temperature T_vol[Nw] "Fluid temperature in the volumes";  
  Real G_vol[Nw] "mass flow flux";  
  SI.Velocity u_vol[Nw] "Local velocity of fluid";
  Medium.DynamicViscosity mu_vol[Nw] "Dynamic viscosity";  
  Medium.ThermalConductivity k_vol[Nw] "Thermal conductivity";
  Medium.Density rho_vol[Nw] "Averaged Density of fluid";  
  SI.SpecificHeatCapacityAtConstantPressure cp_vol[Nw] "HeatCapacityAtConstantPressure";  
  
  // calculated by local averaged properties
  parameter Boolean heating = true "true if fluid is heated, false if fluid is cooled";  
  Real Pr[Nw];
  Real Pr_cmp[Nw] "Dummy variables, for align the comparision";
  Real f[Nw] "Fanning Friction Factor - used to calculate pressure drop";
  Real dp[Nw];
  SI.PerUnit hc[Nw] "Local heat transfer coefficient";
  SI.PerUnit Re[Nw] "Reynolds number";
  SI.PerUnit Nu[Nw] "Nussult numbers";
  SI.PerUnit Re_l[Nw] "Reynolds number limited to validity range";
  SI.PerUnit Re_cmp[Nw] "dummpy Reynolds number per volume, for comparison (with Other HT models) purpose only";  
  SI.CoefficientOfHeatTransfer gamma[Nw] "Heat transfer coefficients at the volumes";  
  parameter Real gamma_max = 15000, gamma_min = 1500 "max and min valid gamma, to prevent unfeasible gamma";  

  // Real kc_T[Nw];    
  Medium.AbsolutePressure p_vol[Nw] "volume averaged pressure";
  SI.SpecificEnthalpy h_vol[Nw];
  
  // quantities evaluated at corresponding wall volume temperature
  Medium.ThermodynamicState state_w[Nw] "ThermoDynamicState";
  SI.Density rho_w[Nw] "density";
  SI.SpecificEnthalpy h_w[Nw];
  SI.SpecificHeatCapacityAtConstantPressure cp_bar[Nw] "mean specific heat in Eq. 12 in Liao's paper";
  Real Gr[Nw] "Grasholf number";
  
  Real q1[Nw], q2[Nw], q3[Nw] "intermediate values for debug/correlation modelling purpose";
  
  // node properties    
  SI.Velocity u[Nf] "Local velocity of fluid";
  Real G[Nf] "mass flow flux";
  Medium.DynamicViscosity mu[Nf] "Dynamic viscosity";
  Medium.ThermalConductivity k[Nf] "Thermal conductivity";
  Medium.Density rho[Nf] "Density of fluid";  
  SI.SpecificHeatCapacityAtConstantPressure cp[Nf] "HeatCapacityAtConstantPressure";  
  SI.SpecificEnthalpy h[Nf] ;

  Modelica.SIunits.Length t_wall = (2 - Modelica.Constants.pi / 4) * (d_c / 2) "thickness of wall between two neighboring hot and cold";
  parameter SI.Length pitch "pitch length";
  parameter Real phi "pitch angle, degree";
  parameter String name_material = "inconel 750";
  // parameter Real kc_dp = 2 "Pressure drop correction coeffecient";
  // parameter Real kc_cf =if abs(phi - 0) < Modelica.Constants.eps then 1 else 12; 
  // DO NOT set values for these two parameters since use_rho_bar will be used as flags in if-statement. 
  parameter Real use_rho_bar "> 0, use rho_bar for dp calculation. error in passing a boolean from OMPython so Real type variable is used here";
  parameter Real rho_bar "Averaged rho, >0: valid rho and will be used for dp calculation";   

  parameter Real Cf_C1 = 1.0;
  parameter Real Cf_C2 = 1.0;
  parameter Real Cf_C3 = 1.0;
  
  // Real C1[Nw] "calibration coefficient in Eq. 1 in [Marchionni 2019], = 1 / kc_dp for straight Channel and = 12 / kc_dp for zigzag Channel, which is conducted by numerical simulation against CFD result in [Meshram 2012]";
  //Real C2[Nw] "calibration coefficient in Eq. 4 in [Marchionni 2019], = 1.0 / C1 ,which is conducted by numerical simulation against CFD result in [Meshram 2012]";
  parameter Real table_th_conductivity[:, :] = [149, 16.9; 316, 20.5; 538, 26.5; 649, 28.7; 760, 31.4; 871, 35.3];
  Modelica.Blocks.Tables.CombiTable1D th_conductivity(tableOnFile = false, table = table_th_conductivity, tableName = "conductivity", fileName = "conductivity", smoothness = Modelica.Blocks.Types.Smoothness.ContinuousDerivative) annotation(
    Placement(transformation(extent = {{60, -80}, {80, -60}}, rotation = 0)));    

  SI.ThermalConductivity k_wall[Nw] "wall thermal conductivity - determined by material of wall and local temperature";

  
equation
  assert(Nw == Nf - 1, "Number of volumes Nw on wall side should be equal to number of volumes fluid side Nf - 1");
   
// Fluid properties at the nodes    
  for j in 1:Nf loop
    G[j] = noEvent(abs(w[j])/ A); 
    u[j] = noEvent(abs(w[j]) / A / rho[j]);
    mu[j] = noEvent(Medium.dynamicViscosity(fluidState[j]));
    k[j] = noEvent(Medium.thermalConductivity(fluidState[j]));
    rho[j] = noEvent(Medium.density(fluidState[j]));
    cp[j] = noEvent(Medium.specificHeatCapacityCp(fluidState[j]));   
    h[j] = noEvent(Medium.specificEnthalpy(fluidState[j])); 
  end for;   
  
  // required for table initialization and equation balance
  th_conductivity.u[1] = 0.0;
  
  for j in 1:Nw loop  
    // calculate averaged values    
    if useAverageTemperature then
      T_vol[j] = noEvent((T[j] + T[j + 1]) / 2);
      G_vol[j] = noEvent((G[j] + G[j+1]) / 2);
      u_vol[j] = noEvent((u[j] + u[j+1]) / 2);
      mu_vol[j] = noEvent((mu[j] + mu[j+1]) / 2);
      k_vol[j] = noEvent((k[j] + k[j+1]) / 2);
      rho_vol[j] = noEvent((rho[j] + rho[j+1]) / 2);
      cp_vol[j] = noEvent((cp[j] + cp[j+1]) / 2); 
      p_vol[j] = noEvent((fluidState[j].p + fluidState[j+1].p) / 2);
      h_vol[j] = noEvent((h[j] + h[j+1]) / 2);
    else
      // upstream
      T_vol[j] = noEvent(T[j]);
      G_vol[j] = noEvent(G[j]);
      u_vol[j] = noEvent(u[j]);
      mu_vol[j] = noEvent(mu[j]);
      k_vol[j] = noEvent(k[j]);
      rho_vol[j] = noEvent(rho[j]);
      cp_vol[j] = noEvent(cp[j]);
      p_vol[j] = noEvent(fluidState[j].p);
      h_vol[j] = noEvent(h[j]);
    end if;    
    
    // kc_T[j] = if T_vol[j] < 600 then 2.0 else 1.0;  
    
    // calculate thermaldynamic state and relative quantities on wall temperature       
    // Only for Liao's correlation and its update since the calculation is time consuming
    state_w[j] = Medium.setState_pT(p_vol[j], Tw[j]);
    rho_w[j] = Medium.density(state_w[j]);
    h_w[j] = Medium.specificEnthalpy(state_w[j]);    
    // use abs() to prevent negetive Gr
    Gr[j] = abs((rho_w[j] - rho_vol[j]) * rho_vol[j] * g * (d_c ^ 3) / (mu_vol[j]^2));   
    
    if abs(Tw[j] - T_vol[j]) < 1e-3 then
      cp_bar[j] = cp_vol[j]; 
    else
      cp_bar[j] = (h_w[j] - h_vol[j]) / (Tw[j] - T_vol[j]);  
    end if;
    
    // calculate Re and Pr  
    Re[j] =  noEvent(G_vol[j] * Dhyd / mu_vol[j]);    
    Pr[j] = noEvent(cp_vol[j] * mu_vol[j] / k_vol[j]);      
          
    // low Reynolds number case
    if Re[j] < Re_min then 
      // for low Reynolds case, use Dittus-Boelter correlations
      // Nu[j] = noEvent(0.023 * Re[j] ^ 0.8 * Pr[j] ^ 0.4);   
      // Yoon-Sun's Correlations for (200 < Re < 4700)
      // Nu[j] = 5.05 + (0.02 * from_deg(phi) + 0.003 ) * Re[j] * (Pr[j] ^ 0.6);
      
      // Kim's Correlations for (0 < Re < 2500)
      if abs(w[j]) < 1e-6 then
        Nu[j] = 1e-10; // ,1 extrem case, no mass flow or heat transfer.
      else
        Nu[j] = 4.089 + 0.00365 * Re[j] * (Pr[j] ^ 0.58);
      end if;
        
      // dummy variables  
      Re_l[j] = Re_min; 
   
      f[j] = 0;
      
      q1[j] = 0;
      q2[j] = 0;
      q3[j] = 0; 
      dp[j] = 0;  
    else      
      Re_l[j] = noEvent(Functions.smoothSat(Re[j], Re_min, 1e9, Re_min / 2));               
  
      q1[j] = noEvent(Gr[j] / (Re_l[j] ^ 2));      
      q2[j] = noEvent(rho_w[j] / rho_vol[j]);      
      q3[j] = noEvent(cp_bar[j] / cp_vol[j]);  

      // Nusselt number by Liao-Zhao correlation
      Nu[j] = 0.14 * (Re_l[j] ^ 0.69) * (Pr[j] * 0.66); // for T_b > T_pc which is true for sco2 RCBC
      // Nu[j] = 0.013 * (Re_l[j] ^ 1.0) * (Pr[j] * -0.05) * ( rho_pc / rho[j]); // for T_b < T_pc which is true for sco2 RCBC
      // Nu[j] = noEvent(0.124 * (Re_l[j] ^ 0.8) * (Pr[j] ^ 0.4) * q1[j] ^ 0.203 * q2[j] ^ 0.842 * q3[j] ^ 0.384);   

      // No f[j] in Liao's correlation, so use Ngo's correlation as default
      f[j] = noEvent(0.1924 * (Re_l[j] ^ (-0.091)));    
      
      if use_rho_bar > 0 then    
        // dp[j] = noEvent(kc_dp * f[j] * l * rho_bar * u_vol[j] ^ 2 / Dhyd);
        dp[j] = noEvent(f[j] * l * rho_bar * u_vol[j] ^ 2 / Dhyd);
      else
        // dp[j] = noEvent(kc_dp * f[j] * l * rho[j] * u_vol[j] ^ 2 / Dhyd);
        dp[j] = noEvent(f[j] * l * rho[j] * u_vol[j] ^ 2 / Dhyd);
      end if;  
      //   
    end if;
    Re_cmp[j] = Re_l[j];  
    Pr_cmp[j] = Pr[j];     
    hc[j] = noEvent(Nu[j] * k_vol[j] / Dhyd);   

    //th_conductivity.u[1] = (Tw[j] + T_vol[j]) / 2;
    // k_wall[j] =  MyUtil.metal_conductivity(th_conductivity.tableID, (Tw[j] + T_vol[j]) / 2);
    k_wall[j] =  MyUtil.metal_conductivity(th_conductivity.tableID, Tw[j]);
    
    gamma[j] = noEvent(1 / (1 / hc[j] + t_wall / k_wall[j]));
    //gamma[j] = noEvent(hc[j]);
    
    /*
    if abs(hc[j]) < 1e-6 then
      gamma[j] = 0;
    else
      gamma[j] = noEvent(1 / (1 / hc[j] + t_wall / k_wall[j]));
    end if;
    
    
    MyUtil.myAssert(
    debug = false, 
    val_test = Tw[j], min = 0, max = 1e6,
    name_val = "Tw[j]", 
    val_ref = {kc, gamma[j], hc[j], k_vol[j]}, 
    name_val_ref = {"kc", "gamma[j]", "hc[j]", "k_vol[j]"});  
    
    MyUtil.myAssert(
    debug = false, 
    val_test = T_vol[j], min = 0, max = 1e6,
    name_val = "T_vol[j]", 
    val_ref = {kc, gamma[j], hc[j], k_vol[j]}, 
    name_val_ref = {"kc", "gamma[j]", "hc[j]", "k_vol[j]"});        
    */
    
    Qw[j] = noEvent((Tw[j] - T_vol[j]) * kc * omega * l * gamma[j] * Nt);
  end for;
/*        
algorithm 
    // ******************************************** 
    // this algorithm section is used for debug purpose
    // if a variable become invalid (zero, inf or weired)
    // 1. make a copy of its equation
    // 2. move the copy into this section, rewrite the equation into an algorithm ('=' -> ':=')
    // 3. comment the origin one
    // 4. update the calling of function MyAseert accordingly to print value in console
    //
    // once the bug fixed, DO REMEMBER to rewind the changed lines reversely. 
    // ********************************************
    for j in 1:Nw loop
      Qw[j] := (Tw[j] - T_vol[j]) * kc * omega * l * gamma[j] * Nt;    
      
        
      MyUtil.myAssert(
      debug = false, 
      val_test = Tw[j], min = 0, max = 1e6,
      name_val = "Tw[j]", 
      val_ref = {kc, gamma[j], hc[j], k_wall[j], k_vol[j]}, 
      name_val_ref = {"kc", "gamma[j]", "hc[j]", "k_wall[j]", "k_vol[j]"});    
      
    
      MyUtil.myAssert(
      debug = false, 
      val_test = T_vol[j], min = 0, max = 1e6,
      name_val = "T_vol[j]", 
      val_ref = {kc, gamma[j], hc[j], k_wall[j], k_vol[j]}, 
      name_val_ref = {"kc", "gamma[j]", "hc[j]", "k_wall[j]", "k_vol[j]"});    
    end for;
*/
end SHYoonHeatTransferFV;
