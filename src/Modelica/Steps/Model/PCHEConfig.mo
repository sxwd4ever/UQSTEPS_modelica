within Steps.Model;

model PCHEConfig
  "Configuration for PCHE single component test."
  
  import Modelica.SIunits.Conversions.{from_degC, from_deg};
  import Modelica.SIunits.{Temperature, Pressure, SpecificEnthalpy};
  import Util = Utilities.Util;
  import Steps.Utilities.CoolProp.PropsSI; 
  import Steps.Components.PCHEGeoParam;  

  replaceable package medium_main   = Steps.Media.SCO2 constrainedby Modelica.Media.Interfaces.PartialPureSubstance;  
  replaceable package medium_heater = Steps.Media.MoltenSalt.MoltenSalt_pT constrainedby Modelica.Media.Interfaces.PartialPureSubstance;
  
  constant Real pi = Modelica.Constants.pi;
  
  // **** Boundary Conditions as Start values for recuperators - start ****
  // Core INPUT parameters for boundary conditions **** 
  // following values are mainly retrieved from Pinjarra Hill HPL the test case
  parameter Modelica.SIunits.Pressure p_main   = 12.567478e6;
  parameter Modelica.SIunits.Pressure p_heater = from_bar(100);
  
  parameter Modelica.SIunits.Temperature T_amb      = from_degC(35);
  parameter Modelica.SIunits.Temperature T_hot_in   = from_degC(103.222748);
  parameter Modelica.SIunits.Temperature T_hot_out  = from_degC(96.145935);
  parameter Modelica.SIunits.Temperature T_cold_in  = from_degC(28.910231);
  parameter Modelica.SIunits.Temperature T_cold_out = from_degC(99.666342);

  parameter Modelica.SIunits.MassFlowRate mdot_main   = 0.0854299999999999;
  parameter Modelica.SIunits.MassFlowRate mdot_heater = 1.218944;

  parameter Integer N_seg = 10 "default number of discretized segments in one tube";
  
  // Thermal-hydraulic properties
   // Stainless 316, 316L, 317, 317L
  parameter Modelica.SIunits.Density rho_wall             = 8030 "density of wall, kg/m3";
  parameter Modelica.SIunits.SpecificHeatCapacity cp_wall = 485 "cp of wall, J/kg-K";
  // Thermal conductivity of heater's metal wall 

  // **** START -  Definition of Key Point's state *****
  // points around PCHE
  parameter Model.ThermoState st_hot_in(
    p    = p_heater,
    T    = T_hot_in,
    h    = medium_heater.specificEnthalpy_pT(p = st_hot_in.p, T = st_hot_in.T),
    mdot = mdot_heater    
  ); 

  parameter Model.ThermoState st_hot_out(
    p    = p_heater,
    T    = T_hot_out,
    h    = medium_heater.specificEnthalpy_pT(p = st_hot_out.p, T = st_hot_out.T),
    mdot = mdot_heater    
  ); 

  parameter Model.ThermoState st_cold_in(
    p    = p_main,
    T    = T_cold_in,
    h    = medium_main.specificEnthalpy_pT(p = st_cold_in.p, T = st_cold_in.T),
    mdot = mdot_main    
  ); 

  parameter Model.ThermoState st_cold_out(
    p    = p_main,
    T    = T_cold_out,
    h    = medium_main.specificEnthalpy_pT(p = st_cold_out.p, T = st_cold_out.T),
    mdot = mdot_main    
  );  
  // **** END   -  Definition of Key Point's state *****


  // **** START - Definition of components' config *****
  // PCHE
  // PCHE's size of heat exchanger gas - gas
  // N_ch groups of fluid(cold, inner)-tube-gas(hot, outter) tubes 
  parameter Modelica.SIunits.Radius r_i     = 1.72e-3 / 2 "mm tube's internal radius";
  parameter Modelica.SIunits.Radius r_t     = 1.72e-3 / 2 "tube's external radius";
  parameter Modelica.SIunits.Radius r_o     = 1.72e-3 / 2 "radius of external side of one group";
  parameter Modelica.SIunits.Radius w_sd    = 2.3e-3 "Width of the solid domain";
  parameter Modelica.SIunits.Radius h_sd    = 4.17e-3 "Height of the solid domain, containing one cold tube and one hot tube";
  parameter Modelica.SIunits.Length t_ch    = 0.51e-3 "thickness between two channels";
  parameter Modelica.SIunits.Length l_pitch = 12.3e-3 "pitch length";
  parameter Modelica.SIunits.Length a_phi   = 35 "pitch angle";

  parameter Modelica.SIunits.Length L      = 250 * 1e-3 "m";
  parameter Modelica.SIunits.Length l_wall = 420e-3 "Length of wall, not necessarily equals to length of flow path";  
  parameter Integer N_ch                   = integer(1200);

  // cfg for heater's cold/fluid side
  parameter AreaGeometry ga_hot  = SetAreaGeometry_SemiCircle(r = r_i);
  parameter AreaGeometry ga_cold = SetAreaGeometry_SemiCircle(r = r_o);
  parameter AreaGeometry ga_wall = SetAreaGeometry_Wall(r = r_i, w = w_sd, h = h_sd, p1 = t_ch);
  parameter PathGeometry gp_hot  = SetPathGeometry(geo_area = ga_hot, L = L, N_seg = N_seg);
  parameter PathGeometry gp_cold = SetPathGeometry(geo_area = ga_cold, L = L, N_seg = N_seg);
  parameter PathGeometry gp_wall = SetPathGeometry(geo_area = ga_wall, L = L, N_seg = N_seg, p1 = l_wall);
  
  parameter Model.HeatExchangerConfig cfg_PCHE(
    cfg_hot(
      st_in    = st_hot_in,
      st_out   = st_hot_out,
      geo_area = ga_hot,
      geo_path = gp_hot,
      N_ch     = N_ch,
      u        = 0,
      l_pitch  = l_pitch,
      a_phi    = a_phi
    ),
    cfg_cold(
      st_in    = st_cold_in,
      st_out   = st_cold_out,
      geo_area = ga_cold,
      geo_path = gp_cold,
      N_ch     = N_ch,
      u        = 0,
      l_pitch  = l_pitch,
      a_phi    = a_phi
    ),
    cfg_fluid = cfg_PCHE.cfg_cold,
    cfg_gas   = cfg_PCHE.cfg_hot,
    cfg_wall(
      st_init  = st_hot_in,
      geo_area = ga_wall,
      geo_wall = gp_wall,
      // table_k  = table_k_wall,
      rho_mcm  = rho_wall * cp_wall,
      lambda   = 200
    )
  );  
  // **** END - Definition of Components' Config *****  

  
end PCHEConfig;
