within Steps.Model;

model RCBCycleConfig
  "full set of params set for Off-design RCBC - DON'T Alter it directly, use override instead, see Model PBConfigs"
  
  import Modelica.SIunits.Conversions.{from_degC, from_deg};
  import Modelica.SIunits.{Temperature, Pressure, SpecificEnthalpy};
  import Util = Utilities.Util;
  import Steps.Utilities.CoolProp.PropsSI; 
  import Steps.Components.PCHEGeoParam;  

  replaceable package medium_main   = Steps.Media.SCO2 constrainedby Modelica.Media.Interfaces.PartialPureSubstance;
  replaceable package medium_cooler = ThermoPower.Water.StandardWater constrainedby Modelica.Media.Interfaces.PartialPureSubstance;
  replaceable package medium_heater = Steps.Media.MoltenSalt.MoltenSalt_pT constrainedby Modelica.Media.Interfaces.PartialPureSubstance;
  
  constant Real pi = Modelica.Constants.pi;
  // efficiency of main compressor, bypass_compressor and turbine
  parameter Real eta_comp   = 0.89;
  parameter Real eta_recomp = 0.89;
  parameter Real eta_turb   = 0.89;
  parameter Real Ns_comp    = 2100;
  parameter Real Ns_recomp  = Ns_comp;  
  parameter Real Ns_turb    = Ns_comp;
  
  // **** Boundary Conditions as Start values for recuperators - start ****
  // Core INPUT parameters for boundary conditions **** 
  // following values are calculated by sscar for 10 MW power block( mdot_main = 125, T_amb = 35oC, split_ratio = 0.675)
  parameter Modelica.SIunits.Pressure p_comp_in  = 9e6;
  parameter Modelica.SIunits.Pressure p_comp_out = 20e6;
  parameter Modelica.SIunits.Pressure p_amb      = 1.01325e6;
  parameter Modelica.SIunits.Pressure p_heater   = 20e6;
  
  parameter Modelica.SIunits.Temperature T_amb             = from_degC(35);
  parameter Modelica.SIunits.Temperature T_HTR_hot_in      = from_degC(636.95057734);
  parameter Modelica.SIunits.Temperature T_HTR_cold_out    = from_degC(605.011512655);
  parameter Modelica.SIunits.Temperature T_HTR_hot_out     = from_degC(167.14450349);
  parameter Modelica.SIunits.Temperature T_HTR_cold_in     = from_degC(162.14458875);
  parameter Modelica.SIunits.Temperature T_LTR_cold_in     = from_degC(85.84344872);
  parameter Modelica.SIunits.Temperature T_LTR_hot_out     = from_degC(90.84344686);
  parameter Modelica.SIunits.Temperature T_heater_hot_in   = from_degC(800);
  parameter Modelica.SIunits.Temperature T_heater_hot_out  = from_degC(600);
  parameter Modelica.SIunits.Temperature T_heater_cold_out = from_degC(710);
  parameter Modelica.SIunits.Temperature T_cooler_cold_out = from_degC(61.7479462700001);
  parameter Modelica.SIunits.Temperature T_cooler_hot_out  = from_degC(45);
  parameter Modelica.SIunits.Temperature T_recom_out       = from_degC(227.7022889);

  parameter Modelica.SIunits.MassFlowRate mdot_main   = 128.774;
  parameter Modelica.SIunits.MassFlowRate mdot_comp   = 84.375;
  parameter Modelica.SIunits.MassFlowRate mdot_heater = 40;
  parameter Modelica.SIunits.MassFlowRate mdot_cooler = 40;

  parameter Integer N_seg = 10 "default number of discretized segments in one tube";
  
  // parameters based on input parameters.   
  parameter Modelica.SIunits.MassFlowRate mdot_bypass = mdot_main - mdot_comp "mdot in by pass path";
  parameter Real splitter_split_ratio                 = mdot_bypass/mdot_main "mass split ratio of splitter";


  // Thermal-hydraulic properties
  // default cp and rho for alloy X-750
  parameter Modelica.SIunits.Density rho_wall             = 8280 "density of wall, kg/m3";
  parameter Modelica.SIunits.SpecificHeatCapacity cp_wall = 431 "cp of wall, J/kg-K";
  // Thermal conductivity of LTR's metal wall 
  // material inconel_750
  //parameter Real table_k_LTR_wall[:, :]  = [149, 16.9; 316, 20.5; 538, 26.5; 649, 28.7; 760, 31.4; 871, 35.3];
  //parameter Real table_k_HTR_wall[:, :]  = [149, 16.9; 316, 20.5; 538, 26.5; 649, 28.7; 760, 31.4; 871, 35.3];
  // parameter Real table_k_metalwall[:,:]  = [293.15, 12.1; 373.15, 16.3; 773.15, 21.5];


  // **** START -  Definition of Key Point's state *****
  // names are mainly referred to heat exchangers such as LTR, HTR, Cooler and Heater  
  // except for st_recomp_in, st_recomp_out, which are refered to recompressor
  // points around cooler
  parameter Model.ThermoState st_cooler_hot_in(
    p    = p_comp_in,
    T    = T_LTR_hot_out,
    h    = medium_main.specificEnthalpy_pT(p = st_cooler_hot_in.p, T = st_cooler_hot_in.T),
    mdot = mdot_comp    
  ); 

  parameter Model.ThermoState st_cooler_hot_out(
    p    = p_comp_in,
    T    = T_cooler_hot_out,
    h    = medium_main.specificEnthalpy_pT(p = st_cooler_hot_out.p, T = st_cooler_hot_out.T),
    mdot = mdot_comp    
  );

  parameter Model.ThermoState st_cooler_cold_in(
    p    = p_amb,
    T    = T_amb,
    h    = medium_cooler.specificEnthalpy_pT(p = st_cooler_cold_in.p, T = st_cooler_cold_in.T),
    mdot = mdot_cooler    
  ); 
  
  parameter Model.ThermoState st_cooler_cold_out(
    p    = p_amb,
    T    = T_cooler_cold_out,
    h    = medium_cooler.specificEnthalpy_pT(p = st_cooler_cold_out.p, T = st_cooler_cold_out.T),
    mdot = mdot_cooler    
  );

  // points around LTR
  parameter Model.ThermoState st_LTR_hot_in(
    p    = p_comp_in,
    T    = T_HTR_hot_out,
    h    = medium_main.specificEnthalpy_pT(p = st_LTR_hot_in.p, T = st_LTR_hot_in.T),
    mdot = mdot_main  
  );

  parameter Model.ThermoState st_LTR_hot_out(
    p    = p_comp_in,
    T    = T_LTR_hot_out,
    h    = medium_main.specificEnthalpy_pT(p = st_LTR_hot_out.p, T = st_LTR_hot_out.T),
    mdot = mdot_main 
  );

  parameter Model.ThermoState st_LTR_cold_in(
    p    = p_comp_out,
    T    = T_LTR_cold_in,
    h    = medium_main.specificEnthalpy_pT(p = st_LTR_cold_in.p, T = st_LTR_cold_in.T),
    mdot = mdot_comp  
  );

  parameter Model.ThermoState st_LTR_cold_out(
    p    = p_comp_out,
    T    = T_HTR_cold_in,
    h    = medium_main.specificEnthalpy_pT(p = st_LTR_cold_out.p, T = st_LTR_cold_out.T),
    mdot = mdot_comp 
  );

  // points around HTR
  parameter Model.ThermoState st_HTR_hot_in(
    p    = p_comp_in,
    T    = T_HTR_hot_in,
    h    = medium_main.specificEnthalpy_pT(p = st_HTR_hot_in.p, T = st_HTR_hot_in.T),
    mdot = mdot_main  
  );

  parameter Model.ThermoState st_HTR_hot_out(
    p    = p_comp_in,
    T    = T_HTR_hot_out,
    h    = medium_main.specificEnthalpy_pT(p = st_HTR_hot_out.p, T = st_HTR_hot_out.T),
    mdot = mdot_main
  );

  parameter Model.ThermoState st_HTR_cold_in(
    p    = p_comp_out,
    T    = T_HTR_cold_in,
    h    = medium_main.specificEnthalpy_pT(p = st_HTR_cold_in.p, T = st_HTR_cold_in.T),
    mdot = mdot_main
  );

  parameter Model.ThermoState st_HTR_cold_out(
    p    = p_comp_out,
    T    = T_HTR_cold_out,
    h    = medium_main.specificEnthalpy_pT(p = st_HTR_cold_out.p, T = st_HTR_cold_out.T),
    mdot = mdot_main
  );

  parameter Model.ThermoState st_HTR_wall(
    p = p_comp_out,
    T = (T_HTR_cold_in + T_HTR_hot_in) / 2,
    h = medium_main.specificEnthalpy_pT(p = st_HTR_wall.p, T = st_HTR_wall.T),
    mdot = mdot_main
  );

  // points around heater
  parameter Model.ThermoState st_heater_hot_in(
    p    = p_heater,
    T    = T_heater_hot_in,
    h    = medium_heater.specificEnthalpy_pT(p = st_heater_hot_in.p, T = st_heater_hot_in.T),
    mdot = mdot_heater    
  ); 

  parameter Model.ThermoState st_heater_hot_out(
    p    = p_heater,
    T    = T_heater_hot_out,
    h    = medium_heater.specificEnthalpy_pT(p = st_heater_hot_out.p, T = st_heater_hot_out.T),
    mdot = mdot_heater    
  );

  parameter Model.ThermoState st_heater_cold_in(
    p    = p_comp_out,
    T    = T_HTR_cold_out,
    h    = medium_main.specificEnthalpy_pT(p = st_heater_cold_in.p, T = st_heater_cold_in.T),
    mdot = mdot_main    
  ); 
  
  parameter Model.ThermoState st_heater_cold_out(
    p    = p_comp_out,
    T    = T_heater_cold_out,
    h    = medium_main.specificEnthalpy_pT(p = p_comp_out, T = T_heater_cold_out),
    mdot = mdot_main    
  );

  parameter ThermoState st_recomp_in(
    p    = p_comp_in,
    T    = T_LTR_hot_out,
    h    = medium_main.specificEnthalpy_pT(p = p_comp_in, T = T_LTR_hot_out),
    mdot = mdot_bypass
    );  

  parameter ThermoState st_recomp_out(
    p    = p_comp_out,
    T    = T_HTR_cold_in,
    h    = medium_main.specificEnthalpy_pT(p = p_comp_out, T = T_recom_out),
    mdot = mdot_bypass
    );      
  // **** END   -  Definition of Key Point's state *****


  // **** START - Definittion of components' config *****
  // cooler
  // cooler's size of heat exchanger  
  // N_ch groups of fluid(cold, inner)-tube vs 1 gas(hot, outter) tube  

  // for ThermoPower.HeatExchanger
  parameter Modelica.SIunits.Radius r_i_cooler = 0.5e-3 "mm, tube's internal radius(single tube)";
  parameter Modelica.SIunits.Radius r_t_cooler = 0.7e-3 "tube's external radius(single tube)";
  parameter Modelica.SIunits.Radius r_o_cooler = 1e-3 "radius of external side of heat exchanger";
  parameter Modelica.SIunits.Length L_cooler   = 1 "m length of the tube/cooler";
  parameter Integer N_ch_cooler                = 50000;
  parameter Integer N_seg_cooler               = N_seg;

  parameter Real gamma_cold_cooler = (7728.47 + 14290.1) / 2; // calculuated with Test.Transient.TestDyn_cooler on 10MWe point
  parameter Real gamma_hot_cooler = (2046.73 + 3388) / 2; 
  
  //parameter Real UA_nom_cooler = 1.938761018e6;  // Arbitrary UA value, works for off-design RCBC sim

  // cfg for cooler's cold/fluid side
  parameter AreaGeometry ga_cooler_hot  = SetAreaGeometry_Circle(r = r_o_cooler);
  parameter AreaGeometry ga_cooler_cold = SetAreaGeometry_Circle(r = r_i_cooler);
  parameter AreaGeometry ga_cooler_wall = SetAreaGeometry_Tube(r_int = r_i_cooler, r_ext = r_t_cooler);
  parameter PathGeometry gp_cooler_hot  = SetPathGeometry(geo_area = ga_cooler_hot, L = L_cooler, N_seg = N_seg_cooler);
  parameter PathGeometry gp_cooler_cold = SetPathGeometry(geo_area = ga_cooler_cold, L = L_cooler, N_seg = N_seg_cooler);
  parameter PathGeometry gp_cooler_wall = SetPathGeometry(geo_area = ga_cooler_wall, L = L_cooler, N_seg = N_seg_cooler);
  
  // works for shell-tube heater, one fluid tube only
  parameter Model.HeatExchangerConfig cfg_cooler(
    cfg_hot(
      st_in    = st_cooler_hot_in,
      st_out   = st_cooler_hot_out,
      geo_area = ga_cooler_hot,
      geo_path = gp_cooler_hot,
      N_ch     = N_ch_cooler,
      u        = 0,
      UA_nom   = gamma_hot_cooler * gp_cooler_hot.A_surf,
      gamma_HE = gamma_hot_cooler 
    ),
    cfg_cold(
      st_in    = st_cooler_cold_in,
      st_out   = st_cooler_cold_out,
      geo_area = ga_cooler_cold,
      geo_path = gp_cooler_cold,
      N_ch     = N_ch_cooler,
      u        = 0,
      UA_nom   = gamma_cold_cooler * gp_cooler_cold.A_surf,
      gamma_HE = gamma_cold_cooler       
    ),
    cfg_fluid = cfg_cooler.cfg_cold,
    cfg_gas   = cfg_cooler.cfg_hot,
    cfg_wall(
      st_init  = st_cooler_hot_in,
      geo_area = ga_cooler_wall,
      geo_wall = gp_cooler_wall,
      // table_k  = table_k_LTR_wall,
      rho_mcm  = rho_wall * cp_wall,
      lambda   = 200
    )
  );
  /*
  // for SCOPE.PCHE
  parameter Modelica.SIunits.Radius r_i_cooler     = 1.5e-3 "mm tube's internal radius";
  parameter Modelica.SIunits.Radius r_t_cooler     = 1.5e-3 "tube's external radius";
  parameter Modelica.SIunits.Radius r_o_cooler     = 1.5e-3 "radius of external side of one group";
  parameter Modelica.SIunits.Radius w_sd_cooler    = 2.3e-3 "Width of the solid domain";
  parameter Modelica.SIunits.Radius h_sd_cooler    = 4.17e-3 "Height of the solid domain, containing one cold tube and one hot tube";
  parameter Modelica.SIunits.Length t_ch_cooler    = 0.51e-3 "thinckness between two channels";
  parameter Modelica.SIunits.Length l_pitch_cooler = 12.3e-3 "pitch length";
  parameter Modelica.SIunits.Length a_phi_cooler   = 35 "pitch angle";

  parameter Modelica.SIunits.Length L_cooler   = 2.5 "m";
  parameter Modelica.SIunits.Length l_wall_cooler  = 420e-3 "Length of wall, not necessarily equals to length of flow path";  
  parameter Integer N_seg_cooler               = N_seg;
  parameter Integer N_ch_cooler                = 30000;

  // cfg for cooler's cold/fluid side
  parameter AreaGeometry ga_cooler_hot  = SetAreaGeometry_SemiCircle(r = r_i_cooler);
  parameter AreaGeometry ga_cooler_cold = SetAreaGeometry_SemiCircle(r = r_o_cooler);
  parameter AreaGeometry ga_cooler_wall = SetAreaGeometry_Wall(r = r_i_cooler, w = w_sd_cooler, h = h_sd_cooler, p1 = t_ch_cooler);
  parameter PathGeometry gp_cooler_hot  = SetPathGeometry(geo_area = ga_cooler_hot, L = L_cooler, N_seg = N_seg_cooler);
  parameter PathGeometry gp_cooler_cold = SetPathGeometry(geo_area = ga_cooler_cold, L = L_cooler, N_seg = N_seg_cooler);
  parameter PathGeometry gp_cooler_wall = SetPathGeometry(geo_area = ga_cooler_wall, L = L_cooler, N_seg = N_seg_cooler, p1 = l_wall_cooler);
  
  parameter Model.HeatExchangerConfig cfg_cooler(
    cfg_hot(
      st_in    = st_cooler_hot_in,
      st_out   = st_cooler_hot_out,
      geo_area = ga_cooler_hot,
      geo_path = gp_cooler_hot,
      N_ch     = N_ch_cooler,
      u        = 0,
      l_pitch  = l_pitch_cooler,
      a_phi    = a_phi_cooler
    ),
    cfg_cold(
      st_in    = st_cooler_cold_in,
      st_out   = st_cooler_cold_out,
      geo_area = ga_cooler_cold,
      geo_path = gp_cooler_cold,
      N_ch     = N_ch_cooler,
      u        = 0,
      l_pitch  = l_pitch_cooler,
      a_phi    = a_phi_cooler
    ),
    cfg_wall(
      st_init  = st_cooler_hot_in,
      geo_area = ga_cooler_wall,
      geo_wall = gp_cooler_wall,
      // table_k  = table_k_cooler_wall,
      rho_mcm  = rho_wall * cp_wall,
      lambda   = 200
    )
  );
    */
  // LTR
  // LTR's's size of heat exchanger gas - gas
  // N_ch_LTR groups of fluid(cold, inner)-tube-gas(hot, outter) tubes 
  parameter Modelica.SIunits.Radius r_i_LTR     = 1.5e-3 "mm tube's internal radius";
  parameter Modelica.SIunits.Radius r_t_LTR     = 1.5e-3 "tube's external radius";
  parameter Modelica.SIunits.Radius r_o_LTR     = 1.5e-3 "radius of external side of one group";
  parameter Modelica.SIunits.Radius w_sd_LTR    = 2.3e-3 "Width of the solid domain";
  parameter Modelica.SIunits.Radius h_sd_LTR    = 4.17e-3 "Height of the solid domain, containing one cold tube and one hot tube";
  parameter Modelica.SIunits.Length t_ch_LTR    = 0.51e-3 "thinckness between two channels";
  parameter Modelica.SIunits.Length l_pitch_LTR = 12.3e-3 "pitch length";
  parameter Modelica.SIunits.Length a_phi_LTR   = 35 "pitch angle";

  parameter Modelica.SIunits.Length L_LTR   = 2.5 "m";
  parameter Modelica.SIunits.Length l_wall_LTR  = 420e-3 "Length of wall, not necessarily equals to length of flow path";  
  parameter Integer N_seg_LTR               = N_seg;
  parameter Integer N_ch_LTR                = 30000;

  // cfg for LTR's cold/fluid side
  parameter AreaGeometry ga_LTR_hot  = SetAreaGeometry_SemiCircle(r = r_i_LTR);
  parameter AreaGeometry ga_LTR_cold = SetAreaGeometry_SemiCircle(r = r_o_LTR);
  parameter AreaGeometry ga_LTR_wall = SetAreaGeometry_Wall(r = r_i_LTR, w = w_sd_LTR, h = h_sd_LTR, p1 = t_ch_LTR);
  parameter PathGeometry gp_LTR_hot  = SetPathGeometry(geo_area = ga_LTR_hot, L = L_LTR, N_seg = N_seg_LTR);
  parameter PathGeometry gp_LTR_cold = SetPathGeometry(geo_area = ga_LTR_cold, L = L_LTR, N_seg = N_seg_LTR);
  parameter PathGeometry gp_LTR_wall = SetPathGeometry(geo_area = ga_LTR_wall, L = L_LTR, N_seg = N_seg_LTR, p1 = l_wall_LTR);
  
  parameter Model.HeatExchangerConfig cfg_LTR(
    cfg_hot(
      st_in    = st_LTR_hot_in,
      st_out   = st_LTR_hot_out,
      geo_area = ga_LTR_hot,
      geo_path = gp_LTR_hot,
      N_ch     = N_ch_LTR,
      u        = 0,
      l_pitch  = l_pitch_LTR,
      a_phi    = a_phi_LTR
    ),
    cfg_cold(
      st_in    = st_LTR_cold_in,
      st_out   = st_LTR_cold_out,
      geo_area = ga_LTR_cold,
      geo_path = gp_LTR_cold,
      N_ch     = N_ch_LTR,
      u        = 0,
      l_pitch  = l_pitch_LTR,
      a_phi    = a_phi_LTR
    ),
    cfg_fluid = cfg_LTR.cfg_cold,
    cfg_gas   = cfg_LTR.cfg_hot,
    cfg_wall(
      st_init  = st_LTR_hot_in,
      geo_area = ga_LTR_wall,
      geo_wall = gp_LTR_wall,
      // table_k  = table_k_LTR_wall,
      rho_mcm  = rho_wall * cp_wall,
      lambda   = 200
    )
  );

  // HTR
  // HTR's's size of heat exchanger gas - gas
  // N_ch_HTR groups of fluid(cold, inner)-tube-gas(hot, outter) tubes 
  parameter Modelica.SIunits.Radius r_i_HTR     = 1.5e-3 "mm tube's internal radius";
  parameter Modelica.SIunits.Radius r_t_HTR     = 1.5e-3 "tube's external radius";
  parameter Modelica.SIunits.Radius r_o_HTR     = 1.5e-3 "radius of external side of one group";
  parameter Modelica.SIunits.Radius w_sd_HTR    = 2.3e-3 "Width of the solid domain";
  parameter Modelica.SIunits.Radius h_sd_HTR    = 4.17e-3 "Height of the solid domain, containing one cold tube and one hot tube";
  parameter Modelica.SIunits.Length t_ch_HTR    = 0.51e-3 "thinckness between two channels";
  parameter Modelica.SIunits.Length l_pitch_HTR = 12.3e-3 "pitch length";
  parameter Modelica.SIunits.Length a_phi_HTR   = 35 "pitch angle";

  parameter Modelica.SIunits.Length L_HTR      = 2.5 "m";
  parameter Modelica.SIunits.Length l_wall_HTR = 420e-3 "Length of wall, not necessarily equals to length of flow path";
  parameter Integer N_seg_HTR                  = N_seg;
  parameter Integer N_ch_HTR                   = 30000;

  // cfg for HTR's cold/fluid side
  parameter AreaGeometry ga_HTR_hot  = SetAreaGeometry_SemiCircle(r = r_i_HTR);
  parameter AreaGeometry ga_HTR_cold = SetAreaGeometry_SemiCircle(r = r_o_HTR);
  parameter AreaGeometry ga_HTR_wall = SetAreaGeometry_Wall(r = r_i_HTR, w = w_sd_HTR, h = h_sd_HTR, p1 = t_ch_HTR);
  parameter PathGeometry gp_HTR_hot  = SetPathGeometry(geo_area = ga_HTR_hot, L = L_HTR, N_seg = N_seg_HTR);
  parameter PathGeometry gp_HTR_cold = SetPathGeometry(geo_area = ga_HTR_cold, L = L_HTR, N_seg = N_seg_HTR);
  parameter PathGeometry gp_HTR_wall = SetPathGeometry(geo_area = ga_HTR_wall, L = L_HTR, N_seg = N_seg_HTR, p1 = l_wall_HTR);
  
  parameter Model.HeatExchangerConfig cfg_HTR(
    cfg_hot(
      st_in    = st_HTR_hot_in,
      st_out   = st_HTR_hot_out,
      geo_area = ga_HTR_hot,
      geo_path = gp_HTR_hot,
      N_ch     = N_ch_HTR,
      u        = 0,
      l_pitch  = l_pitch_HTR,
      a_phi    = a_phi_HTR
    ),
    cfg_cold(
      st_in    = st_HTR_cold_in,
      st_out   = st_HTR_cold_out,
      geo_area = ga_HTR_cold,
      geo_path = gp_HTR_cold,
      N_ch     = N_ch_HTR,
      u        = 0,
      l_pitch  = l_pitch_HTR,
      a_phi    = a_phi_HTR
    ),
    cfg_fluid = cfg_HTR.cfg_cold,
    cfg_gas   = cfg_HTR.cfg_hot,
    cfg_wall(
      st_init  = st_HTR_wall,
      geo_area = ga_HTR_wall,
      geo_wall = gp_HTR_wall,
      // table_k  = table_k_HTR_wall,
      rho_mcm  = rho_wall * cp_wall,
      lambda   = 200
    )
  );

  // heater
  // heater's's size of heat exchanger gas - gas
  // N_ch_heater groups of fluid(cold, inner)-tube-gas(hot, outter) tubes 
  parameter Modelica.SIunits.Radius r_i_heater = 1e-3 "mm tube's internal radius";
  parameter Modelica.SIunits.Radius r_t_heater = 2e-3 "tube's external radius";
  parameter Modelica.SIunits.Radius r_o_heater = 3e-3 "radius of external side of one group";
  parameter Modelica.SIunits.Length L_heater   = 1 "m";
  parameter Integer N_seg_heater               = N_seg;
  parameter Integer N_ch_heater                = 10000;

  // cfg for heater's cold/fluid side  
  parameter AreaGeometry ga_heater_hot  = SetAreaGeometry_Circle(r = r_i_heater);
  parameter AreaGeometry ga_heater_cold = SetAreaGeometry_Circle(r = r_o_heater);
  parameter AreaGeometry ga_heater_wall = SetAreaGeometry_Tube(r_int = r_i_heater, r_ext = r_o_heater);
  parameter PathGeometry gp_heater_hot  = SetPathGeometry(geo_area = ga_heater_hot, L = L_heater, N_seg = N_seg_heater);
  parameter PathGeometry gp_heater_cold = SetPathGeometry(geo_area = ga_heater_cold, L = L_heater, N_seg = N_seg_heater);
  parameter PathGeometry gp_heater_wall = SetPathGeometry(geo_area = ga_heater_wall, L = L_heater, N_seg = N_seg_heater);  
  
  parameter Real gamma_cold_heater = (20407.3 + 20514.8) / 2; // calculuated with Test.Transient.TestDyn_heater on 10MWe point
  parameter Real gamma_hot_heater = (1633.5 + 1650.54) / 2;  

  parameter Model.HeatExchangerConfig cfg_heater(
    cfg_hot(
      st_in    = st_heater_hot_in,
      st_out   = st_heater_hot_out,
      geo_area = ga_heater_hot,
      geo_path = gp_heater_hot,
      N_ch     = N_ch_heater,
      u        = 0,
      UA_nom   = gamma_hot_heater * gp_heater_hot.A_surf,
      gamma_HE = gamma_hot_heater
    ),
    cfg_cold(
      st_in    = st_heater_cold_in,
      st_out   = st_heater_cold_out,
      geo_area = ga_heater_cold,
      geo_path = gp_heater_cold,
      N_ch     = N_ch_heater,
      u        = 0,
      UA_nom   = gamma_cold_heater * gp_heater_cold.A_surf,
      gamma_HE = gamma_cold_heater
    ),
    cfg_gas = cfg_heater.cfg_cold,
    cfg_fluid = cfg_heater.cfg_hot,
    cfg_wall(
      st_init  = st_heater_hot_in,
      geo_area = ga_heater_wall,
      geo_wall = gp_heater_wall,
      // table_k  = table_k_LTR_wall,
      rho_mcm  = rho_wall * cp_wall,
      lambda   = 200
    )
  );

  parameter TurbomachineryConfig cfg_comp(
    st_in  = st_cooler_hot_out,
    st_out = st_LTR_cold_in,
    N      = Ns_comp,
    T_nom  = 0,
    eta    = eta_comp
  );

  parameter TurbomachineryConfig cfg_recomp(
    st_in  = st_recomp_in,
    st_out = st_recomp_out,
    N      = Ns_recomp,
    T_nom  = 0,
    eta    = eta_recomp
  );

  parameter Model.TurbomachineryConfig cfg_turb(
    st_in  = st_heater_cold_out,
    st_out = st_HTR_hot_in,
    N      = Ns_turb,
    T_nom  = 0,
    eta    = eta_turb
  );  

  parameter Model.SplitterConfig cfg_splitter(
    st_in    = st_LTR_hot_out,
    st_out   = st_cooler_hot_in,
    st_split = st_recomp_in,
    gamma    = splitter_split_ratio
  );

  parameter Model.SplitterConfig cfg_merger(
    st_in    = st_LTR_cold_out,
    st_out   = st_HTR_cold_in,
    st_split = st_recomp_out,
    gamma    = splitter_split_ratio
  );

  // **** END - Definittion of Components' Config *****  

  
end RCBCycleConfig;
