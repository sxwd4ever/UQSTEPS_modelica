within Steps.Model;

model PBConfig_PCHE
  "Preset param set for Off-design power block test with PCHE as recuperator"
  extends PBConfiguration(
    cfg_LTR_cold(
      geo(
        V = V_LTR, // * N_ch_LTR, 
        A_ex = A_LTR, // * N_ch_LTR, exchange surface between fluid-tube
        L = L_LTR, 
        d = r_LTR * 2
        ), // calculate by steps the python code
      thermo(gamma_he = 1.938761018e6/A_LTR "200")
    ),
    cfg_LTR_tube(
      geo(
        V = V_LTR * 2, //r_t_LTR^2 * pi * L_LTR * N_ch_LTR - cfg_LTR_cold.geo.V,
        A_ex = A_LTR, //  * N_ch_LTR, // assume thickness of tube approximately 0
        L = L_LTR, 
        d = r_LTR * 2),
      thermo(rho_mcm = 8.28e6 * 0.431, lambda = 12) // refer http://www.matweb.com/search/datasheet_print.aspx?matguid=9612aa3272134531b8b33eb80e61a1af&n=1 for INCONEL® Alloy X-750
    ),
    cfg_LTR_hot(
      geo(
        V = V_LTR, // r_o_LTR^2 * pi * L_LTR * N_ch_LTR - cfg_LTR_tube.geo.V , 
        A_ex = A_LTR, // * N_ch_LTR, 
        L = L_LTR, 
        d = r_LTR * 2),
      thermo(gamma_he = 1.938761018e6/A_LTR "200")
    ),
    N_ch_LTR = 332449, 
    cfg_HTR_cold(
      geo(
        V = V_HTR, // * N_ch_HTR, 
        A_ex = A_HTR, // * N_ch_HTR, exchange surface between fluid-tube
        L = L_HTR, 
        d = r_HTR * 2), // calculate by steps the python code
      thermo(gamma_he = 1.938761018e6/A_HTR "200")
    ),
    cfg_HTR_tube(
      geo(
        V = V_HTR * 2, //r_t_HTR^2 * pi * L_HTR * N_ch_HTR - cfg_HTR_cold.geo.V,
        A_ex = A_HTR, //  * N_ch_HTR, // assume thickness of tube approximately 0
        L = L_HTR, 
        d = r_HTR * 2),        
      thermo(rho_mcm = 8.28e6 * 0.431, lambda = 12) // refer http://www.matweb.com/search/datasheet_print.aspx?matguid=9612aa3272134531b8b33eb80e61a1af&n=1 for INCONEL® Alloy X-750
    ),
    cfg_HTR_hot(
      geo(
        V = V_HTR, // r_o_HTR^2 * pi * L_HTR * N_ch_HTR - cfg_HTR_tube.geo.V , 
        A_ex = A_HTR, // * N_ch_HTR, 
        L = L_HTR, 
        d = r_HTR * 2),
      thermo(gamma_he = 1.938761018e6/A_HTR "200")
    ),
    N_ch_HTR = 422585      
  );

  import Modelica.SIunits.Conversions.{from_degC, from_deg};
  import Modelica.SIunits.{Temperature, Pressure, SpecificEnthalpy};
  import Util = Utilities.Util;
  import Steps.Utilities.CoolProp.PropsSI; 
  import Steps.Components.PCHEGeoParam;  

  parameter PCHEGeoParam geo_HTR(
    // pitch length, m
    pitch = 12e-3,
    // pitch angle
    phi = from_deg((180 - 108) /2),
    // length of pche, m
    L = 1000e-3, // 2860e-3,
    // Diameter of semi_circular, m
    d = 2e-3,
    // number of channels
    N_ch = integer(94e3),
    // number of segments
    N_seg = 50);
    
   parameter PCHEGeoParam geo_LTR(
    // pitch length, m
    pitch = 12e-3,
    // pitch angle
    phi = from_deg((180 - 108) /2),
    // length of pche, m
    L = 3270e-3,
    // Diameter of semi_circular, m
    d = 2e-3,
    // number of channels
    N_ch = integer(125e3),
    // number of segments
    N_seg = 50);
 
  // test configuration for hot/cold/tube side of heatexchanger
  
  // In following calculation, V, A_ex are account for single tube/channel, not for total
  // check the Calculation in ThemoPower.PowerPlants.HRSG.Components.HE to understand the meaning of 
  // exsurface_G/F and extSurfaceTub, *Vol     
  parameter Real r_HTR = 1e-3;  
  parameter Real p_HTR = (pi + 2) * r_HTR;
  parameter Real V_HTR = r_HTR^2 * pi * L_HTR / 2;
  parameter Real A_HTR = p_HTR * L_HTR;

  parameter Real r_LTR = 1e-3;
  parameter Real p_LTR = (pi + 2) * r_LTR;
  parameter Real V_LTR = r_LTR^2 * pi * L_LTR / 2;
  parameter Real A_LTR = p_LTR * L_LTR;

end PBConfig_PCHE;
