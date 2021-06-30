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
      thermo(rho_mcm = rho_mcm, lambda = 12) // refer http://www.matweb.com/search/datasheet_print.aspx?matguid=9612aa3272134531b8b33eb80e61a1af&n=1 for INCONEL® Alloy X-750
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
      thermo(rho_mcm = rho_mcm, lambda = 12) // refer http://www.matweb.com/search/datasheet_print.aspx?matguid=9612aa3272134531b8b33eb80e61a1af&n=1 for INCONEL® Alloy X-750
    ),
    cfg_HTR_hot(
      geo(
        V = V_HTR, // r_o_HTR^2 * pi * L_HTR * N_ch_HTR - cfg_HTR_tube.geo.V , 
        A_ex = A_HTR, // * N_ch_HTR, 
        L = L_HTR, 
        d = r_HTR * 2),
      thermo(gamma_he = 1.938761018e6/A_HTR "200")
    ),
    N_ch_HTR = 422585,
    cfg_heater_cold(
      geo(
        V = V_heater, // * N_ch_heater, 
        A_ex = A_heater, // * N_ch_heater, exchange surface between fluid-tube
        L = L_h, 
        d = r_heater * 2), // calculate by steps the python code
      thermo(gamma_he = 1.938761018e6/A_heater "200")
    ),
    cfg_heater_tube(
      geo(
        V = V_heater * 2, //r_t_heater^2 * pi * L_h * N_ch_heater - cfg_heater_cold.geo.V,
        A_ex = A_heater, //  * N_ch_heater, // assume thickness of tube approximately 0
        L = L_h, 
        d = r_heater * 2),        
      thermo(rho_mcm = rho_mcm, lambda = 12) // refer http://www.matweb.com/search/datasheet_print.aspx?matguid=9612aa3272134531b8b33eb80e61a1af&n=1 for INCONEL® Alloy X-750
    ),
    cfg_heater_hot(
      geo(
        V = V_heater, // r_o_heater^2 * pi * L_h * N_ch_heater - cfg_heater_tube.geo.V , 
        A_ex = A_heater, // * N_ch_heater, 
        L = L_h, 
        d = r_heater * 2),
      thermo(gamma_he = 1.938761018e6/A_heater "200")
    ),
    N_ch_h = 422585           
  );

  import Modelica.SIunits.Conversions.{from_degC, from_deg};
  import Modelica.SIunits.{Temperature, Pressure, SpecificEnthalpy};
  import Util = Utilities.Util;
  import Steps.Utilities.CoolProp.PropsSI; 
  import Steps.Components.PCHEGeoParam;  
  
 
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
  
  parameter Real r_heater     = 1e-3;  
  parameter Real peri_heater  = (pi + 2) * r_heater;
  parameter Real V_heater     = r_heater^2 * pi * L_h / 2;
  parameter Real A_heater     = peri_heater * L_h;  
  
  parameter Modelica.SIunits.Length pitch = 12.3e-3 "pitch length";
  parameter Real phi = 36 "pitch angle, degree";  
  
  // default cp and rho for alloy X-750
  parameter Modelica.SIunits.Density rho_wall = 8280 "density of wall, kg/m3";
  parameter Modelica.SIunits.SpecificHeatCapacity cp_wall = 431 "cp of wall, J/kg-K";
  
  parameter Real rho_mcm = rho_wall * cp_wall "Metal heat capacity per unit volume [J/m^3.K]";

end PBConfig_PCHE;
