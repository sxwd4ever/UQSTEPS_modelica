within Steps;

package Model "Data model definition for consistent parameter configuration"
  extends Modelica.Icons.Package;
  import ThermoPower.Choices.Flow1D.FFtypes;    

  // ****** START - record named EntityXXXX will be descrepted *****
  record EntityConfig "Record for the ThermoState, align with the struct definition in C"
    parameter ThermoState st_start "thermo state for start point";
    parameter ThermoState st_nom "nominal state values";
    parameter EntityGeoParam geo "Geometry values";
    parameter EntityThermoParam thermo "ThermoDynamic parameters";
  end EntityConfig;

  record EntityThermoParam "Record of entity's thermodynamic parameters, such as conductance, pressure drop, etc."
    parameter Real mdot "mass flow rate";
    parameter ThermoPower.Choices.Flow1D.FFtypes FFtype = ThermoPower.Choices.Flow1D.FFtypes.NoFriction "Friction Factor Type";
    parameter Real dp "pressure drop";
    parameter Real Kf "Nominal hydraulic resistance coefficient";
    parameter Real rho "density";
    parameter Real C_f "Fanning friction factor";
    parameter Real K_fc "Friction factor correction coefficient";
    parameter Real K_fl "Linear friction factor";
    parameter Real rho_mcm "heat capacity per unit volume";
    parameter Real lambda "Thermal conductivity of the entity (density by specific heat capacity)";
    parameter Real gamma_he = 200 "Constant heat transfer coefficient";
    parameter Real UAnom = 1.0 "Universal/Total heat transfer coefficient";
  end EntityThermoParam;

  record EntityGeoParam "Record Entity's Geo parameter"
    parameter Real L "length";
    parameter Real V "volume";
    parameter Real H "Elevation of outlet over inlet";
    parameter Real A_ex "Area of Heat Exchange surface";
    parameter Real peri_ex "Perimeter of heat transfer surface (single tube)";
    parameter Real peri_ex_hyd "Wet perimeter (single tube)";
    parameter Real d "diameter for a cylinder shaped entity - pipe, tube, etc.";
    parameter Real d_hyd "Hydraulic Diameter (single tube)";
    parameter Integer N_ch = 1 "(optional) number of channels";
    parameter Integer N_seg = 10 "number of segments";
    parameter Real e_rel "Relative roughness (ratio roughness/diameter)";
  end EntityGeoParam;

  // ****** END - record named EntityXXXX will be descrepted *****

  record ThermoState "Record for the ThermoState, align with the struct definition in C"
    Real T "Temperature, K";
    Real p "Pressure, pa";
    Modelica.SIunits.SpecificEnthalpy h "specific enthalpy, J/kg";
    Real mdot "mass flow rate kg/s";
    // optional id, to specify which point this state is assigned to
    // because of the mapping constraintns of Modelica and external functions
    // integer id is used instead of string(char *) name
    // more detailed naming specificaiton should be provided
    // hot_in=0, cold_in=1, hot_out=2, cold_out=3
    Integer id = -1;
  end ThermoState;

  record AreaGeometry
    Modelica.SIunits.Area A "Area of cross section";
    Modelica.SIunits.Length peri "perimeter";
    Modelica.SIunits.Length peri_wet "wet perimeter";
    Modelica.SIunits.Length d "Diameter, valid for circle, semi-circle, otherwise 0";
    Modelica.SIunits.Length w "Width, valid for square rectangle or triangle, otherwise 0";
    Modelica.SIunits.Length h "Height, valid for square rectangle or triangle, otherwise 0";
    Modelica.SIunits.Length d_hyd "Hydraulic Diameter";
    Real p1 "customerized parameters, app determined";
    Real p2 "customerized parameters, app determined";    
  end AreaGeometry;

  record PathGeometry
    Modelica.SIunits.Length l "path length";
    Modelica.SIunits.Volume V "volume of the path";
    Modelica.SIunits.Area A_cr "cross-sectional area";
    Modelica.SIunits.Area A_surf "surface area";
    Real e_rel "(Relative) Roughness";
    Integer N_seg "Number of segments, valid for a discritized path";
    Real p1 "customerized parameters, app determined";
    Real p2 "customerized parameters, app determined";      
  end PathGeometry;

  record FlowConfig
    ThermoState st_in "inlet state";
    ThermoState st_out "outlet state";
    AreaGeometry geo_area "Cross sectional area geometry, for area of single path";
    PathGeometry geo_path "Path geometry, for single path";
    Integer N_ch "number of channels";
    Modelica.SIunits.Velocity u "flow velocity";
    Real l_pitch "pitch length";
    Real a_phi "pitch angle";
  end FlowConfig;

  record WallConfig
    ThermoState st_init "global averaged config for wall initialization";
    AreaGeometry geo_area "Wall cross section geometry";
    PathGeometry geo_wall "Wall geometry along the length";
    // material inconel_750
    Real table_k[:, :] = [149, 16.9; 316, 20.5; 538, 26.5; 649, 28.7; 760, 31.4; 871, 35.3] "Thermal conductivity of metal wall";
    Real rho_mcm = 7900 * 578.05;
    Real lambda;
  end WallConfig;

  record TurbomachineryConfig
    ThermoState st_in "inlet state";
    ThermoState st_out "outlet state";
    Real N "rotational speed";
    Real T_nom "nominal torque";
    Real eta "efficiency";
  end TurbomachineryConfig;

  record HeatExchangerConfig
    FlowConfig cfg_hot "hot side config";
    FlowConfig cfg_cold "cold side config";
    WallConfig cfg_wall "wall configuration";
  end HeatExchangerConfig;

  record SplitterConfig
    ThermoState st_in "inlet config";
    ThermoState st_out "outlet config";
    ThermoState st_split "split inlet/outlet config";
    Real gamma "split ratio";
  end SplitterConfig;

  partial function SetAreaGeometry
    input Real r;
    input Real L;
    output AreaGeometry geo;
    constant Real pi = Modelica.Constants.pi;
  end SetAreaGeometry;

  function SetAreaGeometry_Circle
    extends SetAreaGeometry;
  protected 
    Real peri_in;    
  algorithm
    peri_in := (pi + 2) * r;
    geo := AreaGeometry(
      d        = r * 2,
      peri     = peri_in,
      peri_wet = peri_in,
      A        = peri_in * L,
      w        = 0,
      h        = 0,
      d_hyd    = r * 2,
      p1       = 0,
      p2       = 0
    );
  end SetAreaGeometry_Circle;

  function SetAreaGeometry_Wall
    "Specific function to set the wall geometry of a PCHE"
    input Real w;
    input Real h;
    input Real p1 = 0;
    input Real p2 = 0;    
    extends SetAreaGeometry;
  protected 
    Real peri_in;    
  algorithm
    peri_in := (pi + 2) * r;
    geo := AreaGeometry(
      d        = r * 2,
      peri     = peri_in,
      peri_wet = peri_in,
      A        = peri_in * L,
      w        = w,
      h        = h,
      d_hyd    = r * 2,
      p1       = p1,
      p2       = p2
    );
  end SetAreaGeometry_Wall;

  function SetPathGeometry
    input Real r;
    input Real L;
    input Integer N_seg;
    input Real e_rel = 0;
    input Real p1    = 0;
    input Real p2    = 0;
    output PathGeometry geo;
    constant Real pi = Modelica.Constants.pi;
  protected
    Real peri_in;
  algorithm
    peri_in := (pi + 2) * r;
    geo := PathGeometry(
      l      = L,
      V      = r ^ 2 * pi * L / 2,
      A_cr   = r ^ 2 * pi,
      A_surf = r * 2 * pi * L,
      e_rel  = e_rel,
      N_seg  = N_seg,
      p1     = p1,
      p2     = p2
    );
  end SetPathGeometry;  

end Model;
