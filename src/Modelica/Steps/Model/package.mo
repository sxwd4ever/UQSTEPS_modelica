within Steps;

package Model "Data model definition for consistent parameter configuration"
  extends Modelica.Icons.Package;
  import ThermoPower.Choices.Flow1D.FFtypes;

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

  record ThermoState "Record for the ThermoState, align with the struct definition in C"
    Real T "Temperature, K";
    Real p "Pressure, pa";
    Real h "specific enthalpy, J/kg";
    Real mdot "mass flow rate kg/s";
    // optional id, to specify which point this state is assigned to
    // because of the mapping constraintns of Modelica and external functions
    // integer id is used instead of string(char *) name
    // more detailed naming specificaiton should be provided
    // hot_in=0, cold_in=1, hot_out=2, cold_out=3
    Integer id = -1;
  end ThermoState;










end Model;
