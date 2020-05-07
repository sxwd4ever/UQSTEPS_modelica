within Steps.Media;
package SCO2 "supercritical CO2"

  extends Modelica.Media.Interfaces.PartialMedium(ThermoStates = Modelica.Media.Interfaces.Choices.IndependentVariables.ph, final mediumName = "CO2", final substanceNames = {"CO2"}, final singleState = false, final reducedX = true, final fixedX = true, Temperature(min = 20, max = 873.15, start = 500));
  import CP = Steps.Utilities.CoolProp;  
  import SI = Modelica.SIunits;
  
  type Temperature = Real (
      final quantity="ThermodynamicTemperature",
      final unit="K",
      min=273.15,
      max=1200,
      start=300,
      nominal=300,
      displayUnit="degC")
      "Temperature of the medium; add of min and max constraints"
      annotation(absoluteValue=true);
  
  redeclare record extends ThermodynamicState "A selection of variables that uniquely defines the thermodynamic state"
      AbsolutePressure p "Absolute pressure of medium";
      Temperature T "Temperature of SCO2";
  end ThermodynamicState;
  
  model CO2_pT "Base properties of medium"
    extends BaseProperties(final standardOrderComponents = true);
  equation
      state.p = p;
      state.T = T;     
      d = CP.PropsSI("D", "P", state.p, "T", state.T, mediumName);
      h = CP.PropsSI("H", "P", state.p, "T", state.T, mediumName);
      u = h - p / d;
      MM = 0.044;
      R = 8.3144 / MM; 
  end CO2_pT;
    
  redeclare function setState_pTX "Return thermodynamic state from p, T, and X or Xi"
    extends Modelica.Icons.Function;
    input AbsolutePressure p "Pressure";
    input Temperature T "Temperature";
    input MassFraction X[:] = reference_X "Mass fractions";
    output ThermodynamicState state "Thermodynamic state record";
  algorithm
    state := ThermodynamicState(p = p, T = T);
  end setState_pTX;

  redeclare function setState_phX "Return thermodynamic state from p, h, and X or Xi"
    extends Modelica.Icons.Function;
    input AbsolutePressure p "Pressure";
    input SpecificEnthalpy h "Specific enthalpy";
    input MassFraction X[:] = reference_X "Mass fractions";
    output ThermodynamicState state "Thermodynamic state record";
  algorithm
    state := ThermodynamicState(p = p, T = CP.PropsSI("T", "P", p, "H", h, mediumName));
  end setState_phX;

  redeclare function setState_psX "Return thermodynamic state from p, s, and X or Xi"
    extends Modelica.Icons.Function;
    input AbsolutePressure p "Pressure";
    input SpecificEntropy s "Specific entropy";
    input MassFraction X[:] = reference_X "Mass fractions";
    output ThermodynamicState state "Thermodynamic state record";
  algorithm
    state := ThermodynamicState(p = p, T = CP.PropsSI("T", "P", p, "S", s, mediumName));
  end setState_psX;

  redeclare function extends pressure "Return pressure"
  algorithm
      p := state.p;
  end pressure;

  redeclare function extends temperature "Return temperature"
  algorithm
      T := state.T;
  end temperature;

  redeclare function extends specificEnthalpy "Return specific enthalpy"
  algorithm
      h := CP.PropsSI("H", "P", state.p, "T", state.T, mediumName);
  end specificEnthalpy;

  redeclare function extends density "Return density"
  protected
      outer Modelica.Blocks.Types.ExternalCombiTable2D tableIDd_p_h;
  algorithm
      d := CP.PropsSI("D", "P", state.p, "T", state.T, mediumName);
  end density;

  redeclare function extends specificInternalEnergy "Return specific internal energy"
  protected
      outer Modelica.Blocks.Types.ExternalCombiTable2D tableIDu_p_h;
  algorithm
      u := CP.PropsSI("U", "P", state.p, "T", state.T, mediumName);
  end specificInternalEnergy;

  redeclare function extends dynamicViscosity "Return dynamic viscosity"
  algorithm
      eta := CP.PropsSI("V", "P", state.p, "T", state.T, mediumName);
  end dynamicViscosity;

  redeclare function extends thermalConductivity "Return thermal conductivity"
  algorithm
      lambda := CP.PropsSI("CONDUCTIVITY", "P", state.p, "T", state.T, mediumName);
  end thermalConductivity;


  redeclare function extends specificEntropy "Return specific entropy"
  algorithm
      s := CP.PropsSI("S", "P", state.p, "T", state.T, mediumName);
  end specificEntropy;

  redeclare function extends specificHelmholtzEnergy "Return specific Helmholtz energy"
  algorithm
      f := CP.PropsSI("HELMHOLTZMASS", "P", state.p, "T", state.T, mediumName);
  end specificHelmholtzEnergy;

  redeclare function extends specificHeatCapacityCp "Return specific heat capacity at constant pressure"
  algorithm
      cp := CP.PropsSI("CP0MASS", "P", state.p, "T", state.T, mediumName);
  end specificHeatCapacityCp;

  redeclare function extends specificHeatCapacityCv "Return specific heat capacity at constant volume"
  algorithm
      cv := CP.PropsSI("CVMASS", "P", state.p, "T", state.T, mediumName);
  end specificHeatCapacityCv;  
end SCO2;

