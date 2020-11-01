within Steps.Media;
package SCO2 "supercritical CO2"
  /*
  extends Modelica.Media.IdealGases.Common.SingleGasNasa(
  mediumName="Carbon Dioxide",
  data=Modelica.Media.IdealGases.Common.SingleGasesData.CO2,
  fluidConstants={Modelica.Media.IdealGases.Common.FluidData.CO2});
  */
  
  extends Modelica.Media.Interfaces.PartialPureSubstance(
    ThermoStates=Modelica.Media.Interfaces.Choices.IndependentVariables.pT,
    FluidConstants = Modelica.Media.Interfaces.Types.IdealGas.FluidConstants,
    mediumName="CO2",
    substanceNames={"CO2"},
    singleState=false,
    final reducedX = true, 
    final fixedX = true, 
    Temperature(min=200, max=6000, start=500, nominal=500),
    SpecificEnthalpy(start=0, nominal=1.0e5),
    Density(start=10, nominal=10),
    AbsolutePressure(start=10e5, nominal=10e5));  
  
  constant Modelica.Media.IdealGases.Common.DataRecord CO2(
    name="CO2",
    MM=0.0440095,
    Hf=-8941478.544405185,
    H0=212805.6215135368,
    Tlimit=1000,
    alow={49436.5054,-626.411601,5.30172524,0.002503813816,-2.127308728e-007,-7.68998878e-010,
        2.849677801e-013},
    blow={-45281.9846,-7.04827944},
    ahigh={117696.2419,-1788.791477,8.29152319,-9.22315678e-005,4.86367688e-009,
        -1.891053312e-012,6.330036589999999e-016},
    bhigh={-39083.5059,-26.52669281},
    R=188.9244822140674);
  
  
  
  import CP = Steps.Utilities.CoolProp;  
  import SI = Modelica.SIunits;
  
  /*
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
  */
  

  redeclare record extends ThermodynamicState "A selection of variables that uniquely defines the thermodynamic state"
      AbsolutePressure p "Absolute pressure of medium";
      Temperature T "Temperature of medium";
  end ThermodynamicState;

  
  /*
  record FullThermodynamicState
    "Full state info for a medium, for data transfer purpose"
    
    Modelica.SIunits.Temperature T;
    Modelica.SIunits.AbsolutePressure p;
    Modelica.SIunits.SpecificEnthalpy h;    
    Modelica.SIunits.SpecificEntropy s;
    Modelica.SIunits.Density d;
    
  end FullThermodynamicState;
  */
  
  redeclare replaceable model extends BaseProperties(final standardOrderComponents = true) "Base properties of medium"
    Modelica.SIunits.SpecificEntropy s;
  equation
      //state.p = p;
      // state.T = T;     
      state = newState_pT(p, T);
      d = CP.PropsSI("D", "P", p, "T", T, mediumName);
      h = CP.PropsSI("H", "P", p, "T", T, mediumName);
      s = CP.PropsSI("S", "P", p, "T", T, mediumName);
      u = h - p / d;
      MM = 0.044;
      R = 8.3144 / MM; 
  end BaseProperties;
    
  function getFullState
    "get the full state of this medium"
    input BaseProperties prop;
    output FullThermodynamicState state;
    algorithm
    
    state.T := prop.T;
    state.p := prop.p;
    state.h := prop.h;
    state.d := prop.d;
    state.s := prop.s;
    
  end getFullState;
  
  function newState_pT
    "factory method using to create new state, for debug purpose"
    
    input Real p;
    input Real T;
    
    output ThermodynamicState state;
   
  algorithm
    state.p := p;
    state.T := T; 

    print("p=" + String(p) + " T=" + String(T) +"\n");
  end newState_pT;
    
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
  
  function setState_ph
    "Return thermodynamic state from p and h"
    extends Modelica.Icons.Function;
    input AbsolutePressure p "Pressure";
    input SpecificEnthalpy h "Specific enthalpy";
    input FixedPhase phase=0
      "2 for two-phase, 1 for one-phase, 0 if not known";
    output ThermodynamicState state "Thermodynamic state record";
  algorithm
    state := setState_phX(
            p,
            h,
            fill(0, 0),
            phase);
  end setState_ph;

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
  
  redeclare function extends density_derp_T
    "Returns the partial derivative of density with respect to pressure at constant temperature"
  algorithm
    ddpT := 1/(state.T*CO2.R);
    annotation(Inline=true,smoothOrder=2);
  end density_derp_T;

  redeclare function extends density_derT_p
    "Returns the partial derivative of density with respect to temperature at constant pressure"
  algorithm
    ddTp := -state.p/(state.T*state.T*CO2.R);
    annotation(Inline=true,smoothOrder=2);
  end density_derT_p;  
  
  redeclare function extends density_derX
    "Returns the partial derivative of density with respect to mass fractions at constant pressure and temperature"
  algorithm
    dddX := fill(0,nX);
    annotation(Inline=true,smoothOrder=2);
  end density_derX;  
end SCO2;
