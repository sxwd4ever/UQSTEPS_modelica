within Steps.Media;
package CO2 "Ideal gas \"CO2\" from NASA Glenn coefficients"
  /*extends Modelica.Media.IdealGases.Common.SingleGasNasa(
    mediumName=data.name,
    data=Modelica.Media.IdealGases.Common.SingleGasesData.CO2,
    fluidConstants={Modelica.Media.IdealGases.Common.FluidData.CO2});*/
  
  extends Modelica.Media.Interfaces.PartialPureSubstance(
     ThermoStates=Modelica.Media.Interfaces.Choices.IndependentVariables.pT,
     redeclare final record FluidConstants =
        Modelica.Media.Interfaces.Types.IdealGas.FluidConstants,
     mediumName=data.name,
     substanceNames={data.name},
     singleState=false,
     Temperature(min=200, max=6000, start=500, nominal=500),
     SpecificEnthalpy(start=if Functions.referenceChoice==ReferenceEnthalpy.ZeroAt0K then data.H0 else
        if Functions.referenceChoice==ReferenceEnthalpy.UserDefined then Functions.h_offset else 0, nominal=1.0e5),
     Density(start=10, nominal=10),
     AbsolutePressure(start=10e5, nominal=10e5));
  /*
  extends Modelica.Media.Interfaces.PartialPureSubstance(
    ThermoStates=Modelica.Media.Interfaces.Choices.IndependentVariables.pT,
    FluidConstants = Modelica.Media.Interfaces.Types.IdealGas.FluidConstants,
    mediumName="CO2",
    substanceNames={"CO2"},
    singleState=false,
    //final reducedX = true, 
    final fixedX = true, 
    Temperature(min=200, max=6000, start=500, nominal=500),
    SpecificEnthalpy(start=0, nominal=1.0e5),
    Density(start=10, nominal=10),
    AbsolutePressure(start=10e5, nominal=10e5));  
  */
  
  import CP = Steps.Utilities.CoolProp;  
  import SI = Modelica.SIunits;
  import Modelica.Math;
  import Modelica.Media.Interfaces.Choices.ReferenceEnthalpy;  
  import Modelica.Media.IdealGases.Common;
  import Modelica.Media.IdealGases.Common.Functions;
  import Modelica.Media.Air.ReferenceAir.Air_Utilities;
  
  constant Modelica.Media.IdealGases.Common.DataRecord data(
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
    R=188.9244822140674) "Data record of ideal gas substance";
  
  redeclare record extends ThermodynamicState "A selection of variables that uniquely defines the thermodynamic state"
      AbsolutePressure p "Absolute pressure of medium";
      Temperature T "Temperature of medium";
  end ThermodynamicState;

  redeclare replaceable model extends BaseProperties(
   T(stateSelect=if preferredMediumStates then StateSelect.prefer else StateSelect.default),
   p(stateSelect=if preferredMediumStates then StateSelect.prefer else StateSelect.default)) "Base properties of medium"
    Modelica.SIunits.SpecificEntropy s;
  equation
      //state.p = p;
      // state.T = T; 
    MM = data.MM;
    R = data.R;
    /*
    h = Modelica.Media.IdealGases.Common.Functions.h_T(
            data, T,
            Modelica.Media.IdealGases.Common.Functions.excludeEnthalpyOfFormation,
            Modelica.Media.IdealGases.Common.Functions.referenceChoice,
            Modelica.Media.IdealGases.Common.Functions.h_offset);
    u = h - R*T;
    */
    h = specificEnthalpy(state);
    u = h - p / d;
    // Has to be written in the form d=f(p,T) in order that static
    // state selection for p and T is possible
    d = density(state);
    s = specificEntropy(state);
    // connect state with BaseProperties
    state.T = T;
    state.p = p;
    /*      
      state = newState_pT(p, T);
      d = CP.PropsSI("D", "P", p, "T", T, mediumName);
      h = CP.PropsSI("H", "P", p, "T", T, mediumName);
      s = CP.PropsSI("S", "P", p, "T", T, mediumName);
      u = h - p / d;
      MM = 0.044;
      R = 8.3144 / MM; */
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
    annotation(Inline=true,smoothOrder=2);
  end setState_pTX;

  function setState_ph
    "Return thermodynamic state from p and h"
    extends Modelica.Icons.Function;
    input AbsolutePressure p "Pressure";
    input SpecificEnthalpy h "Specific enthalpy";
    input FixedPhase phase=0
      "2 for two-phase, 1 for one-phase, 0 if not known";
    output ThermodynamicState state "Thermodynamic state record";
  algorithm
    state := setState_phX(p, h, fill(0, 0), phase);
  end setState_ph;  
   
  redeclare function setState_phX
  "Return thermodynamic state as function of p, h and composition X"
    extends Modelica.Icons.Function;
    input AbsolutePressure p "Pressure";
    input SpecificEnthalpy h "Specific enthalpy";
    input MassFraction X[:]=reference_X "Mass fractions";
    output ThermodynamicState state;
  algorithm
    state := ThermodynamicState(p=p,T=CP.PropsSI("T", "P", p, "H", h, mediumName));
    annotation(Inline=true,smoothOrder=2);
  end setState_phX;  
  
  redeclare function setState_psX
  "Return thermodynamic state as function of p, s and composition X"
    extends Modelica.Icons.Function;
    input AbsolutePressure p "Pressure";
    input SpecificEntropy s "Specific entropy";
    input MassFraction X[:]=reference_X "Mass fractions";
    output ThermodynamicState state;
  algorithm
    state := ThermodynamicState(p=p,T=CP.PropsSI("T", "P", p, "S", s, mediumName));
    annotation(Inline=true,smoothOrder=2);
  end setState_psX;

  redeclare function setState_dTX
  "Return thermodynamic state as function of d, T and composition X"
    extends Modelica.Icons.Function;
    input Density d "Density";
    input Temperature T "Temperature";
    input MassFraction X[:]=reference_X "Mass fractions";
    output ThermodynamicState state;
  algorithm
    state := ThermodynamicState(p=d*data.R*T,T=T);
    annotation(Inline=true,smoothOrder=2);
  end setState_dTX;

  redeclare function extends setSmoothState
  "Return thermodynamic state so that it smoothly approximates: if x > 0 then state_a else state_b"
  algorithm
    state := ThermodynamicState(p=Media.Common.smoothStep(x, state_a.p, state_b.p, x_small),
                                T=Media.Common.smoothStep(x, state_a.T, state_b.T, x_small));
    annotation(Inline=true,smoothOrder=2);
  end setSmoothState;

  redeclare function extends pressure "Return pressure of ideal gas"
  algorithm
    p := state.p;
    annotation(Inline=true,smoothOrder=2);
  end pressure;

  redeclare function extends temperature "Return temperature of ideal gas"
  algorithm
    T := state.T;
    annotation(Inline=true,smoothOrder=2);
  end temperature;

  redeclare function extends density "Return density of ideal gas"
  algorithm
    d := CP.PropsSI("D", "P", state.p, "T", state.T, mediumName);
    //d := state.p/(data.R*state.T);
    annotation(Inline=true,smoothOrder=2);
  end density;

  redeclare function extends specificEnthalpy "Return specific enthalpy"
    extends Modelica.Icons.Function;
  algorithm
    // h := Modelica.Media.IdealGases.Common.Functions.h_T(data,state.T);    
    
    h := CP.PropsSI("H", "P", state.p, "T", state.T, mediumName);
    annotation(Inline=true,smoothOrder=2);
  end specificEnthalpy;
  
  function specificEnthalpy_pT 
    input AbsolutePressure p "Pressure";
    input Temperature T "Temperature";
    output SpecificEnthalpy h "Specific enthalpy";
	algorithm
    h := specificEnthalpy(setState_pT(p, T));
	end specificEnthalpy_pT;

  redeclare function extends specificInternalEnergy
    "Return specific internal energy"
    extends Modelica.Icons.Function;
  algorithm
    //u := CP.PropsSI("U", "P", state.p, "T", state.T, mediumName);
    u := state.h - state.p / state.d;
    annotation(Inline=true,smoothOrder=2);
  end specificInternalEnergy;

  redeclare function extends specificEntropy "Return specific entropy"
    extends Modelica.Icons.Function;
  algorithm
    s := CP.PropsSI("V", "P", state.p, "T", state.T, mediumName);
    //s := Modelica.Media.IdealGases.Common.Functions.s0_T(
    //         data, state.T) - data.R*Modelica.Math.log(state.p/reference_p);
    annotation(Inline=true,smoothOrder=2);    

  end specificEntropy;

  redeclare function extends specificGibbsEnergy "Return specific Gibbs energy"
    extends Modelica.Icons.Function;
  algorithm
    g := CP.PropsSI("GMASS", "P", state.p, "T", state.T, mediumName);
    annotation(Inline=true,smoothOrder=2);
  end specificGibbsEnergy;

  redeclare function extends specificHelmholtzEnergy
    "Return specific Helmholtz energy"
    extends Modelica.Icons.Function;
  algorithm
    f := CP.PropsSI("HELMHOLTZMASS", "P", state.p, "T", state.T, mediumName);
    annotation(Inline=true,smoothOrder=2);
  end specificHelmholtzEnergy;

  redeclare function extends specificHeatCapacityCp
    "Return specific heat capacity at constant pressure"
  algorithm
    cp := CP.PropsSI("CP0MASS", "P", state.p, "T", state.T, mediumName);
    annotation(Inline=true,smoothOrder=2);
  end specificHeatCapacityCp;

  redeclare function extends specificHeatCapacityCv
    "Compute specific heat capacity at constant volume from temperature and gas data"
  algorithm
    cv := CP.PropsSI("CVMASS", "P", state.p, "T", state.T, mediumName);
    annotation(Inline=true,smoothOrder=2);
  end specificHeatCapacityCv;

  redeclare function extends isentropicExponent "Return isentropic exponent"
  algorithm
    gamma := specificHeatCapacityCp(state)/specificHeatCapacityCv(state);
    annotation(Inline=true,smoothOrder=2);
  end isentropicExponent;

  redeclare function extends velocityOfSound "Return velocity of sound"
    extends Modelica.Icons.Function;
  algorithm
    a := sqrt(max(0,data.R*state.T*Modelica.Media.IdealGases.Common.Functions.cp_T(
                                        data, state.T)/specificHeatCapacityCv(state)));
    annotation(Inline=true,smoothOrder=2);
  end velocityOfSound;
/*
  function isentropicEnthalpyApproximation
    "Approximate method of calculating h_is from upstream properties and downstream pressure"
    extends Modelica.Icons.Function;
    input SI.Pressure p2 "Downstream pressure";
    input ThermodynamicState state "Properties at upstream location";
    input Boolean exclEnthForm=Functions.excludeEnthalpyOfFormation
      "If true, enthalpy of formation Hf is not included in specific enthalpy h";
    input ReferenceEnthalpy refChoice=Functions.referenceChoice
      "Choice of reference enthalpy";
    input SpecificEnthalpy h_off=Functions.h_offset
      "User defined offset for reference enthalpy, if referenceChoice = UserDefined";
    output SI.SpecificEnthalpy h_is "Isentropic enthalpy";
  protected
    IsentropicExponent gamma =  isentropicExponent(state) "Isentropic exponent";
  algorithm
    h_is := Modelica.Media.IdealGases.Common.Functions.h_T(
                data,state.T,exclEnthForm,refChoice,h_off) +
      gamma/(gamma - 1.0)*state.p/density(state)*((p2/state.p)^((gamma - 1)/gamma) - 1.0);
    annotation(Inline=true,smoothOrder=2);
  end isentropicEnthalpyApproximation;
*/
  redeclare function extends isentropicEnthalpy "Return isentropic enthalpy"
  input Boolean exclEnthForm=Functions.excludeEnthalpyOfFormation
      "If true, enthalpy of formation Hf is not included in specific enthalpy h";
  input ReferenceEnthalpy refChoice=Functions.referenceChoice
      "Choice of reference enthalpy";
  input SpecificEnthalpy h_off=Functions.h_offset
      "User defined offset for reference enthalpy, if referenceChoice = UserDefined";
  protected
    Modelica.SIunits.SpecificEntropy s;
  algorithm
    s := specificEntropy(refState);
    h_is := specificEntropy(setState_psX(p=p_downstream, s=s)); 
    annotation(Inline=true,smoothOrder=2);
  end isentropicEnthalpy;

  redeclare function extends isobaricExpansionCoefficient
    "Returns overall the isobaric expansion coefficient beta"
  algorithm
    beta := CP.PropsSI("ISOBARIC_EXPANSION_COEFFICIENT", "P", state.p, "T", state.T, mediumName); // 1/state.T;
    annotation(Inline=true,smoothOrder=2);
  end isobaricExpansionCoefficient;

  redeclare function extends isothermalCompressibility
    "Returns overall the isothermal compressibility factor"
  algorithm
    kappa := CP.PropsSI("ISOTHERMAL_COMPRESSIBILITY", "P", state.p, "T", state.T, mediumName); // 1.0/state.p;
    annotation(Inline=true,smoothOrder=2);
  end isothermalCompressibility;

  redeclare function extends density_derp_T
    "Returns the partial derivative of density with respect to pressure at constant temperature"
  algorithm
    ddpT := CP.PropsSI("d(DMASS)/d(P)|T", "P", state.p, "T", state.T, mediumName); //1/(state.T*data.R);
    annotation(Inline=true,smoothOrder=2);
  end density_derp_T;

  redeclare function extends density_derT_p
    "Returns the partial derivative of density with respect to temperature at constant pressure"
  algorithm
    ddTp := CP.PropsSI("d(DMASS)/d(T)|P", "P", state.p, "T", state.T, mediumName); //-state.p/(state.T*state.T*data.R);
    annotation(Inline=true,smoothOrder=2);
  end density_derT_p;

  redeclare function extends density_derX
    "Returns the partial derivative of density with respect to mass fractions at constant pressure and temperature"
  algorithm
    dddX := fill(0,nX);
    annotation(Inline=true,smoothOrder=2);
  end density_derX;

  redeclare replaceable function extends dynamicViscosity "Dynamic viscosity"
  algorithm
    eta := CP.PropsSI("V", "P", state.p, "T", state.T, mediumName);
    annotation (smoothOrder=2);
  end dynamicViscosity;

  redeclare replaceable function extends thermalConductivity
    "Thermal conductivity of gas"
  //  input IdealGases.Common.DataRecord data "Ideal gas data";
  //  input Integer method=Functions.methodForThermalConductivity "1: Eucken Method, 2: Modified Eucken Method";
  algorithm
    lambda := CP.PropsSI("CONDUCTIVITY", "P", state.p, "T", state.T, mediumName);
    annotation (smoothOrder=2);
  end thermalConductivity;

  redeclare function extends molarMass "Return the molar mass of the medium"
  algorithm
    MM := data.MM;
    annotation(Inline=true,smoothOrder=2);
  end molarMass;
  
  function extends density_derh_p
    "Density derivative by specific enthalpy"
  algorithm    
    ddhp := CP.PropsSI("d(DMASS)/d(Hmass)|P", "P", state.p, "Hmass", specificEnthalpy_pT(state.p, state.T), mediumName);
    annotation (Inline=true);
  end density_derh_p;

  function extends density_derp_h "Density derivative by pressure"
  algorithm
    ddph := CP.PropsSI("d(DMASS)/d(P)|Hmass", "P", state.p, "Hmass", specificEnthalpy_pT(state.p, state.T), mediumName);
    annotation (Inline=true);
  end density_derp_h;  
/*
  function T_h "Compute temperature from specific enthalpy"
    extends Modelica.Icons.Function;
    input SpecificEnthalpy h "Specific enthalpy";
    output Temperature T "Temperature";

  protected
  package Internal
      "Solve h(data,T) for T with given h (use only indirectly via temperature_phX)"
    extends Modelica.Media.Common.OneNonLinearEquation;
    redeclare record extends f_nonlinear_Data
        "Data to be passed to non-linear function"
      extends Modelica.Media.IdealGases.Common.DataRecord;
    end f_nonlinear_Data;

    redeclare function extends f_nonlinear
    algorithm
        y := Modelica.Media.IdealGases.Common.Functions.h_T(
                 f_nonlinear_data,x);
    end f_nonlinear;

    // Dummy definition has to be added for current Dymola
    redeclare function extends solve
    end solve;
  end Internal;

  algorithm
    T := Internal.solve(h, 200, 6000, 1.0e5, {1}, data);
  end T_h;

  function T_ps "Compute temperature from pressure and specific entropy"
    extends Modelica.Icons.Function;
    input AbsolutePressure p "Pressure";
    input SpecificEntropy s "Specific entropy";
    output Temperature T "Temperature";

  protected
  package Internal
      "Solve h(data,T) for T with given h (use only indirectly via temperature_phX)"
    extends Modelica.Media.Common.OneNonLinearEquation;
    redeclare record extends f_nonlinear_Data
        "Data to be passed to non-linear function"
      extends Modelica.Media.IdealGases.Common.DataRecord;
    end f_nonlinear_Data;

    redeclare function extends f_nonlinear
    algorithm
        y := Modelica.Media.IdealGases.Common.Functions.s0_T(
                  f_nonlinear_data,x)- data.R*Modelica.Math.log(p/reference_p);
    end f_nonlinear;

    // Dummy definition has to be added for current Dymola
    redeclare function extends solve
    end solve;
  end Internal;

  algorithm
    T := Internal.solve(s, 200, 6000, p, {1}, data);
  end T_ps;
*/
/*
// the functions below are not strictly necessary, there are just here for compatibility reasons

  function dynamicViscosityLowPressure
    "Dynamic viscosity of low pressure gases"
    extends Modelica.Icons.Function;
    input SI.Temp_K T "Gas temperature";
    input SI.Temp_K Tc "Critical temperature of gas";
    input SI.MolarMass M "Molar mass of gas";
    input SI.MolarVolume Vc "Critical molar volume of gas";
    input Real w "Acentric factor of gas";
    input DipoleMoment mu "Dipole moment of gas molecule";
    input Real k =  0.0 "Special correction for highly polar substances";
    output SI.DynamicViscosity eta "Dynamic viscosity of gas";
  protected
    parameter Real Const1_SI=40.785*10^(-9.5)
      "Constant in formula for eta converted to SI units";
    parameter Real Const2_SI=131.3/1000.0
      "Constant in formula for mur converted to SI units";
    Real mur=Const2_SI*mu/sqrt(Vc*Tc)
      "Dimensionless dipole moment of gas molecule";
    Real Fc=1 - 0.2756*w + 0.059035*mur^4 + k
      "Factor to account for molecular shape and polarities of gas";
    Real Tstar "Dimensionless temperature defined by equation below";
    Real Ov "Viscosity collision integral for the gas";

  algorithm
    eta := Functions.dynamicViscosityLowPressure(T,Tc,M,Vc,w,mu,k);
    annotation (smoothOrder=2,
                Documentation(info="<html>
<p>
The used formula are based on the method of Chung et al (1984, 1988) referred to in ref [1] chapter 9.
The formula 9-4.10 is the one being used. The Formula is given in non-SI units, the following conversion constants were used to
transform the formula to SI units:
</p>

<ul>
<li> <b>Const1_SI:</b> The factor 10^(-9.5) =10^(-2.5)*1e-7 where the
     factor 10^(-2.5) originates from the conversion of g/mol->kg/mol + cm^3/mol->m^3/mol
      and the factor 1e-7 is due to conversion from microPoise->Pa.s.</li>
<li>  <b>Const2_SI:</b> The factor 1/3.335641e-27 = 1e-3/3.335641e-30
      where the factor 3.335641e-30 comes from debye->C.m and
      1e-3 is due to conversion from cm^3/mol->m^3/mol</li>
</ul>

<h4>References:</h4>
<p>
[1] Bruce E. Poling, John E. Prausnitz, John P. O'Connell, \"The Properties of Gases and Liquids\" 5th Ed. Mc Graw Hill.
</p>

<h4>Author</h4>
<p>T. Skoglund, Lund, Sweden, 2004-08-31</p>

</html>"));
  end dynamicViscosityLowPressure;

  function thermalConductivityEstimate
    "Thermal conductivity of polyatomic gases(Eucken and Modified Eucken correlation)"
    extends Modelica.Icons.Function;
    input Interfaces.PartialMedium.SpecificHeatCapacity Cp
      "Constant pressure heat capacity";
    input Interfaces.PartialMedium.DynamicViscosity eta "Dynamic viscosity";
    input Integer method(min=1,max=2)=1
      "1: Eucken Method, 2: Modified Eucken Method";
    input IdealGases.Common.DataRecord data "Ideal gas data";
    output Interfaces.PartialMedium.ThermalConductivity lambda
      "Thermal conductivity [W/(m.k)]";
  algorithm
    lambda := Functions.thermalConductivityEstimate(Cp,eta,method,data);
    annotation (smoothOrder=2);
  end thermalConductivityEstimate;

  annotation (*/
  
end CO2;
