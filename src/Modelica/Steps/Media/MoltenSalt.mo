within Steps.Media;

package MoltenSalt
  package MoltenSalt_pT "Molten Salt (60% NaNO3, 40% KNO3 by weight), explicit in p and T"
    /* For a new medium, make a copy of this package and remove
                                                                        the "partial" keyword from the package definition above.
                                                                        The statement below extends from PartialMedium and sets some
                                                                        package constants. Provide values for these constants
                                                                        that are appropriate for your medium model. Note that other
                                                                        constants (such as nX, nXi) are automatically defined by
                                                                        definitions given in the base class Interfaces.PartialMedium"
                                                                    */
    /* IMPORTANT extends from PartialPureSubstance instead of PartialMedium - Xin Sui 20210208 */
    extends Modelica.Media.Interfaces.PartialPureSubstance(
    ThermoStates = Modelica.Media.Interfaces.Choices.IndependentVariables.pT, 
    final mediumName = "MoltenSalt", 
    final substanceNames = {"NaNO3", "KNO3"}, 
    final singleState = false, final reducedX = true, 
    final fixedX = true, 
    Temperature(min = 573.15, max = 873.15, start = 800));
    
    import Steps.Media.MoltenSalt.MoltenSalt_utilities.*;
    // Provide medium constants here
    //constant SpecificHeatCapacity cp_const=123456
    // "Constant specific heat capacity at constant pressure";
    /* 
    The vector substanceNames is mandatory, as the number of
    substances is determined based on its size. Here we assume
    a single-component medium.
    singleState is true if u and d do not depend on pressure, but only
    on a thermal variable (temperature or enthalpy). Otherwise, set it
    to false.
    For a single-substance medium, just set reducedX and fixedX to true, and there's
    no need to bother about medium compositions at all. Otherwise, set
    final reducedX = true if the medium model has nS-1 independent mass
    fraction, or reducedX = false if the medium model has nS independent
    mass fractions (nS = number of substances).
    If a mixture has a fixed composition set fixedX=true, otherwise false.
    The modifiers for reducedX and fixedX should normally be final
    since the other equations are based on these values.
  
    It is also possible to redeclare the min, max, and start attributes of
    Medium types, defined in the base class Interfaces.PartialMedium
    (the example of Temperature is shown here). Min and max attributes
    should be set in accordance to the limits of validity of the medium
    model, while the start attribute should be a reasonable default value
    for the initialization of nonlinear solver iterations 
    */
    /* 
    Provide an implementation of model BaseProperties,
    that is defined in PartialMedium. Select two independent
    variables from p, T, d, u, h. The other independent
    variables are the mass fractions "Xi", if there is more
    than one substance. Provide 3 equations to obtain the remaining
    variables as functions of the independent variables.
    It is also necessary to provide two additional equations to set
    the gas constant R and the molar mass MM of the medium.
    Finally, the thermodynamic state vector, defined in the base class
    Interfaces.PartialMedium.BaseProperties, should be set, according to
    its definition (see ThermodynamicState below).
    The computation of vector X[nX] from Xi[nXi] is already included in
    the base class Interfaces.PartialMedium.BaseProperties, so it should not
    be repeated here.
    The code fragment below is for a single-substance medium with
    p,T as independent variables.
    */
    /*
        redeclare record extends ThermodynamicState
            "A selection of variables that uniquely defines the thermodynamic state"
            AbsolutePressure p "Absolute pressure of medium";
            SpecificEnthalpy h "Specific enthalpy";
            annotation (Documentation(info="<html>
      
                </html>"));
        end ThermodynamicState;
        */

    redeclare replaceable record extends ThermodynamicState "A selection of variables that uniquely defines the thermodynamic state"
        AbsolutePressure p "Absolute pressure of medium";
        Temperature T "Temperature of medium";
      annotation(
        Documentation(info = "<html>
                
                </html>"));
    end ThermodynamicState;

    redeclare model extends BaseProperties(final standardOrderComponents = true) "Base properties of medium"
      equation
        d = rho_T(T);
        h = h_T(state.T);
        u = h - p / d;
        MM = 0.091438;
        R = 8.3144 / MM;
        state.p = p;
        T = state.T;
    end BaseProperties;

    /* Provide implementations of the following optional properties.
    If not available, delete the corresponding function.
    The record "ThermodynamicState" contains the input arguments
    of all the function and is defined together with the used
    type definitions in PartialMedium. The record most often contains two of the
    variables "p, T, d, h" (e.g., medium.T)
    */

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
      state := ThermodynamicState(p = p, T = T_h(h));
    end setState_phX;
    
    /*
    replaceable function setState_ph "Return thermodynamic state from p and h"
      extends Modelica.Icons.Function;
      input AbsolutePressure p "Pressure";
      input SpecificEnthalpy h "Specific enthalpy";
      output ThermodynamicState state "Thermodynamic state record";
    algorithm
      state := ThermodynamicState(p = p, T = T_h(h));
    end setState_ph;
    */
    
    redeclare function setState_ph
      "Return thermodynamic state from p and h"
      extends Modelica.Icons.Function;
      input AbsolutePressure p "Pressure";
      input SpecificEnthalpy h "Specific enthalpy";
      input FixedPhase phase=0
        "2 for two-phase, 1 for one-phase, 0 if not known";
      output ThermodynamicState state "Thermodynamic state record";
      //output Modelica.Media.Interfaces.PartialMedium.ThermodynamicState state "Thermodynamic state record";
    algorithm
      state.p := p;
      state.T := T_h(h);
      //state.d := density(p, T=T_h(h)); //setState_phX(p, h, fill(0, 0));
      //state.phase := 1;
    end setState_ph;

    redeclare function setState_psX "Return thermodynamic state from p, s, and X or Xi"
      extends Modelica.Icons.Function;
      input AbsolutePressure p "Pressure";
      input SpecificEntropy s "Specific entropy";
      input MassFraction X[:] = reference_X "Mass fractions";
      output ThermodynamicState state "Thermodynamic state record";
    algorithm
      state := ThermodynamicState(p = p, T = T_s(s));
    end setState_psX;

    redeclare function setState_dTX "Return thermodynamic state from d, T, and X or Xi"
      extends Modelica.Icons.Function;
      input Density d "Pressure";
      input Temperature T "Specific entropy";
      input MassFraction X[:] = reference_X "Mass fractions";
      output ThermodynamicState state "Thermodynamic state record";
    algorithm
      state := ThermodynamicState(p = p_rho(d), T = T);
    end setState_dTX;

    redeclare function extends pressure "Return pressure"
      algorithm
        p := state.p;
      annotation(
        Inline = true);
    end pressure;

    redeclare function extends temperature "Return temperature"
      algorithm
        T := state.T;
      annotation(
        Inline = true);
    end temperature;

    redeclare function extends specificEnthalpy "Return specific enthalpy"
      algorithm
        h := h_T(state.T);
      annotation(
        Inline = true);
    end specificEnthalpy;

    redeclare function extends density "Return density"
      algorithm
        d := rho_T(state.T);
      annotation(
        Inline = true);
    end density;

    redeclare function extends specificInternalEnergy "Return specific internal energy"
      algorithm
        u := state.h - state.p / rho_T(state.T);
      annotation(
        Inline = true);
    end specificInternalEnergy;

    redeclare function extends dynamicViscosity "Return dynamic viscosity"
      algorithm
        eta := eta_T(state.T);
      annotation(
        Documentation(info = "<html>
  
                </html>"));
    end dynamicViscosity;

    redeclare function extends thermalConductivity "Return thermal conductivity"
      algorithm
        lambda := lamda_T(state.T);
      annotation(
        Documentation(info = "<html>
  
                </html>"));
    end thermalConductivity;

    redeclare function extends specificEntropy "Return specific entropy"
      algorithm
        s := s_T(state.T);
      annotation(
        Documentation(info = "<html>
  
                </html>"));
    end specificEntropy;

    redeclare function extends specificGibbsEnergy "Return specific Gibbs energy"
      algorithm
        g := gibbs_T(state.T);
      annotation(
        Documentation(info = "<html>
  
            </html>"));
    end specificGibbsEnergy;

    redeclare function extends specificHelmholtzEnergy "Return specific Helmholtz energy"
      algorithm
        f := helmholtz_pT(T = state.T, p = state.p);
      annotation(
        Documentation(info = "<html>
  
            </html>"));
    end specificHelmholtzEnergy;

    redeclare function extends specificHeatCapacityCp "Return specific heat capacity at constant pressure"
      algorithm
        cp := cp_T(state.T);
      annotation(
        Documentation(info = "<html>
  
                </html>"));
    end specificHeatCapacityCp;

    redeclare function extends specificHeatCapacityCv "Return specific heat capacity at constant volume"
      algorithm
        cv := 0;
      annotation(
        Documentation(info = "<html>
  
                </html>"));
    end specificHeatCapacityCv;

    redeclare function extends isentropicExponent "Return isentropic exponent"
        extends Modelica.Icons.Function;

      algorithm
        gamma := 1;
      annotation(
        Documentation(info = "<html>
  
                </html>"));
    end isentropicExponent;

    redeclare function extends isentropicEnthalpy "Return isentropic enthalpy"
      algorithm
        h_is := 0;
// To be completed
      annotation(
        Documentation(info = "<html>
  
                </html>"));
    end isentropicEnthalpy;

    redeclare function extends velocityOfSound "Return velocity of sound"
        extends Modelica.Icons.Function;

      algorithm
        a := 0;
      annotation(
        Documentation(info = "<html>
  
                </html>"));
    end velocityOfSound;

    redeclare function extends isobaricExpansionCoefficient "Return overall the isobaric expansion coefficient beta"
      algorithm
        beta := beta_T(state.T);
      annotation(
        Documentation(info = "<html>
  
                    </html>"));
    end isobaricExpansionCoefficient;

    redeclare function extends isothermalCompressibility "Return overall the isothermal compressibility factor"
      algorithm
        kappa := kappa_T(state.T);
      annotation(
        Documentation(info = "<html>
  
                </html>"));
    end isothermalCompressibility;

    function enthalpyOfVaporization "Return vaporization enthalpy of condensing fluid"
      extends Modelica.Icons.Function;
      input ThermodynamicState state "Thermodynamic state record";
      output SpecificEnthalpy r0 "Vaporization enthalpy";
    algorithm
      r0 := h_fg_T(state.T);
      annotation(
        Documentation(info = "<html>
  
                </html>"));
    end enthalpyOfVaporization;

    redeclare function extends density_derp_h "Density derivative by pressure"
    algorithm
      ddph := 0 "Incompressible"; 
      // CP.PropsSI("d(DMASS)/d(P)|Hmass", "P", state.p, "Hmass", specificEnthalpy_pT(state.p, state.T), mediumName);
      annotation (Inline=true);
    end density_derp_h; 

    redeclare function extends density_derh_p
      "Return density derivative w.r.t. specific enthalpy at constant pressure"
      algorithm
        ddhp:= drho_dT_T(state.T) * dT_dh_h(h_T(state.T));
    end density_derh_p;

    redeclare function extends density_derT_p "Return density derivative w.r.t. temperature at constant pressure"
      algorithm
        ddTp := drho_dT_T(state.T);
      annotation(
        Documentation(info = "<html>
  
                </html>"));
    end density_derT_p;
    annotation(
      Documentation(info = "<html>
    <p><span style=\"font-family: Arial,sans-serif;\">Calculation of fluid properties for the mixture common molten salt (60&percnt; NaNO<sub>3</sub> and 40&percnt; KNO<sub>3</sub>) in the fluid region of 573.15 to 873.15 Kelvin.The use of this molten salt is usally used in solar thermal system due to its high thermal capacity and the high range of temperature. </span></p>
    <p><span style=\"font-family: Arial,sans-serif;\">This package of thermodynamic properties is explicit for pressure and specific enthalpy, however it is based in functions where the dependency is exclusively in the temperature. </span></p>
    <p><b><span style=\"font-family: Arial,sans-serif; color: #008000;\">Restriction</span></b></p>
    <p><span style=\"font-family: Arial,sans-serif;\">The functions provided by this package shall be used inside of the restricted limits according to the referenced literature. </span></p>
    <ul>
    <li><b><span style=\"font-family: Arial,sans-serif;\">573.15 Kelvin &le; T &le; 873.15 Kelvin </span></b></li>
    <li><b><span style=\"font-family: Arial,sans-serif;\">explicit for pressure and specific enthalpy </span></b></li>
    </ul>
    <p><b><span style=\"font-family: Arial,sans-serif;\">References</span></b> </p>
    <p style=\"margin-left: 30px;\">Zavoico, A. B. (2001). <i>Solar Power Tower - Design Basis Document</i>. <i>Technical Report SAND2001-2100</i>. Alburquerque, USA. Retrieved from http://prod.sandia.gov/techlib/access-control.cgi/2001/012100.pdf</p>
    <p style=\"margin-left: 30px;\">Ferri, R., Cammi, A., &AMP; Mazzei, D. (2008). Molten salt mixture properties in RELAP5 code for thermodynamic solar applications. <i>International Journal of Thermal Sciences</i>, <i>47</i>(12), 1676&ndash;1687. http://doi.org/10.1016/j.ijthermalsci.2008.01.007</p>
    </html>"));
  end MoltenSalt_pT;

  package MoltenSalt_utilities
    function beta_T "Isobaric expansion coefficient of molten salt"
      extends Modelica.Icons.Function;
      input Modelica.SIunits.Temperature T "Temperature";
      output Modelica.Media.Interfaces.Types.IsobaricExpansionCoefficient beta_isb "Isobaric expansion coefficient";
    algorithm
// Ref: Ferri, R., Cammi, A., & Mazzei, D. (2008). Molten salt mixture properties in RELAP5 code for thermodynamic solar applications. International Journal of Thermal Sciences, 47(12), 1676–1687, pp. 1678
// Valid from 533.15K to 873.15K liquid on saturation curve
      beta_isb := 0.636 / (2263.72 - 0.636 * T);
    end beta_T;

    function cp_T "Specific heat capacity of molten salt at constant pressue as a function of temperature"
      extends Modelica.Icons.Function;
      input Modelica.SIunits.Temperature T "Temperature";
      output Modelica.SIunits.SpecificHeatCapacity cp "Specific heat capacity";
    algorithm
//Ref: Zavoico, A. B. (2001). Solar Power Tower - Design Basis Document. Technical Report SAND2001-2100. Alburquerque, USA, pp. 23
// Valid from 533.15K to 873.15K liquid on saturation curve:
      cp := 1396.0182 + 0.172 * T;
    end cp_T;

    function dh_dT_T "Derivative of specific enthalpy of molten salt w.r.t temperature"
      import SolarTherm.Media.Sodium.Sodium_utilities.*;
      extends Modelica.Icons.Function;
      input Modelica.SIunits.Temperature T "Temperature";
      output Real dh "Derivative of specific enthalpy w.r.t temperature";
    algorithm
      dh := cp_T(T);
    end dh_dT_T;

    function dp_sat_dT_T "Derivative of the saturation pressure of molten salt w.r.t temperature"
      import SolarTherm.Media.MoltenSalt.MoltenSalt_utilities.*;
      extends Modelica.Icons.Function;
      input Modelica.SIunits.Temperature T "Temperature";
      output Real dp_sat "Derivative of the saturation pressure w.r.t temperature";
    algorithm
// Ref: Ferri, R., Cammi, A., & Mazzei, D. (2008). Molten salt mixture properties in RELAP5 code for thermodynamic solar applications. International Journal of Thermal Sciences, 47(12), 1676–1687, pp. 1679
      dp_sat := (-5523.9586 / T ^ 2) * p_sat_T(T);
    end dp_sat_dT_T;

    function drho_dT_T "Derivative of density of molten salt w.r.t temperature"
      extends Modelica.Icons.Function;
      input Modelica.SIunits.Temperature T "Temperature";
      output Real drho "Derivative of density w.r.t temperature";
    algorithm
// Ref: Zavoico, A. B. (2001). Solar Power Tower - Design Basis Document. Technical Report SAND2001-2100. Alburquerque, USA, pp. 23
// Valid from 533.15K to 873.15K liquid on saturation curve:
      drho := -0.636;
    end drho_dT_T;    

    function dT_dh_h "Derivative of temperature of molten salt w.r.t specific enthalpy"
      extends Modelica.Icons.Function;
      import SolarTherm.Media.MoltenSalt.MoltenSalt_utilities.*;
      input Modelica.SIunits.SpecificEnthalpy h "Specific enthalpy";
      output Real dT "Derivative of temperature w.r.t specific enthalpy";
    protected
      Modelica.SIunits.Temperature T;
    algorithm
// Ref: Ferri, R., Cammi, A., & Mazzei, D. (2008). Molten salt mixture properties in RELAP5 code for thermodynamic solar applications. International Journal of Thermal Sciences, 47(12), 1676–1687, pp. 1678
// Valid from 533.15K to 873.15K liquid on saturation curve
      T := T_h(h);
      dT := 1 / dh_dT_T(T);
    end dT_dh_h;

    function dT_sat_dp_p "Derivative of the saturation temperature of molten salt w.r.t pressure"
      import SolarTherm.Media.MoltenSalt.MoltenSalt_utilities.*;
      extends Modelica.Icons.Function;
      input Modelica.Media.Interfaces.Types.AbsolutePressure p "Pressure";
      output Real dT_sat "Derivative of the saturation temperature w.r.t pressure";
    algorithm
// Ref: Ferri, R., Cammi, A., & Mazzei, D. (2008). Molten salt mixture properties in RELAP5 code for thermodynamic solar applications. International Journal of Thermal Sciences, 47(12), 1676–1687, pp. 1679
      dT_sat := T_sat_p(p) ^ 2 / (5523.9586 * p);
    end dT_sat_dp_p;

    function eta_T "Dynamic viscosity of molten salt as a function of temperature"
      extends Modelica.Icons.Function;
      input Modelica.SIunits.Temperature K "Temperature";
      output Modelica.SIunits.DynamicViscosity eta "Dynamic viscosity";
    protected
      Modelica.SIunits.Temp_C T = Modelica.SIunits.Conversions.to_degC(K) "Temperature in degCelsius";
    algorithm
// Ref: Zavoico, A. B. (2001). Solar Power Tower - Design Basis Document. Technical Report SAND2001-2100. Alburquerque, USA, pp. 23
// Valid from 533.15K to 873.15K liquid on saturation curve:
      eta := 0.001 * (22.714 - 0.120 * T + 2.281e-4 * T ^ 2 - 1.474e-7 * T ^ 3);
    end eta_T;

    function gibbs_T "Specific Gibbs energy of molten salt"
      import SolarTherm.Media.MoltenSalt.MoltenSalt_utilities.*;
      extends Modelica.Icons.Function;
      input Modelica.SIunits.Temperature T "Temperature";
      output Modelica.SIunits.SpecificGibbsFreeEnergy gibbs "Specific Gibbs energy";
    algorithm
      gibbs := h_T(T) - T * s_T(T);
    end gibbs_T;

    function h_fg_p "Specific enthalpy of vaporization of molten salt as a function of pressure"
      extends Modelica.Icons.Function;
      input Modelica.Media.Interfaces.Types.AbsolutePressure p "Pressure";
      output Modelica.SIunits.SpecificEnthalpy h_fg "Specific enthalpy";
    algorithm
// Ref: Ferri, R., Cammi, A., & Mazzei, D. (2008). Molten salt mixture properties in RELAP5 code for thermodynamic solar applications. International Journal of Thermal Sciences, 47(12), 1676–1687, pp. 1679
      h_fg := 0.0507 * (1e7 - p);
    end h_fg_p;

    function h_rho "Specific enthalpy of molten salt as a function of density"
      extends Modelica.Icons.Function;
      input Modelica.SIunits.Density rho "Density";
      output Modelica.SIunits.SpecificEnthalpy h "Specific enthalpy";
    protected
      constant Real p1 = 39.713718745914576;
      constant Real p2 = 50.407945969307448;
      constant Real p3 = -1.228120420950819e02;
      constant Real p4 = -1.394940386838149e02;
      constant Real p5 = 1.134882498794937e02;
      constant Real p6 = 1.049937053496403e02;
      constant Real p7 = 9.752089378977765e02;
      constant Real p8 = -1.641079585719889e05;
      constant Real p9 = 1.019046845737454e06;
      constant Real rho_mean = 1817;
      constant Real rho_std = 69.27;
      Real x;
    algorithm
// Valid from 533.15K to 873.15K
      x := (rho - rho_mean) / rho_std;
      h := p1 * x ^ 8 + p2 * x ^ 7 + p3 * x ^ 6 + p4 * x ^ 5 + p5 * x ^ 4 + p6 * x ^ 3 + p7 * x ^ 2 + p8 * x + p9;
    end h_rho;

    function h_s "Specific enthalpy of molten salt as a function of Specific entropy"
      import SolarTherm.Media.MoltenSalt.MoltenSalt_utilities.*;
      extends Modelica.Icons.Function;
      input Modelica.SIunits.SpecificEntropy s "Specific entropy";
      output Modelica.SIunits.SpecificEnthalpy h "Specific enthalpy";
    protected
      constant Real p1 = -4.207915256284852e-06;
      constant Real p2 = -4.252960376115911e-04;
      constant Real p3 = -0.007372023879740;
      constant Real p4 = 0.235040918858859;
      constant Real p5 = 18.671039765272749;
      constant Real p6 = 6.228064419983780e02;
      constant Real p7 = 1.292193798421334e04;
      constant Real p8 = 1.646825778116541e05;
      constant Real p9 = 1.008137894767280e06;
      constant Real s_mean = 9250;
      constant Real s_std = 237.8;
      Real x;
    algorithm
// Valid from 533.15K to 873.15K
      x := (s - s_mean) / s_std;
      h := p1 * x ^ 8 + p2 * x ^ 7 + p3 * x ^ 6 + p4 * x ^ 5 + p5 * x ^ 4 + p6 * x ^ 3 + p7 * x ^ 2 + p8 * x + p9;
    end h_s;

    function h_T "Specific enthalpy of molten salt as a function of temperature"
      extends Modelica.Icons.Function;
      input Modelica.SIunits.Temperature T "Temperature";
      output Modelica.SIunits.SpecificEnthalpy h "Specific enthalpy";
    algorithm
// h is obtained by integrating (cp dT). The integration constant was added such that the h value at T = 0 K becomes zero.
// Ref: Ferri, R., Cammi, A., & Mazzei, D. (2008). Molten salt mixture properties in RELAP5 code for thermodynamic solar applications. International Journal of Thermal Sciences, 47(12), 1676–1687, pp. 1678
// Valid from 533.15K to 873.15K liquid on saturation curve
      h := 1396.0182 * T + 0.086 * T ^ 2;
      annotation(
        derivative = h_T_der);
    end h_T;

    function h_T_der "Derivative of specific enthalpy of molten salt w.r.t. time"
      import SolarTherm.Media.MoltenSalt.MoltenSalt_utilities.*;
      extends Modelica.Icons.Function;
      input Modelica.SIunits.Temperature T "Temperature";
      input Real der_T;
      output Real der_h "Time derivative of specific enthalpy";
    algorithm
      der_h := dh_dT_T(T) * der_T;
    end h_T_der;

    function helmholtz_pT "Specific Helmholtz energy of molten salt"
      import SolarTherm.Media.MoltenSalt.MoltenSalt_utilities.*;
      extends Modelica.Icons.Function;
      input Modelica.SIunits.Temperature T "Temperature";
      input Modelica.SIunits.AbsolutePressure p "Pressure";
      output Modelica.SIunits.SpecificHelmholtzFreeEnergy helmholtz "Specific Helmholtz energy";
    algorithm
      helmholtz := h_T(T) - p / rho_T(T) - T * s_T(T);
    end helmholtz_pT;

    function kappa_T "Isothermal compressibility of molten salt"
      extends Modelica.Icons.Function;
      input Modelica.SIunits.Temperature T "Temperature";
      output Modelica.SIunits.IsothermalCompressibility kappa_ist "Isothermal compressibility";
    algorithm
// Ref: Ferri, R., Cammi, A., & Mazzei, D. (2008). Molten salt mixture properties in RELAP5 code for thermodynamic solar applications. International Journal of Thermal Sciences, 47(12), 1676–1687, pp. 1678
// Valid from 533.15K to 873.15K liquid on saturation curve
// A correlation of specific volume versus pressure was not available. Instead, the following correlation as a function of temperature typical of liquid metals was found in the literature:
      kappa_ist := 3.97e-11 * exp(1.37e-3 * (T - 293.15));
    end kappa_T;

    function lamda_T "Thermal conductivity of molten salt as a function of temperature"
      extends Modelica.Icons.Function;
      input Modelica.SIunits.Temperature K "Temperature";
      output Modelica.SIunits.ThermalConductivity lamda "Thermal conductivity";
    protected
      Modelica.SIunits.Temp_C T = Modelica.SIunits.Conversions.to_degC(K) "Temperature in degCelsius";
    algorithm
// Ref: Zavoico, A. B. (2001). Solar Power Tower - Design Basis Document. Technical Report SAND2001-2100. Alburquerque, USA, pp. 23
// Valid from 533.15K to 873.15K liquid on saturation curve:
      lamda := 0.443 + 1.9e-4 * T;
    end lamda_T;

    function p_rho "Saturated vapour pressur of molten salt as a function of density"
      extends Modelica.Icons.Function;
      input Modelica.SIunits.Density rho "Density";
      output Modelica.Media.Interfaces.Types.AbsolutePressure p "Pressure";
    protected
      constant Real p1 = 0.733446745332821;
      constant Real p2 = 0.653577028271395;
      constant Real p3 = -1.254638858135562;
      constant Real p4 = 23.599040313108496;
      constant Real p5 = 16.841697526950028;
      constant Real p6 = -1.832200683871549e03;
      constant Real p7 = 1.001166286167321e04;
      constant Real p8 = -2.197288709140443e04;
      constant Real p9 = 1.800126936525206e04;
      constant Real rho_mean = 1817;
      constant Real rho_std = 69.27;
      Real x;
    algorithm
// Valid from 533.15K to 873.15K
      x := (rho - rho_mean) / rho_std;
      p := p1 * x ^ 8 + p2 * x ^ 7 + p3 * x ^ 6 + p4 * x ^ 5 + p5 * x ^ 4 + p6 * x ^ 3 + p7 * x ^ 2 + p8 * x + p9;
    end p_rho;

    function p_sat_T "Saturation pressure of molten salt as a function of temperature"
      extends Modelica.Icons.Function;
      input Modelica.SIunits.Temperature T "Temperature";
      output Modelica.Media.Interfaces.Types.AbsolutePressure p "Pressure";
    algorithm
// Ref: Ferri, R., Cammi, A., & Mazzei, D. (2008). Molten salt mixture properties in RELAP5 code for thermodynamic solar applications. International Journal of Thermal Sciences, 47(12), 1676–1687, pp. 1679
      p := exp(17.69185 - 5523.9586 / T);
      annotation(
        derivative = p_sat_T_der);
    end p_sat_T;

    function p_sat_T_der "Derivative of the saturation pressure of molten salt w.r.t. time"
      import SolarTherm.Media.MoltenSalt.MoltenSalt_utilities.*;
      extends Modelica.Icons.Function;
      input Modelica.SIunits.Temperature T "Temperature";
      input Real der_T "Derivative of T w.r.t. time";
      output Real der_p_sat "Derivative of saturation pressure w.r.t time";
    algorithm
      der_p_sat := dp_sat_dT_T(T) * der_T;
    end p_sat_T_der;

    function rho_T "Density of molten salt as a function of temperature"
      extends Modelica.Icons.Function;
      input Modelica.SIunits.Temperature K "Temperature";
      output Modelica.SIunits.Density rho "Density";
    protected
      Modelica.SIunits.Temp_C T = Modelica.SIunits.Conversions.to_degC(K) "Temperature in degCelsius";
    algorithm
// Ref: Zavoico, A. B. (2001). Solar Power Tower - Design Basis Document. Technical Report SAND2001-2100. Alburquerque, USA, pp. 23
// Valid from 533.15K to 873.15K liquid on saturation curve:
      rho := 2090 - 0.636 * T;
      annotation(
        derivative = rho_T_der);
    end rho_T;

    function rho_T_der "Derivative of the density of molten salt w.r.t. time"
      import SolarTherm.Media.MoltenSalt.MoltenSalt_utilities.*;
      extends Modelica.Icons.Function;
      input Modelica.SIunits.Temperature T "Temperature";
      input Real der_T "Derivative of temperature w.r.t time";
      output Real der_rho "Derivative of density w.r.t time";
    algorithm
      der_rho := drho_dT_T(T) * der_T;
    end rho_T_der;

    function s_rho "Specific entropy of molten salt as a function of density"
      extends Modelica.Icons.Function;
      input Modelica.SIunits.Density rho "Density";
      output Modelica.SIunits.SpecificEntropy s "Specific entropy";
    protected
      constant Real p1 = 0.074442895162941;
      constant Real p2 = 0.094105565955424;
      constant Real p3 = -0.233574614812304;
      constant Real p4 = -0.286315713014475;
      constant Real p5 = 0.013365595797011;
      constant Real p6 = -1.523616571400823;
      constant Real p7 = -16.749820096832305;
      constant Real p8 = -2.345160845000009e02;
      constant Real p9 = 9.265390905220705e03;
      constant Real rho_mean = 1817;
      constant Real rho_std = 69.27;
      Real x;
    algorithm
// 400K to 2500K liquid on saturation curve:
      x := (rho - rho_mean) / rho_std;
//rho_norm
      s := p1 * x ^ 8 + p2 * x ^ 7 + p3 * x ^ 6 + p4 * x ^ 5 + p5 * x ^ 4 + p6 * x ^ 3 + p7 * x ^ 2 + p8 * x + p9;
    end s_rho;

    function s_T "Specific entropy of molten salt as a function of temperature"
      extends Modelica.Icons.Function;
      input Modelica.SIunits.Temperature T "Temperature";
      output Modelica.SIunits.SpecificEntropy s "Specific entropy";
    algorithm
// Ref: Ferri, R., Cammi, A., & Mazzei, D. (2008). Molten salt mixture properties in RELAP5 code for thermodynamic solar applications. International Journal of Thermal Sciences, 47(12), 1676–1687, pp. 1678
// Valid from 533.15K to 873.15K liquid on saturation curve
      s := 1396.0182 * log(T) + 0.172 * T;
    end s_T;

    function T_h "Temperature of molten salt as a function of specific enthalpy"
      extends Modelica.Icons.Function;
      input Modelica.SIunits.SpecificEnthalpy h "Specific enthalpy";
      output Modelica.SIunits.Temperature T "Temperature";
    protected
      Real delta;
    algorithm
// Ref: Ferri, R., Cammi, A., & Mazzei, D. (2008). Molten salt mixture properties in RELAP5 code for thermodynamic solar applications. International Journal of Thermal Sciences, 47(12), 1676–1687, pp. 1678
// Valid from 533.15K to 873.15K liquid on saturation curve
      delta := abs(1396.0182 ^ 2 + 4 * 0.086 * h);
      T := ((-1396.0182) + sqrt(delta)) / (2 * 0.086);
      annotation(
        derivative = T_h_der);
    end T_h;

    function T_h_der "Derivative of temperature of molten salt w.r.t. time"
      import SolarTherm.Media.MoltenSalt.MoltenSalt_utilities.*;
      extends Modelica.Icons.Function;
      input Modelica.SIunits.SpecificEnthalpy h "Specific enthalpy";
      input Real der_h "Derivative of specific enthalpy w.r.t. time";
      output Real der_T "Derivative of Temperature w.r.t. time";
    algorithm
// Ref: Ferri, R., Cammi, A., & Mazzei, D. (2008). Molten salt mixture properties in RELAP5 code for thermodynamic solar applications. International Journal of Thermal Sciences, 47(12), 1676–1687, pp. 1678
// Valid from 533.15K to 873.15K liquid on saturation curve
      der_T := dT_dh_h(h) * der_h;
    end T_h_der;

    function T_rho "Temperature of molten salt as a function of density"
      extends Modelica.Icons.Function;
      input Modelica.SIunits.Density rho "Density";
      output Modelica.SIunits.Temperature T "Temperature";
    protected
      constant Real p1 = 0.026694736413680;
      constant Real p2 = 0.033883173706630;
      constant Real p3 = -0.082551702426422;
      constant Real p4 = -0.093764993849291;
      constant Real p5 = 0.076284443065688;
      constant Real p6 = 0.070574586765291;
      constant Real p7 = -0.021388297459061;
      constant Real p8 = -1.082235011008968e02;
      constant Real p9 = 6.997982926924665e02;
      constant Real rho_mean = 1817;
      constant Real rho_std = 69.27;
      Real x;
    algorithm
// Valid from 533.15K to 873.15K
      x := (rho - rho_mean) / rho_std;
      T := p1 * x ^ 8 + p2 * x ^ 7 + p3 * x ^ 6 + p4 * x ^ 5 + p5 * x ^ 4 + p6 * x ^ 3 + p7 * x ^ 2 + p8 * x + p9;
    end T_rho;

    function T_s "Temperature of molten salt as a function of Specific entropy"
      extends Modelica.Icons.Function;
      input Modelica.SIunits.SpecificEntropy s "Specific entropy";
      output Modelica.SIunits.Temperature T "Temperature";
    protected
      constant Real p1 = 1.504093438015698e-08;
      constant Real p2 = -1.334046332876201e-07;
      constant Real p3 = -1.257773555936631e-05;
      constant Real p4 = -1.860482246373618e-04;
      constant Real p5 = 0.004942593053630;
      constant Real p6 = 0.314097262254148;
      constant Real p7 = 7.857963064919012;
      constant Real p8 = 1.086909242328156e02;
      constant Real p9 = 6.926012804626644e02;
      constant Real s_mean = 9250;
      constant Real s_std = 237.8;
      Real x;
    algorithm
// Valid from 533.15K to 873.15K
      x := (s - s_mean) / s_std;
      T := p1 * x ^ 8 + p2 * x ^ 7 + p3 * x ^ 6 + p4 * x ^ 5 + p5 * x ^ 4 + p6 * x ^ 3 + p7 * x ^ 2 + p8 * x + p9;
    end T_s;

    function T_sat_p "Saturation temperature of molten salt as a function of saturated pressure"
      extends Modelica.Icons.Function;
      input Modelica.Media.Interfaces.Types.AbsolutePressure p "Pressure";
      output Modelica.SIunits.Temperature T "Temperature";
    algorithm
// Ref: Ferri, R., Cammi, A., & Mazzei, D. (2008). Molten salt mixture properties in RELAP5 code for thermodynamic solar applications. International Journal of Thermal Sciences, 47(12), 1676–1687, pp. 1679
      T := 5523.9586 / (17.69185 - log(p));
      annotation(
        derivative = T_sat_p_der);
    end T_sat_p;

    function T_sat_p_der "Derivative of the saturation pressure of molten salt w.r.t. time"
      import SolarTherm.Media.MoltenSalt.MoltenSalt_utilities.*;
      extends Modelica.Icons.Function;
      input Modelica.Media.Interfaces.Types.AbsolutePressure p "Pressure";
      input Real der_p "Derivative of p w.r.t. time";
      output Real der_T_sat "Derivative of saturation temperature w.r.t time";
    algorithm
      der_T_sat := dT_sat_dp_p(p) * der_p;
    end T_sat_p_der;
  end MoltenSalt_utilities;
end MoltenSalt;
