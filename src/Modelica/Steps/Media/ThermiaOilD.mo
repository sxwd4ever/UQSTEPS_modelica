within Steps.Media;

package ThermiaOilD "Shell Thermia Oil D the heat transfer fluid"
  
  import Modelica.SIunits.Conversions.from_degC;
  import Poly = Modelica.Media.Incompressible.TableBased.Polynomials_Temp;
  
  extends Incompressible.TableBased(
    mediumName="Thermia Oil D",
    T_min = from_degC(0), T_max = from_degC(340),
    TinK = false, T0=273.15,
    tableDensity=
      [0, 894; 20, 882; 40, 870; 100, 834; 150, 803; 200, 773; 250, 743; 300, 712; 340, 688
],
    tableHeatCapacity=
      [0, 1925; 20, 1925; 40, 2007; 100, 2254; 150, 2459; 200, 2665; 250, 2871; 300, 3076; 340, 3241],
    tableConductivity=
      [0, 0.169; 20, 0.165; 40, 0.162; 100, 0.151; 150, 0.142; 200, 0.134; 250, 0.125; 300, 0.116; 340, 0.109
],
    tableViscosity = [0, 1.544; 20, 0.283; 40, 0.085; 100, 0.009; 150, 0.003; 200, 0.002; 250, 0.001; 300, 0.001; 340, 0.001],
    hasVaporPressure = false
    );
  /*  
  redeclare replaceable record ThermodynamicState "A selection of variables that uniquely defines the thermodynamic state"
    extends ExternalMedia.Media.BaseClasses.ExternalTwoPhaseMedium.ThermodynamicState;
      ExternalMedia.Media.BaseClasses.ExternalTwoPhaseMedium.AbsolutePressure p "Absolute pressure of medium";
      ExternalMedia.Media.BaseClasses.ExternalTwoPhaseMedium.Temperature T "Temperature of medium";
  end ThermodynamicState;
*/
  /*
    function setState_ph "Returns state record as function of p and h"
      extends Modelica.Icons.Function;
      input AbsolutePressure p "Pressure";
      input SpecificEnthalpy h "Specific enthalpy";
      input FixedPhase phase=0
      "2 for two-phase, 1 for one-phase, 0 if not known";      
      output ThermodynamicState state "Thermodynamic state";
    algorithm
      state :=ThermodynamicState(p=p,T=T_ph(p,h));
      annotation(Inline=true,smoothOrder=3);
    end setState_ph;
    */
  function setState_ph
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
  
  redeclare replaceable function extends density_derp_T "Density derivative by pressure"
  algorithm
    ddpT := 0 "Incompressible"; 
    // CP.PropsSI("d(DMASS)/d(P)|Hmass", "P", state.p, "Hmass", specificEnthalpy_pT(state.p, state.T), mediumName);
    annotation (Inline=true);
  end density_derp_T; 
  /* 
  redeclare replaceable function extends density_derh_p
      "Return density derivative w.r.t. specific enthalpy at constant pressure"
      algorithm
        ddhp := Poly.derivativeValue(poly_rho, h_T(state.T));
        //ddhp:= drho_dT_T(state.T) * dT_dh_h(h_T(state.T));
    end density_derh_p;    
  */
    annotation (Documentation(info="<html>

</html>"));
end ThermiaOilD;
