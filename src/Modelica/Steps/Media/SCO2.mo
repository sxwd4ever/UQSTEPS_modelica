within Steps.Media;
package SCO2 "supercritical CO2"
  extends ExternalMedia.Media.CoolPropMedium(
      // mediumName = "CarbonDioxide",
      mediumName = "CO2",
      // substanceNames = {"CO2|debug=40"}, // for single test, more detailed output
      substanceNames = {"CO2"}, // for parameters sweep, lesser output with the debug flag
      ThermoStates = Modelica.Media.Interfaces.Choices.IndependentVariables.pT,
      // inputChoice = ExternalMedia.Common.InputChoice.pT,
      singleState=false,
      onePhase = true,
      final reducedX = true, 
      final fixedX = true, 
      // just for Pijarra Hill Experimental data which is below critical point. 
      // Temperature(min=204.128, max=2000, start=204.128, nominal=204.128),
      // Valid range of temperature
      Temperature(min=304.128, max=2000, start=304.128, nominal=304.128),
      //SpecificEnthalpy(start=2e5, nominal=1.0e5),
      AbsolutePressure(min=7.377e6, start=7.377e6, nominal=7.377e6),      
      SpecificEnthalpy(start=if Functions.referenceChoice==ReferenceEnthalpy.ZeroAt0K then data.H0 else
        if Functions.referenceChoice==ReferenceEnthalpy.UserDefined then Functions.h_offset else 0, nominal=1.0e5),      
     Density(start=10, nominal=10),
     fluidConstants = {CO2FluidConstants});
     
  import Modelica.Media.Interfaces.Choices.ReferenceEnthalpy;   
  import Modelica.Media.IdealGases.Common;
  import Modelica.Media.IdealGases.Common.Functions;
  import Modelica.Media.Air.ReferenceAir.Air_Utilities;
  import SI = Modelica.SIunits;
              
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
  /*
    type AbsolutePressure    
      "Type for absolute pressure with medium specific attributes"
    extends SI.AbsolutePressure(
      min=7.377e6,
      max=1.0e8,
      nominal=7.377e6,
      start=7.377e6);
    end AbsolutePressure;
    
    type SpecificEnthalpy 
    extends SI.SpecificEnthalpy(
      start=if Functions.referenceChoice==ReferenceEnthalpy.ZeroAt0K then data.H0 else
        if Functions.referenceChoice==ReferenceEnthalpy.UserDefined then Functions.h_offset else 0, 
      nominal=1.0e6);
    end SpecificEnthalpy;
      
    type Temperature
      "Type for temperature with medium specific attributes"
      extends SI.Temperature(
      min=304.128,
      max=1.0e4,
      nominal=304.128,
      start=304.128);
    end Temperature;
  */          
    constant Modelica.Media.Interfaces.Types.IdealGas.FluidConstants CO2Constants = Modelica.Media.IdealGases.Common.FluidData.CO2;
    
    replaceable constant FluidConstants CO2FluidConstants = FluidConstants(
      iupacName=  CO2Constants.iupacName,
      casRegistryNumber=  CO2Constants.casRegistryNumber,
      chemicalFormula=  CO2Constants.chemicalFormula,
      structureFormula=  "unknown",
      molarMass=  0.0440098, // "value in CoolProp slightly different from CO2Constants.molarMass 0.0440095",
      criticalTemperature= 304.128, // "value in CoolProp slightly different from CO2Constants.molarMass 304.12",
      criticalPressure= 7.3773e6, // "value in CoolProp slightly different from CO2Constants.molarMass 73.74e5",
      criticalMolarVolume= 9.41185e-5, // "value in CoolProp slightly different from CO2Constants.molarMass 94.07e-6",
      acentricFactor=  CO2Constants.acentricFactor,
      triplePointTemperature=  280,
      triplePointPressure=  500,
      meltingPoint=  CO2Constants.meltingPoint,
      normalBoilingPoint=  CO2Constants.normalBoilingPoint,
      dipoleMoment=  CO2Constants.dipoleMoment); 	  

  // It is necessary to redeclare this function 
  // althought these two (here and ExternalMedia's code) are identical. 
  // Or else runtime error will occour in the way that Hmoler, P, T will out of range
  // The cause of this bug is not clear yet. - Xin, 06/05/2021

  redeclare replaceable function setState_pT
    "Return thermodynamic state record from p and T"
    extends Modelica.Icons.Function;
    input AbsolutePressure p "pressure";
    input Temperature T "temperature";
    input FixedPhase phase = 1
      "2 for two-phase, 1 for one-phase, 0 if not known";
    output ThermodynamicState state;
  algorithm  
    // assert((p > 1000 and p < 2.5e8) and (h > 0 and h < 2e8), "Invalid Pressure or Specific Enthalpy (p, h)=(" + String(p) + "," + String(h) + ")");  
    state := setState_pT_lib(p, T, phase);
    
    // external media indicates the computation error by assign negetive values for p, h, T.
    if(state.p < 0 and state.T < 0 and state.h < 0) then
      assert(false, "Error in region computation of ExternalMedia" + "(p = " + String(p) + ", T = " + String(T) + ")");    
    end if;      
    
    
    annotation(Include="#include \"externalmedialib.h\"", Library="ExternalMediaLib", IncludeDirectory="modelica://ExternalMedia/Resources/Include", LibraryDirectory="modelica://ExternalMedia/Resources/Library");
  end setState_pT;


  function setState_pT_lib
    "Return thermodynamic state record from p and T"
    extends Modelica.Icons.Function;
    input AbsolutePressure p "pressure";
    input Temperature T "temperature";
    input FixedPhase phase = 1
      "2 for two-phase, 1 for one-phase, 0 if not known";
    output ThermodynamicState state;
  external "C" TwoPhaseMedium_setState_pT_C_impl(p, T, state, mediumName, libraryName, substanceName)
    annotation(Include="#include \"externalmedialib.h\"", Library="ExternalMediaLib", IncludeDirectory="modelica://ExternalMedia/Resources/Include", LibraryDirectory="modelica://ExternalMedia/Resources/Library");
  end setState_pT_lib;
  
  
/*  
  redeclare replaceable function setState_ph
    input AbsolutePressure p "pressure";
    input SpecificEnthalpy h "specific enthalpy";
    input FixedPhase phase = 1
      "2 for two-phase, 1 for one-phase, 0 if not known";
    output ThermodynamicState state;      
  external "C" TwoPhaseMedium_setState_ph_C_impl(p, h, phase, state, mediumName, libraryName, substanceName)
    annotation(Include="#include \"externalmedialib.h\"", Library="ExternalMediaLib", IncludeDirectory="modelica://ExternalMedia/Resources/Include", LibraryDirectory="modelica://ExternalMedia/Resources/Library");
  end setState_ph; 
*/    

  // To locate the error source during initialization, use this wrapper function to add assert()
  redeclare replaceable function setState_ph
    input AbsolutePressure p "pressure";
    input SpecificEnthalpy h "specific enthalpy";
    input FixedPhase phase = 1
      "2 for two-phase, 1 for one-phase, 0 if not known";
    output ThermodynamicState state;        
   
  algorithm    
  
    // assert((p > 1000 and p < 2.5e8) and (h > 0 and h < 2e8), "Invalid Pressure or Specific Enthalpy (p, h)=(" + String(p) + "," + String(h) + ")");  
    state := setState_ph_lib(p, h, phase);
    
    // external media indicates the computation error by assign negetive values for p, h, T.
    if(state.p < 0 and state.T < 0 and state.h < 0) then
      assert(false, "Error in region computation of ExternalMedia" + "(p = " + String(p) + ", h = " + String(h) + ")");    
    end if;  
    
  end setState_ph;
    
  function setState_ph_lib
    "Return thermodynamic state record from p and h"
    extends Modelica.Icons.Function;
    input AbsolutePressure p "pressure";
    input SpecificEnthalpy h "specific enthalpy";
    input FixedPhase phase = 1
      "2 for two-phase, 1 for one-phase, 0 if not known";
    output ThermodynamicState state;
  external "C" TwoPhaseMedium_setState_ph_C_impl(p, h, phase, state, mediumName, libraryName, substanceName)
    annotation(Include="#include \"externalmedialib.h\"", Library="ExternalMediaLib", IncludeDirectory="modelica://ExternalMedia/Resources/Include", LibraryDirectory="modelica://ExternalMedia/Resources/Library");
  end setState_ph_lib;  

/*  
  redeclare function extends specificEnthalpy "Return specific enthalpy as a function of the thermodynamic state record"
      algorithm
        h := specificEnthalpy_pT(p = state.p, T = state.T);
      annotation(
        Inline = true,
        smoothOrder = 2);
    end specificEnthalpy;  
*/
  /*
  // To locate the error source during initialization, use this wrapper function to add assert()
  redeclare replaceable function setState_pT
    "Return thermodynamic state record from p and T"
    extends Modelica.Icons.Function;
    input AbsolutePressure p "pressure";
    input Temperature T "temperature";
    input FixedPhase phase = 1
      "2 for two-phase, 1 for one-phase, 0 if not known";
    output ThermodynamicState state;
  algorithm
    assert(p >= 7.3773e6, "pressure out of range");
    assert(T >= 304.128 and T <=2000, "temperature out of critical range");
    state := setState_pT_lib(p, T, phase);
  end setState_pT;


  function setState_pT_lib
    "Return thermodynamic state record from p and T"
    extends Modelica.Icons.Function;
    input AbsolutePressure p "pressure";
    input Temperature T "temperature";
    input FixedPhase phase = 1
      "2 for two-phase, 1 for one-phase, 0 if not known";
    output ThermodynamicState state;
  external "C" TwoPhaseMedium_setState_pT_C_impl(p, T, state, mediumName, libraryName, substanceName)
    annotation(Include="#include \"externalmedialib.h\"", Library="ExternalMediaLib", IncludeDirectory="modelica://ExternalMedia/Resources/Include", LibraryDirectory="modelica://ExternalMedia/Resources/Library");
  end setState_pT_lib;
  */   
end SCO2;
