within Steps.Utilities;

model Util

	function ToK
	  input Real DegC "celsious Degree";
	  output Real DegK "Degree in Kelvin";
	  
	protected
	  constant Real CONV_DEG2K = 273.15 "Temperature conversion constant degC -> K"; 
	algorithm

	  DegK := DegC + CONV_DEG2K;
	  
	end ToK;
	
	function thermal_conductivity "cal thermal conductivity of a material"
    extends Modelica.Icons.Function;
    
    input String name "material name";
    input Modelica.SIunits.Temp_C temperature;
    input Modelica.Blocks.Types.ExternalCombiTable1D tableID;
    output Modelica.SIunits.ThermalConductivity k;
    
  algorithm
    if(Modelica.Utilities.Strings.compare(name, "inconel 750") == Modelica.Utilities.Types.Compare.Equal) then
      // FIX IT: icol = 0， 1？ 
      k := Modelica.Blocks.Tables.CombiTable1D.getTableValue(tableID, icol = 1, u = temperature, tableAvailable = 1.0); 
    else
      k := 16.2;
    end if;
    
  end thermal_conductivity;   
	
	function myAssert "customerized assert function to generate more detailed information"
    input Boolean debug = false;
    input Real val_test;
    input Real min;
    input Real max;
    input String name_val;
    
    input Real [:] val_ref;
    input String [:] name_val_ref;
    
    algorithm
    
      if not debug then
        assert(val_test > min and val_test < max, "outside range " + keyvalStr(name_val, val_test) + " at " + debugInfo(name_val_ref, val_ref));
      end if;
      
  end myAssert;
  
  function sameValue "Return if two floats (almost) equal to each other"
    extends Modelica.Icons.Function;
    input Real c1 "compared number 1";
    input Real c2 "compared number 2";    
    output Boolean same;
  algorithm
      same := abs(c1 - c2) <= 1e-3;
  end sameValue;  
  
  function debugInfo "Generate debug info for output"
    extends Modelica.Icons.Function;
    input String keys[:];
    input Real values[:];
    output String info;
    
    protected
    
      String buff;
    
    algorithm
      buff := ""; // initialization required. 
      
      for i in 1:size(keys, 1) loop
        buff := buff + keyvalStr(keys[i],values[i]); 
        if i <> size(keys, 1) then
          buff := buff + ",";
        end if;
      end for;
      
      info := buff;
    
  end debugInfo;
  
  function keyvalStr "generate a key=value string"
    extends Modelica.Icons.Function;
    
    input String key;
    input Real value;
    output String out;
    
    algorithm
      out := key + " = " + String(value); // intrim " = " required  
  end keyvalStr;  
  
  function E2I "Convert the enumeration pipetype into Integer as an index"
    input Steps.Components.PipeType enum;
    output Integer index;
    
    algorithm
      index := Integer(enum);
  end E2I;  
end Util;
