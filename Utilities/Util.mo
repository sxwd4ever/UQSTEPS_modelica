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
end Util;
