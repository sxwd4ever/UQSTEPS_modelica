within Steps.Components;

model MaterialConductivity  "Material Conductivity table"
  import MyUtil = Steps.Utilities.Util;
  import TB = Modelica.Blocks.Tables; 
  
  parameter String name_material = "inconel_750";
  
  //input Modelica.SIunits.Temperature T;
  
  //output Modelica.SIunits.ThermalConductivity k;  
  
  //Modelica.Blocks.Types.ExternalCombiTable1D table_conductivity;
  
//protected 	
  Modelica.Blocks.Types.ExternalCombiTable1D table_th_inconel_750 = Modelica.Blocks.Types.ExternalCombiTable1D(tableName = "inconel_750", fileName = Modelica.Utilities.Files.loadResource("modelica://Steps/Resources/Data/th_conductivity.txt"), table = fill(0.0, 6, 2), smoothness = Modelica.Blocks.Types.Smoothness.LinearSegments, columns = 2:2) "thermal conductivity for inconel_750";	

algorithm

  //table_conductivity := table_th_inconel_750;
  //MyUtil.thermal_conductivity(tableID = table_th_inconel_750, name = name_material, temperature = T);    

end MaterialConductivity;
