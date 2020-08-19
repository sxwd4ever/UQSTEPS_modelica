within Steps.Components;

model KimCorrelations  "Correlation constants in Kim [2012]"
  import MyUtil = Steps.Utilities.Util;
  import TB = Modelica.Blocks.Tables; 
  
  input Modelica.SIunits.Length pitch = 24.6 * 1e-3;
  input Modelica.SIunits.Diameter d_h = 0.922 * 1e-3;  
  input Modelica.SIunits.Angle phi = 0.0 "unit rad";
  
	output Real a;
	output Real b;
	output Real c;	
	output Real d;
	
	protected
    Modelica.Blocks.Types.ExternalCombiTable1D table_4a_a = Modelica.Blocks.Types.ExternalCombiTable1D(tableName = "4a_a", fileName = Modelica.Utilities.Files.loadResource("modelica://Steps/Resources/Data/kim_2012.txt"), table = fill(0.0, 9, 2), smoothness = Modelica.Blocks.Types.Smoothness.LinearSegments, columns = 2:2) "Table 4a - column a in Kim[2012] for pitch=24.6, dh=0.922 (dc=1.3 mm))";

  Modelica.Blocks.Types.ExternalCombiTable1D table_4a_b = Modelica.Blocks.Types.ExternalCombiTable1D(tableName = "4a_b", fileName = Modelica.Utilities.Files.loadResource("modelica://Steps/Resources/Data/kim_2012.txt"), table = fill(0.0, 9, 2), smoothness = Modelica.Blocks.Types.Smoothness.LinearSegments, columns = 2:2) "Table 4a - column b in Kim[2012] for pitch=12.3, dh=0.922 (dc=1.3 mm))";
  
  Modelica.Blocks.Types.ExternalCombiTable1D table_4b_a = Modelica.Blocks.Types.ExternalCombiTable1D(tableName = "4b_a", fileName = Modelica.Utilities.Files.loadResource("modelica://Steps/Resources/Data/kim_2012.txt"), table = fill(0.0, 9, 2), smoothness = Modelica.Blocks.Types.Smoothness.LinearSegments, columns = 2:2) "Table 4b - column a in Kim[2012] for pitch=24.6, dh=1.222 (dc=1.3 mm))";  
  
  Modelica.Blocks.Types.ExternalCombiTable1D table_4b_b = Modelica.Blocks.Types.ExternalCombiTable1D(tableName = "4b_b", fileName = Modelica.Utilities.Files.loadResource("modelica://Steps/Resources/Data/kim_2012.txt"), table = fill(0.0, 9, 2), smoothness = Modelica.Blocks.Types.Smoothness.LinearSegments, columns = 2:2) "Table 4b - column b default table";  
  
  Modelica.Blocks.Types.ExternalCombiTable1D table_4c_a = Modelica.Blocks.Types.ExternalCombiTable1D(tableName = "4c_a", fileName = Modelica.Utilities.Files.loadResource("modelica://Steps/Resources/Data/kim_2012.txt"), table = fill(0.0, 9, 2), smoothness = Modelica.Blocks.Types.Smoothness.LinearSegments, columns = 2:2) "Table 4c - column a in Kim[2012] for pitch=24.6, dh=0.922 (dc=1.3 mm))";

  Modelica.Blocks.Types.ExternalCombiTable1D table_4c_b = Modelica.Blocks.Types.ExternalCombiTable1D(tableName = "4c_b", fileName = Modelica.Utilities.Files.loadResource("modelica://Steps/Resources/Data/kim_2012.txt"), table = fill(0.0, 9, 2), smoothness = Modelica.Blocks.Types.Smoothness.LinearSegments, columns = 2:2) "Table 4c - column b in Kim[2012] for pitch=12.3, dh=0.922 (dc=1.3 mm))";  

  Modelica.Blocks.Types.ExternalCombiTable1D table_5a_c = Modelica.Blocks.Types.ExternalCombiTable1D(tableName = "5a_c", fileName = Modelica.Utilities.Files.loadResource("modelica://Steps/Resources/Data/kim_2012.txt"), table = fill(0.0, 9, 2), smoothness = Modelica.Blocks.Types.Smoothness.LinearSegments, columns = 2:2) "Table 5a - column c in Kim[2012] for pitch=24.6, dh=0.922 (dc=1.3 mm))";

  Modelica.Blocks.Types.ExternalCombiTable1D table_5a_d = Modelica.Blocks.Types.ExternalCombiTable1D(tableName = "5a_d", fileName = Modelica.Utilities.Files.loadResource("modelica://Steps/Resources/Data/kim_2012.txt"), table = fill(0.0, 9, 2), smoothness = Modelica.Blocks.Types.Smoothness.LinearSegments, columns = 2:2) "Table 5a - column d in Kim[2012] for pitch=12.3, dh=0.922 (dc=1.3 mm))";
  
  Modelica.Blocks.Types.ExternalCombiTable1D table_5b_c = Modelica.Blocks.Types.ExternalCombiTable1D(tableName = "5b_c", fileName = Modelica.Utilities.Files.loadResource("modelica://Steps/Resources/Data/kim_2012.txt"), table = fill(0.0, 9, 2), smoothness = Modelica.Blocks.Types.Smoothness.LinearSegments, columns = 2:2) "Table 5b - column c in Kim[2012] for pitch=24.6, dh=1.222 (dc=1.3 mm))";
  
  Modelica.Blocks.Types.ExternalCombiTable1D table_5b_d = Modelica.Blocks.Types.ExternalCombiTable1D(tableName = "5b_d", fileName = Modelica.Utilities.Files.loadResource("modelica://Steps/Resources/Data/kim_2012.txt"), table = fill(0.0, 9, 2), smoothness = Modelica.Blocks.Types.Smoothness.LinearSegments, columns = 2:2) "Table 5b - column d default table";
  
  Modelica.Blocks.Types.ExternalCombiTable1D table_5c_c = Modelica.Blocks.Types.ExternalCombiTable1D(tableName = "5c_c", fileName = Modelica.Utilities.Files.loadResource("modelica://Steps/Resources/Data/kim_2012.txt"), table = fill(0.0, 9, 2), smoothness = Modelica.Blocks.Types.Smoothness.LinearSegments, columns = 2:2) "Table 5c - column c in Kim[2012] for pitch=24.6, dh=0.922 (dc=1.3 mm))";

  Modelica.Blocks.Types.ExternalCombiTable1D table_5c_d = Modelica.Blocks.Types.ExternalCombiTable1D(tableName = "5c_d", fileName = Modelica.Utilities.Files.loadResource("modelica://Steps/Resources/Data/kim_2012.txt"), table = fill(0.0, 9, 2), smoothness = Modelica.Blocks.Types.Smoothness.LinearSegments, columns = 2:2) "Table 5c - column d in Kim[2012] for pitch=12.3, dh=0.922 (dc=1.3 mm))";
    
  Modelica.Blocks.Types.ExternalCombiTable1D table_a = Modelica.Blocks.Types.ExternalCombiTable1D(tableName = "4d_a", fileName = Modelica.Utilities.Files.loadResource("modelica://Steps/Resources/Data/kim_2012.txt"), table = fill(0.0, 9, 2), smoothness = Modelica.Blocks.Types.Smoothness.LinearSegments, columns = 2:2);
  
  Modelica.Blocks.Types.ExternalCombiTable1D table_b = Modelica.Blocks.Types.ExternalCombiTable1D(tableName = "4d_b", fileName = Modelica.Utilities.Files.loadResource("modelica://Steps/Resources/Data/kim_2012.txt"), table = fill(0.0, 9, 2), smoothness = Modelica.Blocks.Types.Smoothness.LinearSegments, columns = 2:2);
  
  Modelica.Blocks.Types.ExternalCombiTable1D table_c = Modelica.Blocks.Types.ExternalCombiTable1D(tableName = "5d_c", fileName = Modelica.Utilities.Files.loadResource("modelica://Steps/Resources/Data/kim_2012.txt"), table = fill(0.0, 9, 2), smoothness = Modelica.Blocks.Types.Smoothness.LinearSegments, columns = 2:2);
  
  Modelica.Blocks.Types.ExternalCombiTable1D table_d = Modelica.Blocks.Types.ExternalCombiTable1D(tableName = "5d_d", fileName = Modelica.Utilities.Files.loadResource("modelica://Steps/Resources/Data/kim_2012.txt"), table = fill(0.0, 9, 2), smoothness = Modelica.Blocks.Types.Smoothness.LinearSegments, columns = 2:2);

algorithm

    // determine fitting constant by pitch and hydraulic diameter
    if(abs(pitch - 12.3e-3) <= abs(pitch - 24.6e-3)) then //close to pitch = 12.3
      
      if(abs(d_h - 0.922e-3) <= abs(d_h - 1.222e-3)) then // clost to d_h 0.922
        table_a := table_4b_a;
        table_b := table_4b_b;
        table_c := table_5b_c;
        table_d := table_5b_d;
      else
        //default value - table_4d, 5d should be used here
        table_a := table_4b_a;
        table_b := table_4b_b;
        table_c := table_5b_c;
        table_d := table_5b_d;        
      end if;
      
    else
    
      if(abs(d_h - 0.922e-3) <= abs(d_h - 1.222e-3)) then
        table_a := table_4a_a;
        table_b := table_4a_b;
        table_c := table_5a_c;
        table_d := table_5a_d;
      else
        table_a := table_4c_a;
        table_b := table_4c_b;
        table_c := table_5c_c;
        table_d := table_5c_d;
      end if;       
      
    end if;
    
    /*
    if (MyUtil.sameValue(pitch, 24.6 * 1e-3) and MyUtil.sameValue(d_h, 0.922 * 1e-3)) then
      table_a := table_4a_a;
      table_b := table_4a_b;
      table_c := table_5a_c;
      table_d := table_5a_d;
    elseif (MyUtil.sameValue(pitch, 12.3 * 1e-3) and MyUtil.sameValue(d_h, 0.922 * 1e-3)) then
      table_a := table_4b_a;
      table_b := table_4b_b;
      table_c := table_5b_c;
      table_d := table_5b_d;
    elseif (MyUtil.sameValue(pitch, 24.6 * 1e-3) and MyUtil.sameValue(d_h, 1.222 * 1e-3)) then
      table_a := table_4c_a;
      table_b := table_4c_b;
      table_c := table_5c_c;
      table_d := table_5c_d;
    end if;      
    */
    a := TB.CombiTable1D.getTableValue(table_a, icol = 1, u = Modelica.SIunits.Conversions.to_deg(phi), tableAvailable = 1.0);
    
    b := TB.CombiTable1D.getTableValue(table_b, icol = 1, u = Modelica.SIunits.Conversions.to_deg(phi), tableAvailable = 1.0);
    
    c := TB.CombiTable1D.getTableValue(table_c, icol = 1, u = Modelica.SIunits.Conversions.to_deg(phi), tableAvailable = 1.0);
    
    d := TB.CombiTable1D.getTableValue(table_d, icol = 1, u = Modelica.SIunits.Conversions.to_deg(phi), tableAvailable = 1.0);  

end KimCorrelations;
