within Steps.Components;

model MockPCHeatExchanger
  "Mock Printed Circuit based Heat Exchanger for debug purpose"
  extends BaseExchanger;
  
  import CP = Steps.Utilities.CoolProp; 
  import TB = Modelica.Blocks.Tables;  
  import UTIL = Modelica.Utilities;
  
  replaceable package PBMedia = Steps.Media.SCO2;    
  
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
  
  function thermal_conductivity "cal thermal conductivity of a material"
    extends Modelica.Icons.Function;
    
    input String name "material name";
    input Modelica.SIunits.Temp_C temperature;
    input Modelica.Blocks.Types.ExternalCombiTable1D tableID;
    output Modelica.SIunits.ThermalConductivity k;
    
  algorithm
    if(UTIL.Strings.compare(name, "inconel 750") == Modelica.Utilities.Types.Compare.Equal) then
      // FIX IT: icol = 0， 1？ 
      k := TB.CombiTable1D.getTableValue(tableID, icol = 1, u = temperature, tableAvailable = 1.0); 
    else
      k := 16.2;
    end if;
    
  end thermal_conductivity;    
 
  function E2I
    input PipeType enum;
    output Integer index;
    
    algorithm
      index := Integer(enum);
  end E2I;
  
  /*
  inner Modelica.Blocks.Types.ExternalCombiTable1D table_4a_a = Modelica.Blocks.Types.ExternalCombiTable1D(tableName = "4a_a", fileName = Modelica.Utilities.Files.loadResource("modelica://Steps/Resources/Data/kim_2012.txt"), table = fill(0.0, 9, 2), smoothness = Modelica.Blocks.Types.Smoothness.LinearSegments, columns = 2:2) "Table 4a - column a in Kim[2012] for pitch=24.6, dh=0.922 (dc=1.3 mm))";

  inner Modelica.Blocks.Types.ExternalCombiTable1D table_4a_b = Modelica.Blocks.Types.ExternalCombiTable1D(tableName = "4a_b", fileName = Modelica.Utilities.Files.loadResource("modelica://Steps/Resources/Data/kim_2012.txt"), table = fill(0.0, 9, 2), smoothness = Modelica.Blocks.Types.Smoothness.LinearSegments, columns = 2:2) "Table 4a - column b in Kim[2012] for pitch=12.3, dh=0.922 (dc=1.3 mm))";
  
  inner Modelica.Blocks.Types.ExternalCombiTable1D table_4b_a = Modelica.Blocks.Types.ExternalCombiTable1D(tableName = "4b_a", fileName = Modelica.Utilities.Files.loadResource("modelica://Steps/Resources/Data/kim_2012.txt"), table = fill(0.0, 9, 2), smoothness = Modelica.Blocks.Types.Smoothness.LinearSegments, columns = 2:2) "Table 4b - column a in Kim[2012] for pitch=24.6, dh=1.222 (dc=1.3 mm))";  
  
  inner Modelica.Blocks.Types.ExternalCombiTable1D table_4b_b = Modelica.Blocks.Types.ExternalCombiTable1D(tableName = "4b_b", fileName = Modelica.Utilities.Files.loadResource("modelica://Steps/Resources/Data/kim_2012.txt"), table = fill(0.0, 9, 2), smoothness = Modelica.Blocks.Types.Smoothness.LinearSegments, columns = 2:2) "Table 4b - column b default table";  
  
  inner Modelica.Blocks.Types.ExternalCombiTable1D table_4c_a = Modelica.Blocks.Types.ExternalCombiTable1D(tableName = "4c_a", fileName = Modelica.Utilities.Files.loadResource("modelica://Steps/Resources/Data/kim_2012.txt"), table = fill(0.0, 9, 2), smoothness = Modelica.Blocks.Types.Smoothness.LinearSegments, columns = 2:2) "Table 4c - column a in Kim[2012] for pitch=24.6, dh=0.922 (dc=1.3 mm))";

  inner Modelica.Blocks.Types.ExternalCombiTable1D table_4c_b = Modelica.Blocks.Types.ExternalCombiTable1D(tableName = "4c_b", fileName = Modelica.Utilities.Files.loadResource("modelica://Steps/Resources/Data/kim_2012.txt"), table = fill(0.0, 9, 2), smoothness = Modelica.Blocks.Types.Smoothness.LinearSegments, columns = 2:2) "Table 4c - column b in Kim[2012] for pitch=12.3, dh=0.922 (dc=1.3 mm))";  

  inner Modelica.Blocks.Types.ExternalCombiTable1D table_5a_c = Modelica.Blocks.Types.ExternalCombiTable1D(tableName = "5a_c", fileName = Modelica.Utilities.Files.loadResource("modelica://Steps/Resources/Data/kim_2012.txt"), table = fill(0.0, 9, 2), smoothness = Modelica.Blocks.Types.Smoothness.LinearSegments, columns = 2:2) "Table 5a - column c in Kim[2012] for pitch=24.6, dh=0.922 (dc=1.3 mm))";

  inner Modelica.Blocks.Types.ExternalCombiTable1D table_5a_d = Modelica.Blocks.Types.ExternalCombiTable1D(tableName = "5a_d", fileName = Modelica.Utilities.Files.loadResource("modelica://Steps/Resources/Data/kim_2012.txt"), table = fill(0.0, 9, 2), smoothness = Modelica.Blocks.Types.Smoothness.LinearSegments, columns = 2:2) "Table 5a - column d in Kim[2012] for pitch=12.3, dh=0.922 (dc=1.3 mm))";
  
  inner Modelica.Blocks.Types.ExternalCombiTable1D table_5b_c = Modelica.Blocks.Types.ExternalCombiTable1D(tableName = "5b_c", fileName = Modelica.Utilities.Files.loadResource("modelica://Steps/Resources/Data/kim_2012.txt"), table = fill(0.0, 9, 2), smoothness = Modelica.Blocks.Types.Smoothness.LinearSegments, columns = 2:2) "Table 5b - column c in Kim[2012] for pitch=24.6, dh=1.222 (dc=1.3 mm))";
  
  inner Modelica.Blocks.Types.ExternalCombiTable1D table_5b_d = Modelica.Blocks.Types.ExternalCombiTable1D(tableName = "5b_d", fileName = Modelica.Utilities.Files.loadResource("modelica://Steps/Resources/Data/kim_2012.txt"), table = fill(0.0, 9, 2), smoothness = Modelica.Blocks.Types.Smoothness.LinearSegments, columns = 2:2) "Table 5b - column d default table";
  
  inner Modelica.Blocks.Types.ExternalCombiTable1D table_5c_c = Modelica.Blocks.Types.ExternalCombiTable1D(tableName = "5c_c", fileName = Modelica.Utilities.Files.loadResource("modelica://Steps/Resources/Data/kim_2012.txt"), table = fill(0.0, 9, 2), smoothness = Modelica.Blocks.Types.Smoothness.LinearSegments, columns = 2:2) "Table 5c - column c in Kim[2012] for pitch=24.6, dh=0.922 (dc=1.3 mm))";

  inner Modelica.Blocks.Types.ExternalCombiTable1D table_5c_d = Modelica.Blocks.Types.ExternalCombiTable1D(tableName = "5c_d", fileName = Modelica.Utilities.Files.loadResource("modelica://Steps/Resources/Data/kim_2012.txt"), table = fill(0.0, 9, 2), smoothness = Modelica.Blocks.Types.Smoothness.LinearSegments, columns = 2:2) "Table 5c - column d in Kim[2012] for pitch=12.3, dh=0.922 (dc=1.3 mm))";
    
  inner Modelica.Blocks.Types.ExternalCombiTable1D table_a = Modelica.Blocks.Types.ExternalCombiTable1D(tableName = "4d_a", fileName = Modelica.Utilities.Files.loadResource("modelica://Steps/Resources/Data/kim_2012.txt"), table = fill(0.0, 9, 2), smoothness = Modelica.Blocks.Types.Smoothness.LinearSegments, columns = 2:2);
  
  inner Modelica.Blocks.Types.ExternalCombiTable1D table_b = Modelica.Blocks.Types.ExternalCombiTable1D(tableName = "4d_b", fileName = Modelica.Utilities.Files.loadResource("modelica://Steps/Resources/Data/kim_2012.txt"), table = fill(0.0, 9, 2), smoothness = Modelica.Blocks.Types.Smoothness.LinearSegments, columns = 2:2);
  
  inner Modelica.Blocks.Types.ExternalCombiTable1D table_c = Modelica.Blocks.Types.ExternalCombiTable1D(tableName = "5d_c", fileName = Modelica.Utilities.Files.loadResource("modelica://Steps/Resources/Data/kim_2012.txt"), table = fill(0.0, 9, 2), smoothness = Modelica.Blocks.Types.Smoothness.LinearSegments, columns = 2:2);
  
  inner Modelica.Blocks.Types.ExternalCombiTable1D table_d = Modelica.Blocks.Types.ExternalCombiTable1D(tableName = "5d_d", fileName = Modelica.Utilities.Files.loadResource("modelica://Steps/Resources/Data/kim_2012.txt"), table = fill(0.0, 9, 2), smoothness = Modelica.Blocks.Types.Smoothness.LinearSegments, columns = 2:2);
  */
  inner Modelica.Blocks.Types.ExternalCombiTable1D table_th_inconel_750 = Modelica.Blocks.Types.ExternalCombiTable1D(tableName = "inconel_750", fileName = Modelica.Utilities.Files.loadResource("modelica://Steps/Resources/Data/th_conductivity.txt"), table = fill(0.0, 6, 2), smoothness = Modelica.Blocks.Types.Smoothness.LinearSegments, columns = 2:2) "thermal conductivity for inconel_750";
  
  parameter Integer N_seg = 10 "Number of segments in a tube";
  
  parameter Modelica.SIunits.Length length_cell = 1e-3 "length of a cell";
  
  parameter Modelica.SIunits.Angle phi = 0.0 "unit rad";
  
  parameter Modelica.SIunits.ReynoldsNumber Re_design = 0.0 "On-design ReynoldsNumber";
  
  parameter Modelica.SIunits.Diameter d_c = 0.0 "Diameter of semi-circular channel";
  // inlet/outlet temperature of hot/cool channel
  parameter Modelica.SIunits.Temp_C T_hot_in = 500;
  //parameter Modelica.SIunits.Temp_C T_hot_out;
  parameter Modelica.SIunits.Temp_C T_cool_in = 300;
  //parameter Modelica.SIunits.Temp_C T_cool_out;
  
  //pressure of hot/cool channel
  parameter Modelica.SIunits.Pressure p_hot = 20 * 1e6;
  parameter Modelica.SIunits.Pressure p_cool = 9 * 1e6;
  
  // mass flow of hot/cool channel
  parameter Modelica.SIunits.MassFlowRate m_dot_hot = 5;
  parameter Modelica.SIunits.MassFlowRate m_dot_cool = 5;
  
  parameter Modelica.SIunits.Length pitch = 10 "pitch length of channel";
  
  parameter String name_material = "inconel 750";
  
  Modelica.SIunits.Diameter d_h "Hydraulic Diameter";
  
  Modelica.SIunits.Length peri_c "perimeter of semi-circular";
  
  Modelica.SIunits.Length t_wall "thickness of wall between two neighboring hot and cold";
  
  Integer N_channel "number of channels";
  
  Modelica.SIunits.Area A_c "Area of semi-circular tube";    
  
  Modelica.SIunits.Area A_flow "Flow area of all channels";
  
  Modelica.SIunits.Area A_fc "Area of cold stream area";
  Modelica.SIunits.Area A_fh "Area of hot stream area";
  Modelica.SIunits.Area A_fmax "Area of maximum stream area comparing A_fc and A_fh: A_fmax = max(A_fc, A_fh)";
  
  Modelica.SIunits.Area A_stack "surface area of all cells in a stack";
  
  Modelica.SIunits.DynamicViscosity mu_c "average dynamic Viscosity in cold channel";
  Modelica.SIunits.DynamicViscosity mu_h "average dynamic Viscosity in hot channel"; 
  
  HXCellState [N_seg] state_cell "state array of each cells";
  
  // length of one pipe in HeatExchanger
  Modelica.SIunits.Length length "unit m";

  
  Real fit_const_a "fitting constant a in Eq[3] of [kim, 2011] ";
  
  Real fit_const_b "fitting constant b in Eq[3] of [kim, 2011] ";
  
  Real fit_const_c "fitting constant c in Eq[4] of [kim, 2011] ";
  
  Real fit_const_d "fitting constant d in Eq[4] of [kim, 2011] ";
  
  Modelica.SIunits.MassFlowRate G_c "Mass flux of each cold channel";  
  Modelica.SIunits.MassFlowRate G_h "Mass flux of each cold channel";    
  
  Integer table_select;   
    
  /*
  constant Integer num_pipes = 2;
  
  PipeType pipes[num_pipes] = {PipeType.Hot, PipeType.Cold};  
  
  PBMedia.CO2_pT mediums[num_pipes];
  
  Modelica.SIunits.MassFlowRate mdot_pipes[num_pipes];

  Modelica.SIunits.MassFlowRate G[num_pipes] "Mass flux of each hot channel";
  
 
  Integer idx;
  */
algorithm

  when initial() then
    
    length := length_cell * N_seg "length of a cell";
    peri_c := d_c * Modelica.Constants.pi /2 + d_c;
    
    A_c := Modelica.Constants.pi * d_c * d_c / 8;
    d_h := 4 * A_c / peri_c;    
    
    mu_c := CP.PropsSI("V", "P", p_cool, "T", T_cool_in, PBMedia.mediumName) "average dynamic Viscosity in cold channel";
    mu_h := CP.PropsSI("V", "P", p_hot, "T", T_hot_in, PBMedia.mediumName) "average dynamic Viscosity in hot channel"; 
    
    A_fc := m_dot_cool * d_h / mu_c / Re_design "Area of cold stream area";
    A_fh := m_dot_hot * d_h / mu_h /Re_design "Area of hot stream area";
    A_fmax := max(A_fc, A_fh) "Area of maximum stream area comparing A_fc and A_fh: A_fmax = max(A_fc, A_fh)";    
    
    t_wall := (2 - Modelica.Constants.pi  / 4) * (d_c / 2);
    
    N_channel := integer(A_fmax / A_c);
    
    A_flow := N_channel * A_c;
    
    A_stack := peri_c * length_cell * N_channel "surface area of all cells in a stack";    
    
    table_select := 4;
    // determine fitting constant by pitch and hydraulic diameter  
    /*
    if (sameValue(pitch, 24.6 * 1e-3) and sameValue(d_h, 0.922 * 1e-3)) then
      table_select := 1;
      table_a := table_4a_a;
      table_b := table_4a_b;
      table_c := table_5a_c;
      table_d := table_5a_d;
    elseif (sameValue(pitch, 12.3 * 1e-3) and sameValue(d_h, 0.922 * 1e-3)) then
      table_select := 2;
      table_a := table_4b_a;
      table_b := table_4b_b;
      table_c := table_5b_c;
      table_d := table_5b_d;
    elseif (sameValue(pitch, 24.6 * 1e-3) and sameValue(d_h, 1.222 * 1e-3)) then
      table_select := 3;
      table_a := table_4c_a;
      table_b := table_4c_b;
      table_c := table_5c_c;
      table_d := table_5c_d;
    end if;
    
    fit_const_a := TB.CombiTable1D.getTableValue(table_a, icol = 1, u = Modelica.SIunits.Conversions.to_deg(phi), tableAvailable = 1.0);
    fit_const_b := TB.CombiTable1D.getTableValue(table_b, icol = 1, u = Modelica.SIunits.Conversions.to_deg(phi), tableAvailable = 1.0);
    fit_const_c := TB.CombiTable1D.getTableValue(table_c, icol = 1, u = Modelica.SIunits.Conversions.to_deg(phi), tableAvailable = 1.0);
    fit_const_d := TB.CombiTable1D.getTableValue(table_d, icol = 1, u = Modelica.SIunits.Conversions.to_deg(phi), tableAvailable = 1.0); 
    */
    
    fit_const_a := 0.51539; //TB.CombiTable1D.getTableValue(table_a, icol = 1, u = Modelica.SIunits.Conversions.to_deg(phi), tableAvailable = 1.0);
    fit_const_b := 0.77944; //TB.CombiTable1D.getTableValue(table_b, icol = 1, u = Modelica.SIunits.Conversions.to_deg(phi), tableAvailable = 1.0);
    fit_const_c := 0.04885; //TB.CombiTable1D.getTableValue(table_c, icol = 1, u = Modelica.SIunits.Conversions.to_deg(phi), tableAvailable = 1.0);
    fit_const_d := 0.67896; //TB.CombiTable1D.getTableValue(table_d, icol = 1, u = Modelica.SIunits.Conversions.to_deg(phi), tableAvailable = 1.0); 

  end when;
    
equation

  medium_hot_in.state = PBMedia.setState_pTX(p = inlet_hot.p, T = inlet_hot.T);
  medium_cool_out.state = PBMedia.setState_pTX(p = outlet_cool.p, T = outlet_cool.T);  
  
algorithm  
  // initialize the state for the hot side and cold side in all the cells
  // using the hot inlet and cool outlet.    

  state_cell[1].T_c := medium_cool_out.T;
  state_cell[1].T_h := medium_hot_in.T;
  
  state_cell[1].p_c := medium_cool_out.p;
  state_cell[1].p_h := medium_hot_in.p;  
  
  state_cell[1].mdot_c := - outlet_cool.m_flow;
  state_cell[1].mdot_h := inlet_hot.m_flow;  

  state_cell[1].medium_name_c := medium_cool_out.mediumName;
  state_cell[1].medium_name_h := medium_hot_in.mediumName;  
  
  state_cell[1].h_mass_c := medium_cool_out.h;
  state_cell[1].h_mass_h := medium_hot_in.h;  
  
  G_c := state_cell[1].mdot_c / N_channel / A_c;
  G_h := state_cell[1].mdot_h / N_channel / A_c;
  
  for i in 1 : N_seg loop // index out of range??
    state_cell[i].length := length_cell;
   
    //Debug from this point
    state_cell[i].mu_c := CP.PropsSI("V", "P", state_cell[i].p_c, "T", state_cell[i].T_c, state_cell[i].medium_name_c);
    state_cell[i].mu_h := CP.PropsSI("V", "P", state_cell[i].p_h, "T", state_cell[i].T_h, state_cell[i].medium_name_h);
      
    state_cell[i].k_c := CP.PropsSI("L", "P", state_cell[i].p_c, "T", state_cell[i].T_c, state_cell[i].medium_name_c);  
    
    assert(state_cell[i].k_c > 0 and state_cell[i].k_c < 1e5, 
    "outside range " + keyvalStr("k_c", state_cell[i].k_c) + 
    " at " + debugInfo({"i", "T", "P"}, {i, state_cell[i].T_c, state_cell[i].p_c}));
      
    state_cell[i].k_h := CP.PropsSI("L", "P", state_cell[i].p_h, "T", state_cell[i].T_h, state_cell[i].medium_name_h);  
      
    state_cell[i].Re_c := G_c * d_h / state_cell[i].mu_c;
    state_cell[i].Re_h := G_h * d_h / state_cell[i].mu_h;
      
    state_cell[i].rho_c := CP.PropsSI("D", "P", state_cell[i].p_c, "T", state_cell[i].T_c, state_cell[i].medium_name_c);
    state_cell[i].rho_h := CP.PropsSI("D", "P", state_cell[i].p_h, "T", state_cell[i].T_h, state_cell[i].medium_name_h);
      
    state_cell[i].u_c := state_cell[i].mdot_c / A_flow / state_cell[i].rho_c;  
    state_cell[i].u_h := state_cell[i].mdot_h / A_flow / state_cell[i].rho_h;     
    
    assert(state_cell[i].Re_c > 0 and state_cell[i].Re_c < 1e6,
    "out of range " + keyvalStr("Re_c", state_cell[i].Re_c) + " at " + 
    debugInfo({"i","G","d_h","mu"},{i, G_c, d_h, state_cell[i].mu_c}));  
    
    state_cell[i].Nu_c := 4.089 + fit_const_c * (state_cell[i].Re_c ^ fit_const_d);
    
    assert(state_cell[i].Re_h > 0 and state_cell[i].Re_h < 1e6,
    "out of range " + keyvalStr("Re_h", state_cell[i].Re_h) + " at " + 
    debugInfo({"i","G","d_h","mu"},{i, G_h, d_h, state_cell[i].mu_h}));       
    
    state_cell[i].Nu_h := 4.089 + fit_const_c * (state_cell[i].Re_h ^ fit_const_d);
      
    state_cell[i].h_c := state_cell[i].Nu_c * state_cell[i].k_c / d_h;
    state_cell[i].h_h := state_cell[i].Nu_h * state_cell[i].k_h / d_h;
      
    state_cell[i].f_c := (15.78 + fit_const_a * state_cell[i].Re_c ^ fit_const_b ) / state_cell[i].Re_c;      
    state_cell[i].f_h := (15.78 + fit_const_a * state_cell[i].Re_h ^ fit_const_b ) / state_cell[i].Re_h;      
    
    state_cell[i].dp_c := state_cell[i].f_c * state_cell[i].length * state_cell[i].rho_c *  (state_cell[i].u_c ^ 2) / d_h;    
    state_cell[i].dp_h := state_cell[i].f_h * state_cell[i].length * state_cell[i].rho_h *  (state_cell[i].u_h ^ 2) / d_h;
    
    // no use of following parameters, use default value
    state_cell[i].Pr_c := 1.0;
    state_cell[i].Pr_h := 1.0;
    
    state_cell[i].k_wall := thermal_conductivity(tableID = table_th_inconel_750, name = name_material, temperature = (state_cell[i].T_c + state_cell[i].T_h) / 2);
    
    state_cell[i].U := 1 / ( 1 / state_cell[i].h_h + 1 / state_cell[i].h_c + t_wall / state_cell[i].k_wall);    

    
    if state_cell[i].T_h > state_cell[i].T_c then
      state_cell[i].q := state_cell[i].U * A_stack * (state_cell[i].T_h - state_cell[i].T_c);      
    else
      state_cell[i].q := 0;
    end if;  
    
    // Calculate state for next cell
    if i == N_seg then
      break;
    end if;
    
    // set parameter for next cell of hot and cold pipe
    state_cell[i + 1].medium_name_c := state_cell[i].medium_name_c;      
    state_cell[i + 1].medium_name_h := state_cell[i].medium_name_h;      
      
    state_cell[i + 1].mdot_c := state_cell[i].mdot_c;
    state_cell[i + 1].mdot_h := state_cell[i].mdot_h;
      
    state_cell[i + 1].h_mass_c := (state_cell[i].h_mass_c * state_cell[i].mdot_c - state_cell[i].q) / state_cell[i].mdot_c;
    /*
    assert(i <> 1,
    "out of range " + keyvalStr("h_mass_c", state_cell[i+1].h_mass_c) + " at " + 
    debugInfo({"i+1","h_mass_c[i]","m_dot_inlet_cold","q[i]"},{i+1,state_cell[i].h_mass_c, state_cell[i].mdot_c, state_cell[i].q})); 
    */
    assert(state_cell[i+1].h_mass_c > 0 or state_cell[i+1].h_mass_c < 1e6,
    "out of range " + keyvalStr("h_mass_c", state_cell[i+1].h_mass_c) + " at " + 
    debugInfo({"i+1","h_mass_c[i]","m_dot_inlet_cold","q[i]"},{i+1,state_cell[i].h_mass_c, state_cell[i].mdot_c, state_cell[i].q})); 
    
    state_cell[i + 1].h_mass_h := (state_cell[i].h_mass_h * state_cell[i].mdot_h - state_cell[i].q) / state_cell[i].mdot_h; 
     
    state_cell[i + 1].p_c := state_cell[i].p_c - state_cell[i].dp_c;
    state_cell[i + 1].p_h := state_cell[i].p_h - state_cell[i].dp_h;
      
    state_cell[i + 1].T_c := CP.PropsSI("T", "P",state_cell[i + 1].p_c, "H", state_cell[i + 1].h_mass_c, state_cell[i].medium_name_c);
    
    assert(state_cell[i+1].T_c > 0 and state_cell[i+1].T_c < 1e6, 
    "out of range " + keyvalStr("T_c", state_cell[i+1].T_c) + " at " + debugInfo({"i","P","h_mass_c"},{i,state_cell[i + 1].p_c,state_cell[i + 1].h_mass_c}));     
    
    state_cell[i + 1].T_h := CP.PropsSI("T", "P",state_cell[i + 1].p_h, "H", state_cell[i + 1].h_mass_h, state_cell[i].medium_name_h);

  end for;
  
equation 
  
  medium_cool_in.state = PBMedia.setState_phX(p = state_cell[N_seg].p_c, h = state_cell[N_seg].h_mass_c);   
  inlet_cool.T = medium_cool_in.T;
  inlet_cool.p = medium_cool_in.p;
  inlet_cool.h_outflow = state_cell[N_seg].h_mass_c;  
  //inlet_cool.h_outflow = inStream(outlet_cool.h_outflow);  
  //outlet_cool.h_outflow = - medium_cool_out.h;  
  inlet_cool.h_outflow = - outlet_cool.h_outflow;
  outlet_cool.m_flow + inlet_cool.m_flow = 0;
  
  // for outlet_hot
  medium_hot_out.state = PBMedia.setState_phX(p = state_cell[N_seg].p_h, h = state_cell[N_seg].h_mass_h); 
  outlet_hot.T = medium_hot_out.T;
  outlet_hot.p = medium_hot_out.p;
  outlet_hot.h_outflow = - state_cell[N_seg].h_mass_h;    
  inlet_hot.h_outflow = inStream(outlet_hot.h_outflow);  
  outlet_hot.m_flow + inlet_hot.m_flow = 0;

end MockPCHeatExchanger;
