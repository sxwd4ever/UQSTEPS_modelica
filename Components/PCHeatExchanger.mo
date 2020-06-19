within Steps.Components;

model PCHeatExchanger
  "Printed Circuit based Heat Exchanger"
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
  
  function thermal_conductivity "cal thermal conductivity of a material"
    extends Modelica.Icons.Function;
    
    input String name "material name";
    input Modelica.SIunits.Temp_C temperature;
    output Modelica.SIunits.ThermalConductivity k;

protected
    
    inner Modelica.Blocks.Types.ExternalCombiTable1D table_th_inconel_750 = Modelica.Blocks.Types.ExternalCombiTable1D(tableName = "inconel_750", fileName = Modelica.Utilities.Files.loadResource("modelica://Steps/Resources/Data/th_conductivity.txt"), table = fill(0.0, 6, 2), smoothness = Modelica.Blocks.Types.Smoothness.LinearSegments, columns = 2:2) "thermal conductivity for inconel_750";
    
  algorithm
    if(UTIL.Strings.compare(name, "inconel_75") == Modelica.Utilities.Types.Compare.Equal) then
      // FIX IT: icol = 0， 1？ 
      k := TB.CombiTable1D.getTableValue(table_th_inconel_750, icol = 1, u = temperature, tableAvailable = 0.0); 
    else
      k := 16.2;
    end if;
    
  end thermal_conductivity;
  
  inner Modelica.Blocks.Types.ExternalCombiTable2D table_4a = Modelica.Blocks.Types.ExternalCombiTable2D(tableName = "4a", fileName = Modelica.Utilities.Files.loadResource("modelica://Steps/Resources/Data/kim_2012.txt"), table = fill(0.0, 9, 3), smoothness = Modelica.Blocks.Types.Smoothness.LinearSegments) "Table 4a in Kim[2012] for pitch=24.6, dh=0.922 (dc=1.3 mm))";

  inner Modelica.Blocks.Types.ExternalCombiTable2D table_4b = Modelica.Blocks.Types.ExternalCombiTable2D(tableName = "4b", fileName = Modelica.Utilities.Files.loadResource("modelica://Steps/Resources/Data/kim_2012.txt"), table = fill(0.0, 9, 3), smoothness = Modelica.Blocks.Types.Smoothness.LinearSegments) "Table 4b in Kim[2012] for pitch=12.3, dh=0.922 (dc=1.3 mm))";
  
  inner Modelica.Blocks.Types.ExternalCombiTable2D table_4c = Modelica.Blocks.Types.ExternalCombiTable2D(tableName = "4c", fileName = Modelica.Utilities.Files.loadResource("modelica://Steps/Resources/Data/kim_2012.txt"), table = fill(0.0, 9, 3), smoothness = Modelica.Blocks.Types.Smoothness.LinearSegments) "Table 4c in Kim[2012] for pitch=24.6, dh=1.222 (dc=1.3 mm))";  
  
  inner Modelica.Blocks.Types.ExternalCombiTable2D table_4d = Modelica.Blocks.Types.ExternalCombiTable2D(tableName = "4d", fileName = Modelica.Utilities.Files.loadResource("modelica://Steps/Resources/Data/kim_2012.txt"), table = fill(0.0, 9, 3), smoothness = Modelica.Blocks.Types.Smoothness.LinearSegments) "Table 5d default table";  

  inner Modelica.Blocks.Types.ExternalCombiTable2D table_5a = Modelica.Blocks.Types.ExternalCombiTable2D(tableName = "5a", fileName = Modelica.Utilities.Files.loadResource("modelica://Steps/Resources/Data/kim_2012.txt"), table = fill(0.0, 9, 3), smoothness = Modelica.Blocks.Types.Smoothness.LinearSegments) "Table 5a in Kim[2012] for pitch=24.6, dh=0.922 (dc=1.3 mm))";

  inner Modelica.Blocks.Types.ExternalCombiTable2D table_5b = Modelica.Blocks.Types.ExternalCombiTable2D(tableName = "5b", fileName = Modelica.Utilities.Files.loadResource("modelica://Steps/Resources/Data/kim_2012.txt"), table = fill(0.0, 9, 3), smoothness = Modelica.Blocks.Types.Smoothness.LinearSegments) "Table 5b in Kim[2012] for pitch=12.3, dh=0.922 (dc=1.3 mm))";
  
  inner Modelica.Blocks.Types.ExternalCombiTable2D table_5c = Modelica.Blocks.Types.ExternalCombiTable2D(tableName = "5c", fileName = Modelica.Utilities.Files.loadResource("modelica://Steps/Resources/Data/kim_2012.txt"), table = fill(0.0, 9, 3), smoothness = Modelica.Blocks.Types.Smoothness.LinearSegments) "Table 5c in Kim[2012] for pitch=24.6, dh=1.222 (dc=1.3 mm))";
  
  inner Modelica.Blocks.Types.ExternalCombiTable2D table_5d = Modelica.Blocks.Types.ExternalCombiTable2D(tableName = "5d", fileName = Modelica.Utilities.Files.loadResource("modelica://Steps/Resources/Data/kim_2012.txt"), table = fill(0.0, 9, 3), smoothness = Modelica.Blocks.Types.Smoothness.LinearSegments) "Table 5d default table";
  
  inner Modelica.Blocks.Types.ExternalCombiTable2D table_4;
  inner Modelica.Blocks.Types.ExternalCombiTable2D table_5;  

  // length of one pipe in HeatExchanger
  parameter Modelica.SIunits.Length length = 1.0 "unit m";
  
  parameter Integer N_seg = 1 "Number of segments in a tube";
  
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
  
  Modelica.SIunits.Length length_cell "length of a cell";
  
  Real fit_const_a "fitting constant a in Eq[3] of [kim, 2011] ";
  
  Real fit_const_b "fitting constant b in Eq[3] of [kim, 2011] ";
  
  Real fit_const_c "fitting constant c in Eq[4] of [kim, 2011] ";
  
  Real fit_const_d "fitting constant d in Eq[4] of [kim, 2011] ";
  
  Modelica.SIunits.MassFlowRate G_h "Mass flux of each hot channel";
  Modelica.SIunits.MassFlowRate G_c "Mass flux of each cold channel";   
  
algorithm

  when initial() then
    
    peri_c := d_c * Modelica.Constants.pi /2 + d_c;
    
    A_c := Modelica.Constants.pi * d_c * d_c / 8;
    
    mu_c := CP.PropsSI("V", "P", p_cool, "T", T_cool_in, PBMedia.mediumName) "average dynamic Viscosity in cold channel";
    mu_h := CP.PropsSI("V", "P", p_hot, "T", T_hot_in, PBMedia.mediumName) "average dynamic Viscosity in hot channel"; 
    
    A_fc := m_dot_cool * d_h / mu_c / Re_design "Area of cold stream area";
    A_fh := m_dot_hot * d_h / mu_h /Re_design "Area of hot stream area";
    A_fmax := max(A_fc, A_fh) "Area of maximum stream area comparing A_fc and A_fh: A_fmax = max(A_fc, A_fh)";
  
    d_h := 4 * A_c / peri_c;
    
    t_wall := (2 - Modelica.Constants.pi  / 4) * (d_c / 2);
    
    N_channel := integer(A_fmax / A_c);
    
    A_flow := N_channel * A_c;
    
    A_stack := peri_c * length_cell * N_channel "surface area of all cells in a stack";
    
    length_cell := length / N_seg "length of a cell";
    
    // determine fitting constant by pitch and hydraulic diameter  
    if (sameValue(pitch, 24.6 * 1e-3) and sameValue(d_h, 0.922 * 1e-3)) then
      table_4 := table_4a;
      table_5 := table_5a;
    elseif (sameValue(pitch, 12.3 * 1e-3) and sameValue(d_h, 0.922 * 1e-3)) then
      table_4 := table_4b;
      table_5 := table_5b;
    elseif (sameValue(pitch, 24.6 * 1e-3) and sameValue(d_h, 1.222 * 1e-3)) then
      table_4 := table_4c;
      table_5 := table_5c;
    else  
      table_4 := table_4d;
      table_5 := table_5d;
    end if;
    
    fit_const_a := TB.CombiTable2D.getTableValue(table_4, u1 = phi, u2 = 1, tableAvailable = 0.0);
    fit_const_b := TB.CombiTable2D.getTableValue(table_4, u1 = phi, u2 = 2, tableAvailable = 0.0);
    fit_const_c := TB.CombiTable2D.getTableValue(table_5, u1 = phi, u2 = 1, tableAvailable = 0.0);
    fit_const_d := TB.CombiTable2D.getTableValue(table_5, u1 = phi, u2 = 2, tableAvailable = 0.0);  
  
  end when;
    
equation

  medium_hot_in.state = PBMedia.setState_pTX(p = inlet_hot.p, T = inlet_hot.T);
  medium_cool_in.state = PBMedia.setState_pTX(p = inlet_cool.p, T = inlet_cool.T);   
  
algorithm  

  G_c := inlet_cool.m_flow / N_channel / A_c;
  G_h := inlet_hot.m_flow / N_channel / A_c;
  
  // Temperature
  state_cell[1].T_h := medium_hot_in.T;
  state_cell[1].T_c := medium_cool_in.T;
  
  // Pressure
  state_cell[1].p_h := medium_hot_in.p;
  state_cell[1].p_c := medium_cool_in.p;
  
  // specific enthalpy
  state_cell[1].h_mass_h := CP.PropsSI("H", "P", medium_hot_in.p, "T", medium_hot_in.T, medium_hot_in.mediumName);
  state_cell[1].h_mass_c := CP.PropsSI("H", "P", medium_cool_in.p, "T", medium_cool_in.T, medium_cool_in.mediumName);  
  
  for i in 1 : N_seg loop // index out of range??
    state_cell[i].length := length_cell;
      
    state_cell[i].mu_h := CP.PropsSI("V", "P", state_cell[i].p_h, "T", state_cell[i].T_h, medium_hot_in.mediumName);
    state_cell[i].mu_c := CP.PropsSI("V", "P", state_cell[i].p_c, "T", state_cell[i].T_c, medium_cool_in.mediumName);
    
    state_cell[i].k_h := CP.PropsSI("L", "P", state_cell[i].p_h, "T", state_cell[i].T_h, medium_hot_in.mediumName);
    state_cell[i].k_c := CP.PropsSI("L", "P", state_cell[i].p_c, "T", state_cell[i].T_c, medium_cool_in.mediumName);  
    
    state_cell[i].Re_h := G_h * d_h / state_cell[i].mu_h;
    state_cell[i].Re_c := G_c * d_h / state_cell[i].mu_c; 
    
    state_cell[i].rho_h := CP.PropsSI("D", "P", state_cell[i].p_h, "T", state_cell[i].T_h, medium_hot_in.mediumName);
    state_cell[i].rho_c := CP.PropsSI("D", "P", state_cell[i].p_c, "T", state_cell[i].T_c, medium_cool_in.mediumName);   
    
    state_cell[i].u_h := inlet_hot.m_flow / A_stack / state_cell[i].rho_h;   
    state_cell[i].u_c := inlet_cool.m_flow / A_stack / state_cell[i].rho_c; 
    
    state_cell[i].Nu_h := 4.089 + fit_const_c * (state_cell[i].Re_h^ fit_const_d);
    state_cell[i].Nu_c := 4.089 + fit_const_c * (state_cell[i].Re_c^ fit_const_d);
    
    state_cell[i].k_wall := thermal_conductivity(name_material, (state_cell[i].T_h + state_cell[i].T_c) / 2);
    
    state_cell[i].h_h := state_cell[i].Nu_h * state_cell[i].k_h / d_h;
    state_cell[i].h_c := state_cell[i].Nu_c * state_cell[i].k_c / d_h;
    
    state_cell[i].U := 1 / ( 1 / state_cell[i].h_h + 1 / state_cell[i].h_c + t_wall / state_cell[i].k_wall);
    
    state_cell[i].f_h := (15.78 + fit_const_a * state_cell[i].Re_h ^ fit_const_b ) / state_cell[i].Re_h;
    state_cell[i].f_c := (15.78 + fit_const_a * state_cell[i].Re_c ^ fit_const_b ) / state_cell[i].Re_c;
    
    if state_cell[i].T_h > state_cell[i].T_c then
      state_cell[i].q := state_cell[i].U * A_stack * (state_cell[i].T_h - state_cell[i].T_c);      
    else
      state_cell[i].q := 0;
    end if;
    
    state_cell[i].dp_h := state_cell[i].f_h * state_cell[i].length * state_cell[i].rho_h *  (state_cell[i].u_h ^ 2) / d_h;
    state_cell[i].dp_c := state_cell[i].f_c * state_cell[i].length * state_cell[i].rho_c *  (state_cell[i].u_c ^ 2) / d_c;
    
    // no use of following parameters, use default value
    state_cell[i].Pr_c := 0;
    state_cell[i].Pr_h := 0;    
    
    // Calculate state for next cell
    if i == N_seg then
      break;
    end if;
    
    state_cell[i + 1].h_mass_h := (state_cell[i].h_mass_h * inlet_hot.m_flow - state_cell[i].q) / inlet_hot.m_flow;
    state_cell[i + 1].h_mass_c := (state_cell[i].h_mass_c * inlet_cool.m_flow - state_cell[i].q) / inlet_cool.m_flow;
    
    state_cell[i + 1].p_h := state_cell[i].p_h - state_cell[i].dp_h;
    state_cell[i + 1].p_c := state_cell[i].p_c - state_cell[i].dp_c;
    
    state_cell[i + 1].T_h := CP.PropsSI("T", "P",state_cell[i + 1].p_h, "H", state_cell[i + 1].h_mass_h, medium_hot_in.mediumName);
    state_cell[i + 1].T_c := CP.PropsSI("T", "P",state_cell[i + 1].p_c, "H", state_cell[i + 1].h_mass_c, medium_cool_in.mediumName);
      
  end for;
  
equation 
  
  // for out_cool
  medium_cool_out.state = PBMedia.setState_phX(p = state_cell[N_seg].p_c,  h = state_cell[N_seg].h_mass_h);
  outlet_cool.T = medium_cool_out.T;
  outlet_cool.p = medium_cool_out.p;
  outlet_cool.h_outflow = medium_cool_out.h;  
  
  inlet_cool.h_outflow = inStream(outlet_cool.h_outflow);    
  outlet_cool.m_flow + inlet_cool.m_flow = 0;
  
  // for outlet_hot
  medium_hot_out.state = PBMedia.setState_phX(p = state_cell[N_seg].p_h, h = state_cell[N_seg].h_mass_c); 
  outlet_hot.T = medium_hot_out.T;
  outlet_hot.p = medium_hot_out.p;
  outlet_hot.h_outflow = medium_hot_out.h;  
  
  inlet_hot.h_outflow = inStream(outlet_hot.h_outflow);  
  outlet_hot.m_flow + inlet_hot.m_flow = 0;

end PCHeatExchanger;
