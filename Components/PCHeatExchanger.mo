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
  
  /*
  function calCellState "calculate state of a cell-pair in an exchanger"
    extends Modelica.Icons.Function;    
    input HXCellState cellState ;
    input String fluid_name_hot;
    input String fluid_name_cold;
    
    output Modelica.SIunits.Temp_C T_cool_next;
    output Modelica.SIunits.Temp_C T_hot_next;
    output Modelica.SIunits.Pressure p_cool_next;
    output Modelica.SIunits.Pressure p_hot_next;
  algorithm    
  
    cellState.mu_h := CP.PropsSI("V", "P", cellState.p_h, "T", cellState.T_h, fluid_name_hot);
    cellState.mu_c := CP.PropsSI("V", "P", cellState.p_h, "T", cellState.T_h, fluid_name_cold);
    
  end calCellState;*/
  
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
  
  parameter Integer N_seg = 10 "Number of segments in a tube";
  
  parameter Modelica.SIunits.Angle phi = 0.0 "unit rad";
  
  parameter Modelica.SIunits.ReynoldsNumber Re_design = 0.0 "On-design ReynoldsNumber";
  
  parameter Modelica.SIunits.Diameter d_c = 0.0 "Diameter of semi-circular channel";
  // inlet/outlet temperature of hot/cool channel
  parameter Modelica.SIunits.Temp_C T_hot_in;
  //parameter Modelica.SIunits.Temp_C T_hot_out;
  parameter Modelica.SIunits.Temp_C T_cool_in;
  //parameter Modelica.SIunits.Temp_C T_cool_out;
  
  //pressure of hot/cool channel
  parameter Modelica.SIunits.Pressure p_hot;
  parameter Modelica.SIunits.Pressure p_cool;
  
  // mass flow of hot/cool channel
  parameter Modelica.SIunits.MassFlowRate m_dot_hot;
  parameter Modelica.SIunits.MassFlowRate m_dot_cool;
  
  parameter Modelica.SIunits.Length pitch "pitch length of channel";
  
  parameter String name_material = "inconel 750";
  
  Modelica.SIunits.Diameter d_h = 0.0 "Hydraulic Diameter";
  
  Modelica.SIunits.Length peri_c = 0.0 "perimeter of semi-circular";
  
  Modelica.SIunits.Length t_wall "thickness of wall between two neighboring hot and cold";
  
  Integer N_channel = 10 "number of channels";
  
  Modelica.SIunits.Area A_fc = 1.0 "Area of cold stream area";
  Modelica.SIunits.Area A_fh = 1.0 "Area of hot stream area";
  Modelica.SIunits.Area A_fmax = 1.0 "Area of maximum stream area comparing A_fc and A_fh: A_fmax = max(A_fc, A_fh)";
  
  Modelica.SIunits.Area A_c = 1.0 "Area of semi-circular tube";    
  
  Modelica.SIunits.Area A_flow = 1.0 "Flow area of all channels";
  
  Modelica.SIunits.Area A_stack = 1.0 "surface area of all cells in a stack";
  
  Modelica.SIunits.DynamicViscosity mu_c "average dynamic Viscosity in cold channel";
  Modelica.SIunits.DynamicViscosity mu_h "average dynamic Viscosity in hot channel"; 
  
  HXCellState [N_seg] state_cell "state array of each cells";
  
  Modelica.SIunits.Length length_cell = 1.0 "length of a cell";
  
  Real fit_const_a = 0.0 "fitting constant a in Eq[3] of [kim, 2011] ";
  
  Real fit_const_b = 0.0 "fitting constant b in Eq[3] of [kim, 2011] ";
  
  Real fit_const_c = 0.0 "fitting constant c in Eq[4] of [kim, 2011] ";
  
  Real fit_const_d = 0.0 "fitting constant d in Eq[4] of [kim, 2011] ";
  
  Modelica.SIunits.MassFlowRate G_h "Mass flux of each hot channel";
  Modelica.SIunits.MassFlowRate G_c "Mass flux of each cold channel"; 

  HXCellState state_cur "state of current cell";
  HXCellState state_next "state of next cell";
  
initial algorithm

  A_c := Modelica.Constants.pi * d_c * d_c / 8; 
  peri_c := d_c * Modelica.Constants.pi /2 + d_c;
  d_h := 4 * A_c / peri_c;
  t_wall := (2 - Modelica.Constants.pi  / 4) * (d_c / 2);
  
  mu_c := CP.PropsSI("V", "P", p_cool, "T", T_cool_in, medium_cool_in.mediumName);
  mu_h := CP.PropsSI("V", "P", p_hot, "T", T_hot_in, medium_hot_in.mediumName);
  
  A_fc := m_dot_cool * d_h / mu_c / Re_design;
  A_fh := m_dot_hot * d_h / mu_h /Re_design;
  
  A_fmax := max(A_fc, A_fh);
  N_channel := integer(A_fmax / A_c);
  A_flow := N_channel * A_c;
  length_cell := length / N_seg;
  A_stack := peri_c * length_cell * N_channel;
  
  for i in 0 : N_seg loop
    state_cell[i].length := length / N_seg;
  end for;
  
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
  
equation
  
  medium_hot_in.state = PBMedia.setState_pTX(p = inlet_hot.p, T = inlet_hot.T);
  medium_cool_in.state = PBMedia.setState_pTX(p = inlet_cool.p, T = inlet_cool.T);   

algorithm   
  G_c := inlet_cool.m_flow / N_channel / A_c;
  G_h := inlet_hot.m_flow / N_channel / A_c;
  
  // Temperature
  state_cell[0].T_h := inlet_hot.T;
  state_cell[0].T_c := inlet_cool.T;
  
  // Pressure
  state_cell[0].p_h := inlet_hot.p;
  state_cell[0].p_c := inlet_cool.p;
  
  // specific enthalpy
  state_cell[0].h_mass_h := CP.PropsSI("H", "P", inlet_hot.p, "T", inlet_hot.T, medium_hot_in.mediumName);
  state_cell[0].h_mass_c := CP.PropsSI("H", "P", inlet_cool.p, "T", inlet_cool.T, medium_cool_in.mediumName);  
  
  for i in 0 : N_seg loop // index out of range??
    state_cur := state_cell[i];
    
    state_cur.length := length_cell;
    
    state_cur.mu_h := CP.PropsSI("V", "P", state_cur.p_h, "T", state_cur.T_h, medium_hot_in.mediumName);
    state_cur.mu_c := CP.PropsSI("V", "P", state_cur.p_c, "T", state_cur.T_c, medium_cool_in.mediumName);
    
    state_cur.k_h := CP.PropsSI("L", "P", state_cur.p_h, "T", state_cur.T_h, medium_hot_in.mediumName);
    state_cur.k_c := CP.PropsSI("L", "P", state_cur.p_c, "T", state_cur.T_c, medium_cool_in.mediumName);  
    
    state_cur.Re_h := G_h * d_h / state_cur.mu_h;
    state_cur.Re_c := G_c * d_h / state_cur.mu_c; 
    
    state_cur.rho_h := CP.PropsSI("D", "P", state_cur.p_h, "T", state_cur.T_h, medium_hot_in.mediumName);
    state_cur.rho_c := CP.PropsSI("D", "P", state_cur.p_c, "T", state_cur.T_c, medium_cool_in.mediumName);   
    
    state_cur.u_h := inlet_hot.m_flow / A_stack / state_cur.rho_h;   
    state_cur.u_c := inlet_cool.m_flow / A_stack / state_cur.rho_c; 
    
    state_cur.Nu_h := 4.089 + fit_const_c * (state_cur.Re_h^ fit_const_d);
    state_cur.Nu_c := 4.089 + fit_const_c * (state_cur.Re_c^ fit_const_d);
    
    state_cur.k_wall := thermal_conductivity(name_material, (state_cur.T_h + state_cur.T_c) / 2);
    
    state_cur.h_h := state_cur.Nu_h * state_cur.k_h / d_h;
    state_cur.h_c := state_cur.Nu_c * state_cur.k_c / d_h;
    
    state_cur.U := 1 / ( 1 / state_cur.h_h + 1 / state_cur.h_c + t_wall / state_cur.k_wall);
    
    state_cur.f_h := (15.78 + fit_const_a * state_cur.Re_h ^ fit_const_b ) / state_cur.Re_h;
    state_cur.f_c := (15.78 + fit_const_a * state_cur.Re_c ^ fit_const_b ) / state_cur.Re_c;
    
    if state_cur.T_h > state_cur.T_c then
      state_cur.q := state_cur.U * A_stack * (state_cur.T_h - state_cur.T_c);      
    else
      state_cur.q := 0;
    end if;
    
    state_cur.dp_h := state_cur.f_h * state_cur.length * state_cur.rho_h *  (state_cur.u_h ^ 2) / d_h;
    state_cur.dp_c := state_cur.f_c * state_cur.length * state_cur.rho_c *  (state_cur.u_c ^ 2) / d_c;
    
    // Calculate state for next cell
    state_next.h_mass_h := (state_cur.h_mass_h * m_dot_hot - state_cur.q) / m_dot_hot;
    state_next.h_mass_c := (state_cur.h_mass_c * m_dot_cool - state_cur.q) / m_dot_cool;
    
    state_next.p_h := state_cur.p_h - state_cur.dp_h;
    state_next.p_c := state_cur.p_c - state_cur.dp_c;
    
    state_next.T_h := CP.PropsSI("T", "P",state_next.p_h, "H", state_next.h_mass_h, medium_hot_in.mediumName);
    state_next.T_c := CP.PropsSI("T", "P",state_next.p_c, "H", state_next.h_mass_c, medium_cool_in.mediumName);
    
  end for;
  
equation 
  
  // for out_cool
  medium_cool_out.state = PBMedia.setState_pTX(p = state_cell[N_seg].p_h, T = state_cell[N_seg].T_h);
  outlet_cool.T = medium_cool_out.T;
  outlet_cool.p = medium_cool_out.p;
  outlet_cool.m_flow + inlet_cool.m_flow = 0;
  outlet_cool.h_outflow = medium_cool_out.h;  
  inlet_cool.h_outflow = inStream(outlet_cool.h_outflow);
  
  // for outlet_hot = outlet_main
  medium_hot_out.state = PBMedia.setState_pTX(p = state_cell[N_seg].p_c, T = state_cell[N_seg].T_c); 
  outlet_hot.T = medium_hot_out.T;
  outlet_hot.p = medium_hot_out.p;
  outlet_hot.m_flow + inlet_hot.m_flow = 0;
  outlet_hot.h_outflow = medium_hot_out.h;  
  inlet_hot.h_outflow = inStream(outlet_hot.h_outflow);

end PCHeatExchanger;
