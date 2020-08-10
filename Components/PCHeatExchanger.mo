within Steps.Components;

model PCHeatExchanger
  "Printed Circuit based Heat Exchanger"
  extends BaseExchanger;  
  
  import CP = Steps.Utilities.CoolProp; 
  import TB = Modelica.Blocks.Tables;  
  import UTIL = Modelica.Utilities;
  import MyUtil = Steps.Utilities.Util;
  
  replaceable package PBMedia = Steps.Media.SCO2;  
  
  Modelica.Blocks.Types.ExternalCombiTable1D table_th_inconel_750 = Modelica.Blocks.Types.ExternalCombiTable1D(tableName = "inconel_750", fileName = Modelica.Utilities.Files.loadResource("modelica://Steps/Resources/Data/th_conductivity.txt"), table = fill(0.0, 6, 2), smoothness = Modelica.Blocks.Types.Smoothness.LinearSegments, columns = 2:2) "thermal conductivity for inconel_750";	
  
  parameter Integer N_seg = 10 "Number of segments in a tube";
  
  parameter Modelica.SIunits.Length length_cell = 1e-3 "length of a cell";
  
  parameter Modelica.SIunits.Angle phi = 0.0 "unit rad";
  
  parameter Modelica.SIunits.ReynoldsNumber Re_design = 1200 "On-design ReynoldsNumber";
  
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
  
  parameter Boolean debug_mode = false;
  
protected  
  
  inner KimCorrelations kim_cor(phi = phi, pitch = pitch, d_h = d_h);  
  
  Modelica.SIunits.Diameter d_h = 4 * A_c / peri_c "Hydraulic Diameter";
  
  Modelica.SIunits.Length peri_c = d_c * Modelica.Constants.pi /2 + d_c "perimeter of semi-circular";
  
  Modelica.SIunits.Length t_wall = (2 - Modelica.Constants.pi  / 4) * (d_c / 2) "thickness of wall between two neighboring hot and cold";
  
  Integer N_channel = integer(A_fmax / A_c) "number of channels";
  
  Modelica.SIunits.Area A_c = Modelica.Constants.pi * d_c * d_c / 8 "Area of semi-circular tube";    
  
  Modelica.SIunits.Area A_flow = N_channel * A_c "Flow area of all channels";
  
  Modelica.SIunits.Area A_fc = m_dot_cool * d_h / mu_c / Re_design "Area of cold stream area";
  Modelica.SIunits.Area A_fh = m_dot_hot * d_h / mu_h /Re_design "Area of hot stream area";
  Modelica.SIunits.Area A_fmax = max(A_fc, A_fh) "Area of maximum stream area comparing A_fc and A_fh: A_fmax = max(A_fc, A_fh)";
  
  Modelica.SIunits.Area A_stack = peri_c * length_cell * N_channel "surface area of all cells in a stack";
  
  Modelica.SIunits.DynamicViscosity mu_c = CP.PropsSI("V", "P", p_cool, "T", T_cool_in, PBMedia.mediumName) "average dynamic Viscosity in cold channel";
  Modelica.SIunits.DynamicViscosity mu_h = CP.PropsSI("V", "P", p_hot, "T", T_hot_in, PBMedia.mediumName) "average dynamic Viscosity in hot channel"; 
  
  HXCellState [N_seg] state_cell "state array of each cells";
  
  Modelica.SIunits.Length length_ch = length_cell * N_seg "length of one pipe in HeatExchanger unit m";

  Modelica.SIunits.MassFlowRate G_c "Mass flux of each cold channel";  
  Modelica.SIunits.MassFlowRate G_h "Mass flux of each hot channel";    
     
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
  
  state_cell[1].h_c := medium_cool_out.h;
  state_cell[1].h_h := medium_hot_in.h;  
  
  G_c := state_cell[1].mdot_c / N_channel / A_c;
  G_h := state_cell[1].mdot_h / N_channel / A_c;

  for i in 1 : N_seg loop // index out of range??
    state_cell[i].length := length_cell;
   
    //Debug from this point
    state_cell[i].mu_c := CP.PropsSI("V", "P", state_cell[i].p_c, "T", state_cell[i].T_c, state_cell[i].medium_name_c);
    state_cell[i].mu_h := CP.PropsSI("V", "P", state_cell[i].p_h, "T", state_cell[i].T_h, state_cell[i].medium_name_h);
      
    state_cell[i].k_c := CP.PropsSI("L", "P", state_cell[i].p_c, "T", state_cell[i].T_c, state_cell[i].medium_name_c);  
    
    MyUtil.myAssert(debug = debug_mode, val_test = state_cell[i].k_c, min = 0, max = 1e5, name_val = "k_c", val_ref = {i, state_cell[i].T_c, state_cell[i].p_c}, name_val_ref = {"i", "T", "P"});    

      
    state_cell[i].k_h := CP.PropsSI("L", "P", state_cell[i].p_h, "T", state_cell[i].T_h, state_cell[i].medium_name_h);  
      
    state_cell[i].Re_c := G_c * d_h / state_cell[i].mu_c;
    state_cell[i].Re_h := G_h * d_h / state_cell[i].mu_h;
      
    state_cell[i].rho_c := CP.PropsSI("D", "P", state_cell[i].p_c, "T", state_cell[i].T_c, state_cell[i].medium_name_c);
    state_cell[i].rho_h := CP.PropsSI("D", "P", state_cell[i].p_h, "T", state_cell[i].T_h, state_cell[i].medium_name_h);
      
    state_cell[i].u_c := state_cell[i].mdot_c / A_flow / state_cell[i].rho_c;  
    state_cell[i].u_h := state_cell[i].mdot_h / A_flow / state_cell[i].rho_h;     
    
    
    MyUtil.myAssert(debug = debug_mode, val_test = state_cell[i].Re_c, min = 0, max = 1e6, name_val = "Re_c", 
    name_val_ref = {"i","G","d_h","mu"},
    val_ref = {i, G_c, d_h, state_cell[i].mu_c});
        
    state_cell[i].Nu_c := 4.089 + kim_cor.c * (state_cell[i].Re_c ^ kim_cor.d);
    
    MyUtil.myAssert(debug = debug_mode, val_test = state_cell[i].Re_h, min = 0, max = 1e6, name_val = "Re_h", 
    name_val_ref = {"i","G","d_h","mu"},
    val_ref = {i, G_h, d_h, state_cell[i].mu_h});      
    
    state_cell[i].Nu_h := 4.089 + kim_cor.c * (state_cell[i].Re_h ^ kim_cor.d);
      
    state_cell[i].hc_c := state_cell[i].Nu_c * state_cell[i].k_c / d_h;
    state_cell[i].hc_h := state_cell[i].Nu_h * state_cell[i].k_h / d_h;
      
    state_cell[i].f_c := (15.78 + kim_cor.a * state_cell[i].Re_c ^ kim_cor.b ) / state_cell[i].Re_c;      
    state_cell[i].f_h := (15.78 + kim_cor.a * state_cell[i].Re_h ^ kim_cor.b ) / state_cell[i].Re_h;      
    
    state_cell[i].dp_c := state_cell[i].f_c * state_cell[i].length * state_cell[i].rho_c *  (state_cell[i].u_c ^ 2) / d_h;    
    state_cell[i].dp_h := state_cell[i].f_h * state_cell[i].length * state_cell[i].rho_h *  (state_cell[i].u_h ^ 2) / d_h;
    
    // no use of following parameters, use default value
    state_cell[i].Pr_c := 1.0;
    state_cell[i].Pr_h := 1.0;
    
    state_cell[i].k_wall := MyUtil.thermal_conductivity(tableID = table_th_inconel_750, name = name_material, temperature = (state_cell[i].T_c + state_cell[i].T_h) / 2);
    
    state_cell[i].U := 1 / ( 1 / state_cell[i].hc_h + 1 / state_cell[i].hc_c + t_wall / state_cell[i].k_wall);    

    
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
      
    state_cell[i + 1].h_c := (state_cell[i].h_c * state_cell[i].mdot_c - state_cell[i].q) / state_cell[i].mdot_c;
    /*
    assert(i <> 1,
    "out of range " + keyvalStr("h_c", state_cell[i+1].h_c) + " at " + 
    debugInfo({"i+1","h_c[i]","m_dot_inlet_cold","q[i]"},{i+1,state_cell[i].h_c, state_cell[i].mdot_c, state_cell[i].q})); 
    */
    
    MyUtil.myAssert(debug = debug_mode, val_test = state_cell[i+1].h_c, min = 0, max = 1e6, name_val = "h_c", 
    name_val_ref = {"i+1","h_c[i]","m_dot_inlet_cold","q[i]"},
    val_ref = {i+1,state_cell[i].h_c, state_cell[i].mdot_c, state_cell[i].q});  
    
    state_cell[i + 1].h_h := (state_cell[i].h_h * state_cell[i].mdot_h - state_cell[i].q) / state_cell[i].mdot_h; 
     
    state_cell[i + 1].p_c := state_cell[i].p_c - state_cell[i].dp_c;
    state_cell[i + 1].p_h := state_cell[i].p_h - state_cell[i].dp_h;
      
    state_cell[i + 1].T_c := CP.PropsSI("T", "P",state_cell[i + 1].p_c, "H", state_cell[i + 1].h_c, state_cell[i].medium_name_c);
    
    MyUtil.myAssert(debug = debug_mode, val_test = state_cell[i+1].T_c, min = 0, max = 1e6, name_val = "T_c", 
    name_val_ref = {"i","P","h_c"},
    val_ref = {i,state_cell[i + 1].p_c,state_cell[i + 1].h_c});    
    
    state_cell[i + 1].T_h := CP.PropsSI("T", "P",state_cell[i + 1].p_h, "H", state_cell[i + 1].h_h, state_cell[i].medium_name_h);

  end for;

equation 
  
  medium_cool_in.state = PBMedia.setState_phX(p = state_cell[N_seg].p_c, h = state_cell[N_seg].h_c);   
  inlet_cool.T = medium_cool_in.T;
  inlet_cool.p = medium_cool_in.p;
  inlet_cool.h_outflow = state_cell[N_seg].h_c;  
  //inlet_cool.h_outflow = inStream(outlet_cool.h_outflow);  
  //outlet_cool.h_outflow = - medium_cool_out.h;  
  inlet_cool.h_outflow = - outlet_cool.h_outflow;
  outlet_cool.m_flow + inlet_cool.m_flow = 0;
  
  // for outlet_hot
  medium_hot_out.state = PBMedia.setState_phX(p = state_cell[N_seg].p_h, h = state_cell[N_seg].h_h); 
  outlet_hot.T = medium_hot_out.T;
  outlet_hot.p = medium_hot_out.p;
  outlet_hot.h_outflow = - state_cell[N_seg].h_h;    
  inlet_hot.h_outflow = inStream(outlet_hot.h_outflow);  
  outlet_hot.m_flow + inlet_hot.m_flow = 0;

end PCHeatExchanger;
