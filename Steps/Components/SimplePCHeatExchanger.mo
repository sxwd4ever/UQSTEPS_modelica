within Steps.Components;

model SimplePCHeatExchanger
  "Simple Printed Circuit based Heat Exchanger for verification purpose"
  //extends BaseExchanger;
  import CP = Steps.Utilities.CoolProp; 
  import TB = Modelica.Blocks.Tables;  
  import UTIL = Modelica.Utilities;
  import MyUtil = Steps.Utilities.Util;
  
  replaceable Steps.Interfaces.PBFluidPort_a inlet_hot(redeclare package Medium = PBMedia, p(start= 1e6)) "Inlet port, previous component";
  replaceable Steps.Interfaces.PBFluidPort_b outlet_hot(redeclare package Medium = PBMedia, p(start= 1e6)) "Outlet port, next component";
  replaceable Steps.Interfaces.PBFluidPort_a inlet_cool(redeclare package Medium = PBMedia, p(start= 1e6)) "Recuperator inlet";
  replaceable Steps.Interfaces.PBFluidPort_b outlet_cool(redeclare package Medium = PBMedia, p(start= 1e6)) "Recuperator outlet";
  
  replaceable package PBMedia = Steps.Media.SCO2;  
  
  parameter Integer N_seg = 6 "Number of segments in a tube";
 
  inner parameter Modelica.SIunits.Length length_cell = 1e-3 "length of a cell";
  
  parameter String name_material = "inconel 750";  
  
  parameter Modelica.SIunits.Angle phi = 0.0 "unit rad";
  
  inner parameter Modelica.SIunits.ReynoldsNumber Re_design = 1200 "On-design ReynoldsNumber";
  
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
  
  parameter Boolean debug_mode = false;
  //protected
  inner KimCorrelations kim_cor(phi = phi, pitch = pitch, d_h = d_h);  

  MaterialConductivity mc(name_material = name_material);
    
  inner Modelica.SIunits.Diameter d_h = 4 * A_c / peri_c "Hydraulic Diameter";
  
  inner Modelica.SIunits.Length peri_c = d_c * Modelica.Constants.pi /2 + d_c "perimeter of semi-circular";
  
  inner Modelica.SIunits.Length t_wall = (2 - Modelica.Constants.pi  / 4) * (d_c / 2) "thickness of wall between two neighboring hot and cold";
  
  inner Integer N_channel = integer(A_fmax / A_c) "number of channels";
  
  inner Modelica.SIunits.Area A_c = Modelica.Constants.pi * d_c * d_c / 8 "Area of semi-circular tube";    
  
  inner Modelica.SIunits.Area A_flow = N_channel * A_c "Flow area of all channels";
  
  inner Modelica.SIunits.Area A_fc = m_dot_cool * d_h / mu_c / Re_design "Area of cold stream area";
  inner Modelica.SIunits.Area A_fh = m_dot_hot * d_h / mu_h /Re_design "Area of hot stream area";
  inner Modelica.SIunits.Area A_fmax = max(A_fc, A_fh) "Area of maximum stream area comparing A_fc and A_fh: A_fmax = max(A_fc, A_fh)";
  
  inner Modelica.SIunits.Area A_stack = peri_c * length_cell * N_channel "surface area of all cells in a stack";
  
  Modelica.SIunits.DynamicViscosity mu_c = CP.PropsSI("V", "P", p_cool, "T", T_cool_in, PBMedia.mediumName) "average dynamic Viscosity in cold channel";
  Modelica.SIunits.DynamicViscosity mu_h = CP.PropsSI("V", "P", p_hot, "T", T_hot_in, PBMedia.mediumName) "average dynamic Viscosity in hot channel";   
  
  Modelica.SIunits.Length length_ch = length_cell * N_seg "length of one pipe in HeatExchanger unit m";
    
  String hot_stream_name, cold_stream_name;

  // two sequences of the hx cells
  HXCell [N_seg] cell_cold;
  HXCell [N_seg] cell_hot;

  model HXCell    
    "One cell in a HX segment"
    
    replaceable package PBMedia = Steps.Media.SCO2; 
    
    replaceable Steps.Interfaces.PBFluidPort_a inlet(redeclare package Medium = Steps.Media.SCO2) "Inlet port, previous component";
    replaceable Steps.Interfaces.PBFluidPort_b outlet(redeclare package Medium = Steps.Media.SCO2) "Outlet port, next component";    
   
    parameter Boolean debug_mode = true;
    
    // Heat Flux
    Modelica.SIunits.Heat Q; 
    
    // mass flux
    Modelica.SIunits.MassFlowRate G;      
    
//  protected
    outer Modelica.SIunits.Diameter d_h "Hydraulic Diameter";  
    
    outer KimCorrelations kim_cor;      
   
    outer Modelica.SIunits.Area A_flow "Flow area of all channels";
    
    outer Modelica.SIunits.Area A_c "Area of semi-circular tube"; 
    
    // length of this cell
    outer Modelica.SIunits.Length length_cell "unit m";  
    
    outer Integer N_channel "number of channels";
    
    //Local temperature
    Modelica.SIunits.Temperature T;
    
    //Local pressure
    Modelica.SIunits.Pressure p;          
  
  equation   
    
    p = (inlet.p + outlet.p) / 2;
    //p = inlet.p;    
    /*
    if inlet.m_flow > 0 then
      T = CP.PropsSI("T", "P", inlet.p, "H", inStream(inlet.h_outflow), PBMedia.mediumName);
    else
      T = CP.PropsSI("T", "P", outlet.p, "H", inStream(outlet.h_outflow), PBMedia.mediumName);
    end if;
    */
    T = (CP.PropsSI("T", "P", inlet.p, "H", inStream(inlet.h_outflow), PBMedia.mediumName) + 
        CP.PropsSI("T", "P", outlet.p, "H", inStream(outlet.h_outflow), PBMedia.mediumName)) / 2 ;  
    
    // mass balance
    inlet.m_flow + outlet.m_flow = 0;
    
    // energy balance
    (outlet.h_outflow - inlet.h_outflow) * inlet.m_flow = Q;     
    
    //outlet.h_outflow * outlet.m_flow = inlet.h_outflow * inlet.m_flow;
    //inlet.h_outflow = inStream(outlet.h_outflow);
    
    //pressure drop
    outlet.p - inlet.p = 0;
      
  end HXCell;
algorithm  
  // initialize the state for the hot side and cold side in all the cells
// using the hot inlet and cool outlet.
  hot_stream_name := PBMedia.mediumName;
  cold_stream_name := PBMedia.mediumName; 

equation  
  // connect all the segments within the heat exchanger, except for the end segment
  for i in 1 : N_seg loop
  
    if i <> N_seg then
      // connect current segment's cool outlet with next segment's cool inlet
      connect(cell_cold[i].outlet, cell_cold[i+1].inlet);
      cell_cold[i+1].inlet.h_outflow = cell_cold[i].outlet.h_outflow;
    end if;
    
    if i <> 1 then
    // connect current segment's hot outlet with previous segment's hot inlet
      connect(cell_hot[i].outlet, cell_hot[i-1].inlet);
      cell_hot[i].outlet.h_outflow = cell_hot[i -1].inlet.h_outflow;
    end if;
        
  end for; 
  
  for i in 1 : N_seg loop
    
    cell_cold[i].G = inlet_cool.m_flow / N_channel / A_c;
    cell_hot[i].G = inlet_hot.m_flow / N_channel / A_c;    
    
    cell_cold[i].Q = 20;
    cell_hot[i].Q = -20;  
  
  end for;
  
  // Now connect the end segement with my inlet and outlet
  connect(inlet_cool, cell_cold[1].inlet);
  connect(cell_cold[N_seg].outlet, outlet_cool);  
  cell_hot[N_seg].inlet.h_outflow = inStream(inlet_hot.h_outflow);
  cell_cold[N_seg].outlet.h_outflow = inStream(outlet_cool.h_outflow);   
    
  connect(outlet_hot, cell_hot[1].outlet);    
  connect(cell_hot[N_seg].inlet, inlet_hot);    
  cell_hot[1].outlet.h_outflow = inStream(outlet_hot.h_outflow);
  cell_cold[1].inlet.h_outflow = inStream(inlet_cool.h_outflow);
  
end SimplePCHeatExchanger;
