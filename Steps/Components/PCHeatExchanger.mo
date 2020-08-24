within Steps.Components;

model PCHeatExchanger
  "Printed Circuit based Heat Exchanger"
  //extends BaseExchanger;
  
  import CP = Steps.Utilities.CoolProp; 
  import TB = Modelica.Blocks.Tables;  
  import UTIL = Modelica.Utilities;
  import MyUtil = Steps.Utilities.Util;
  
  replaceable Steps.Interfaces.PBFluidPort_a inlet_hot(redeclare package Medium = PBMedia, p(start= 9e6)) "Inlet port, previous component";
  replaceable Steps.Interfaces.PBFluidPort_b outlet_hot(redeclare package Medium = PBMedia, p(start= 9e6)) "Outlet port, next component";
  replaceable Steps.Interfaces.PBFluidPort_a inlet_cool(redeclare package Medium = PBMedia, p(start= 20e6)) "Recuperator inlet";
  replaceable Steps.Interfaces.PBFluidPort_b outlet_cool(redeclare package Medium = PBMedia, p(start= 20e6)) "Recuperator outlet";
  
  replaceable package PBMedia = Steps.Media.SCO2;  
  
  parameter Integer N_seg = 6 "Number of segments in a tube";
 
  inner parameter Modelica.SIunits.Length length_cell = 1e-3 "length of a cell";
  
  parameter String name_material = "inconel 750";  
  
  inner parameter Modelica.SIunits.Angle phi = 1.0 "unit rad";
  
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
  
  inner parameter Modelica.SIunits.Length pitch = 10 "pitch length of channel";
  
  parameter Boolean debug_mode = false;
  //protected
  inner KimCorrelations kim_cor(phi = phi, pitch = pitch, d_h = d_h);  

  MaterialConductivity mc(name_material = name_material);
    
  inner Modelica.SIunits.Diameter d_h = 4 * A_c / peri_c "Hydraulic Diameter";
  
  inner Modelica.SIunits.Length peri_c = d_c * Modelica.Constants.pi / 2 + d_c "perimeter of semi-circular";
  
  inner Modelica.SIunits.Length t_wall = (2 - Modelica.Constants.pi  / 4) * (d_c / 2) "thickness of wall between two neighboring hot and cold";
  
  // use ceil to avoid integer(0.01) = 0 so that we have one channel at least
  inner Integer N_channel = integer(ceil(A_fmax / A_c)) "number of channels"; 
  
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
  HXCell [N_seg] cell_cold(
    each cellType = CellType.Cold,
    each ByInlet = true, 
    each inlet.p.start = 8e6, 
    each T.start = Modelica.SIunits.Conversions.from_degC(27.15), 
    id = {i for i in 1 : N_seg}); 
    
  HXCell [N_seg] cell_hot(
    each cellType = CellType.Hot,
    each ByInlet = true, 
    each inlet.p.start = 20e6, 
    each T.start = Modelica.SIunits.Conversions.from_degC(700), 
    id = {i for i in 1 : N_seg});

  // Heat Change
  Modelica.SIunits.Heat Q[N_seg];  
  
  // wall thermal conductivity - determined by material of wall and local temperature
  Modelica.SIunits.ThermalConductivity k_wall[N_seg];
  
  // overall Heat transfer coefficient
  Modelica.SIunits.CoefficientOfHeatTransfer U[N_seg];     
  
  parameter Real testVal = 1; 

  type CellType = enumeration(Cold, Hot);

  model HXCell    
    "One cell in a HX segment"
    
    replaceable package PBMedia = Steps.Media.SCO2; 
    
    replaceable Steps.Interfaces.PBFluidPort_a inlet(redeclare package Medium = Steps.Media.SCO2) "Inlet port, previous component";
    replaceable Steps.Interfaces.PBFluidPort_b outlet(redeclare package Medium = Steps.Media.SCO2) "Outlet port, next component";    
   
    parameter Boolean debug_mode = false;
    
    parameter Boolean ByInlet = true;
    
    parameter Integer id = 0;
    
    parameter CellType cellType = CellType.Cold;
    
    // Heat Flux
    Modelica.SIunits.Heat Q; 
    
    // mass flux
    Modelica.SIunits.MassFlowRate G;      
    
//  protected
    outer Modelica.SIunits.Diameter d_h "Hydraulic Diameter";  
    
    outer KimCorrelations kim_cor;      
   
    outer Modelica.SIunits.Area A_flow "Flow area of all channels";
    
    outer Modelica.SIunits.Area A_c "Area of semi-circular tube"; 
    
    outer Modelica.SIunits.Area A_fc;
    
    outer Modelica.SIunits.Area A_fh;
    
    // length of this cell
    outer Modelica.SIunits.Length length_cell "unit m";  
    
    outer Integer N_channel "number of channels";
    
    outer parameter Modelica.SIunits.ReynoldsNumber Re_design "On-design ReynoldsNumber";
    
    outer parameter Modelica.SIunits.Angle phi;
    
    outer parameter Modelica.SIunits.Length pitch;
    
    //Local temperature
    Modelica.SIunits.Temperature T;
    
    //Local pressure
    Modelica.SIunits.Pressure p;
    
    // Local specific Enthalpy
    Modelica.SIunits.SpecificEnthalpy h;
    
    // Local velocity of fluid
    Modelica.SIunits.Velocity u;
    
    //local parameters of this cell listed as following
    //Local Dynamic Viscosity
    Modelica.SIunits.DynamicViscosity mu;
       
    // Local Conductivity
    Modelica.SIunits.ThermalConductivity k;
    
    // Local Reynolds Number
    Modelica.SIunits.ReynoldsNumber Re(start=Re_design);
    
    // Local Density
    Modelica.SIunits.Density rho;
    
    // Local Nusselt Number
    Modelica.SIunits.NusseltNumber Nu;
    
    // Local PrandtlNumber
    //Modelica.SIunits.PrandtlNumber Pr;
    
    // local Thermal Conductance
    Modelica.SIunits.CoefficientOfHeatTransfer hc;
    
    // Fanning Friction Factor - used to calculate pressure drop
    Real f;          
    
    // Local delta p between neighboring cells    
    Modelica.SIunits.PressureDifference dp;    
    
  equation   
    
    if ByInlet then
      p = inlet.p;
      h = inStream(inlet.h_outflow);      
    else
      p = outlet.p;
      h = inStream(outlet.h_outflow);      
    end if;    
    
    T = CP.PropsSI("T", "P", p, "H", h, PBMedia.mediumName);

    //Debug from this point
    mu = CP.PropsSI("V", "P", p, "T", T, PBMedia.mediumName); 
    
    k = CP.PropsSI("L", "P", p, "T", T, PBMedia.mediumName);  
    
    rho = CP.PropsSI("D", "P", p, "T", T, PBMedia.mediumName);     
  
    Re = G * d_h / mu; 
/*    
algorithm    
    // For debug purpose
    MyUtil.myAssertNotEqual(
    debug = false, 
    val_test = id, compared = 10, 
    name_val = "id", 
    val_ref = {Re, G, d_h, mu, p , T, phi, pitch, kim_cor.c, kim_cor.d}, 
    name_val_ref = {"Re", "G", "d_h", "mu", "p" , "T", "phi", "pitch", "c", "d"}); 

equation
*/    

    u = inlet.m_flow / A_flow / rho;      

algorithm

    Nu := 4.089 + kim_cor.c * (Re ^ kim_cor.d);
    /*
    MyUtil.myAssertNotEqual(
    debug = false, 
    val_test = id, compared = 2, 
    name_val = "Re", 
    val_ref = {id, G, d_h, mu, p , T}, 
    name_val_ref = {"id", "G", "d_h", "mu", "p" , "T"}); 
    */
    // For debug purpose
    MyUtil.myAssert(
    debug = false, 
    val_test = T, min = 0, max = 1e5, 
    name_val = "T", 
    val_ref = {id, G, d_h, mu, p , T, h}, 
    name_val_ref = {"id", "G", "d_h", "mu", "p" , "T", "h"});     
    
equation

    hc = Nu * k / d_h;
    
    f = (15.78 + kim_cor.a * Re ^ kim_cor.b ) / Re;  

    // mass balance
    inlet.m_flow + outlet.m_flow = 0;
    
    // energy balance
    (outlet.h_outflow - inlet.h_outflow) * inlet.m_flow = Q; 

    //pressure drop : kPa * 1000 - > pa
    dp = 2 * f * length_cell * rho *  (u ^ 2) / d_h * 1e3;

    inlet.p - outlet.p = dp;
      
  end HXCell;
  
algorithm  
  // initialize the state for the hot side and cold side in all the cells
  // using the hot inlet and cool outlet.
  hot_stream_name := PBMedia.mediumName;
  cold_stream_name := PBMedia.mediumName; 

equation  
  // connect all the segments within the heat exchanger, except for the end segment
  for i in 1 : N_seg loop
  
    if i <> 1 then
      // connect current segment's cool outlet with next segment's cool inlet
      connect(cell_cold[i].outlet, cell_cold[i-1].inlet);
      cell_cold[i-1].inlet.h_outflow = cell_cold[i].outlet.h_outflow;
    end if;
    
    if i <> N_seg then
    // connect current segment's hot outlet with previous segment's hot inlet
      connect(cell_hot[i].outlet, cell_hot[i+1].inlet);
      cell_hot[i].outlet.h_outflow = cell_hot[i + 1].inlet.h_outflow;
    end if;
        
  end for; 
  
  // Now connect the end segement with my inlet and outlet
  connect(cell_cold[1].outlet, outlet_cool);
  connect(inlet_cool, cell_cold[N_seg].inlet);     
  cell_cold[1].outlet.h_outflow = inStream(outlet_cool.h_outflow);
  inlet_cool.h_outflow = inStream(cell_cold[N_seg].inlet.h_outflow);   
  
  connect(cell_hot[1].inlet, inlet_hot);   
  connect(outlet_hot, cell_hot[N_seg].outlet);         
  cell_hot[1].inlet.h_outflow = inStream(inlet_hot.h_outflow);
  cell_hot[N_seg].outlet.h_outflow = inStream(outlet_hot.h_outflow);  
  
/*  
algorithm

  testVal := 1;


  MyUtil.myAssert(
  debug = true, 
  val_test = testVal, min = -1e5, max = -1, 
  name_val = "testVal", 
  val_ref = {inlet_cool.m_flow, N_channel, A_c}, 
  name_val_ref = {"inlet_cool.m_flow", "N_channel", "A_c"}); 

  testVal := inlet_cool.m_flow / N_channel / A_c;

equation
*/  
  for i in 1 : N_seg loop
  
    k_wall[i] = MyUtil.thermal_conductivity(tableID = mc.table_th_inconel_750, name = name_material, temperature = (cell_cold[i].T + cell_hot[i].T) / 2);
 
    1 / U[i] =  1 / cell_hot[i].hc + 1 / cell_cold[i].hc + t_wall / k_wall[i];   
    
    if cell_hot[i].T > cell_cold[i].T then
      Q[i] = U[i] * A_stack * (cell_hot[i].T - cell_cold[i].T);      
    else
      Q[i] = 0;
    end if;
        
    cell_cold[i].Q = Q[i];
    cell_hot[i].Q = -Q[i];  
    
    cell_cold[i].G = inlet_cool.m_flow / N_channel / A_c;
    cell_hot[i].G = inlet_hot.m_flow / N_channel / A_c;   
  
  end for;



end PCHeatExchanger;
