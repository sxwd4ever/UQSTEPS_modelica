within Steps.Components;


model PCHeatExchanger
  "Printed Circuit based Heat Exchanger"
  //extends BaseExchanger; 
  
  import Steps.Utilities.CoolProp.PropsSI; 
  import CP = Steps.Utilities.CoolProp; 
  import TB = Modelica.Blocks.Tables;  
  import UTIL = Modelica.Utilities;
  import MyUtil = Steps.Utilities.Util;
  import Steps.Utilities.CoolPropExternalObject;
  import Steps.Utilities.CoolProp.QueryProps;
  import Modelica.SIunits.Conversions.{from_degC, from_bar};
  import Steps.Model.HEBoundaryCondition;
  
/*  
  replaceable Steps.Interfaces.PBFluidPort_a inlet_hot(redeclare package Medium = PBMedia, p(start= p_start_hot), h_outflow(start = h_start_hot)) "Inlet port, previous component";
  replaceable Steps.Interfaces.PBFluidPort_b outlet_hot(redeclare package Medium = PBMedia, p(start= p_start_hot), h_outflow(start = h_start_hot)) "Outlet port, next component";
  replaceable Steps.Interfaces.PBFluidPort_a inlet_cold(redeclare package Medium = PBMedia, p(start= p_start_cold), h_outflow(start = h_start_cold)) "Recuperator inlet";
  replaceable Steps.Interfaces.PBFluidPort_b outlet_cold(redeclare package Medium = PBMedia, p(start= p_start_cold), h_outflow(start = h_start_cold)) "Recuperator outlet";
*/
  
  replaceable Steps.Interfaces.PBFluidPort_a inlet_hot(redeclare package Medium = PBMedia, p(start=bc.st_hot_in.p), h_outflow(start=bc.st_hot_in.h), m_flow(start=bc.st_hot_in.mdot)) "Inlet port, previous component";
  replaceable Steps.Interfaces.PBFluidPort_b outlet_hot(redeclare package Medium = PBMedia, p(start=bc.st_hot_out.p), h_outflow(start=bc.st_hot_out.h), m_flow(start=bc.st_hot_out.mdot)) "Outlet port, next component";
  replaceable Steps.Interfaces.PBFluidPort_a inlet_cold(redeclare package Medium = PBMedia, p(start=bc.st_cold_in.p), h_outflow(start=bc.st_cold_in.h), m_flow(start=bc.st_cold_in.mdot)) "Recuperator inlet";
  replaceable Steps.Interfaces.PBFluidPort_b outlet_cold(redeclare package Medium = PBMedia, p(start=bc.st_cold_out.p), h_outflow(start=bc.st_cold_out.h), m_flow(start=bc.st_cold_out.mdot)) "Recuperator outlet";  

  replaceable package PBMedia = Steps.Media.SCO2; 
  
  parameter String name = "PCHE";
  
  parameter String name_material = "inconel 750";    
  
  MaterialConductivity mc(name_material = name_material);
  // On design parameters
  // Geometry parameters  
  inner parameter PCHEGeoParam geo(
    // pitch length
    pitch = 12e-3,
    // pitch angle
    phi = from_deg((180 - 108) /2),
    // length of pche, mm
    length = 2860e-3,
    // Diameter of semi_circular
    d_c = 2e-3,
    // number of channels
    N_ch = integer(80e3),
    // number of segments
    N_seg = 50 // [Kwon2019]'s maximum node number
  ); 
  
  parameter Modelica.SIunits.ReynoldsNumber Re_hot_start = 5e3 "Re off design value in hot stream";

  parameter Modelica.SIunits.ReynoldsNumber Re_cold_start = 5e3 "Re off design value in cold stream";  
  
  // Geometry determined correlation coefficients - a, b, c d
  inner KimCorrelations kim_cor(phi = geo.phi, pitch = geo.pitch, d_h = d_h); 
  
  inner Modelica.SIunits.Length length_cell = geo.length / geo.N_seg "length of a cell";  
  
  inner CoolPropExternalObject cp_wrapper(PBMedia.mediumName, name);
  
  // d_c determined variables, d_h, A_c, peri_c
  inner Modelica.SIunits.Diameter d_h = 4 * A_c / peri_c "Hydraulic Diameter";  
  inner Modelica.SIunits.Area A_c = Modelica.Constants.pi * geo.d_c * geo.d_c / 8 "Area of semi-circular tube"; 
  inner Modelica.SIunits.Length peri_c = geo.d_c * Modelica.Constants.pi / 2 + geo.d_c "perimeter of semi-circular"; 
  inner Modelica.SIunits.Length t_wall = (2 - Modelica.Constants.pi  / 4) * (geo.d_c / 2) "thickness of wall between two neighboring hot and cold"; 
  inner Modelica.SIunits.Area A_stack = peri_c * length_cell * geo.N_ch "surface area of all cells in a stack";
  inner Modelica.SIunits.Area A_flow = geo.N_ch * A_c "Flow area for all channels";
  
  parameter Boolean ByInlet_hot = false "flag indicate if the fluid state is determined by upstream, which can acelerate the convergence of simulation. "; 
  
  parameter Boolean ByInlet_cold = false "flag indicate if inlet states is fixed 1: fixed; 0: free (cold stream is determined by outlet)"; 
  
  parameter HEBoundaryCondition bc(
    st_hot_in(p = from_bar(80), T = from_degC(578.22), h = PropsSI("H", "P", bc.st_hot_in.p, "T", bc.st_hot_in.T, PBMedia.mediumName), mdot = 51.91),   
    st_cold_in(p = from_bar(200), T = from_degC(151.45), h = PropsSI("H", "P", bc.st_cold_in.p, "T", bc.st_cold_in.T, PBMedia.mediumName), mdot = 51.91), 
    st_hot_out(p = from_bar(80), T = from_degC(156.5), h = PropsSI("H", "P", bc.st_hot_out.p, "T", bc.st_hot_out.T, PBMedia.mediumName), mdot = 51.91),
    st_cold_out(p = from_bar(200), T = from_degC(533.5), h = PropsSI("H", "P", bc.st_cold_out.p, "T", bc.st_cold_out.T, PBMedia.mediumName), mdot = 51.91)    
  );  
  
  // output for debug purpose
  Modelica.SIunits.Temp_C T_hot_in = to_degC(cell_hot[1].T);  
  Modelica.SIunits.Temp_C T_cold_in = to_degC(cell_cold[geo.N_seg].T);  
  Modelica.SIunits.Temp_C T_hot_out = to_degC(cell_hot[geo.N_seg].T);  
  Modelica.SIunits.Temp_C T_cold_out = to_degC(cell_cold[1].T);  
  
  // two sequences of the hx cells
  // use id to distinguish cold and hot cell
  // hot_cell.id > 1000 and < 2000
  // cold_cell.id > 2000
  HXCell [geo.N_seg] cell_cold(
    each ByInlet = ByInlet_cold, 
    each inlet.p.start = bc.st_cold_in.p, 
    each inlet.h_outflow.start = bc.st_cold_in.h,
    //each outlet.p.start = p_start_cold,
    //each outlet.h_outflow.start = CP.PropsSI("H", "P", p_start_cold, "T", T_start_cold, PBMedia.mediumName),
    each T.start = bc.st_cold_in.T, 
    each Re.start = Re_cold_start,  
    //each inlet.m_flow.start = mdot_start_cold,  
    id = {i + 2000 for i in 1 : geo.N_seg}); 
    
  HXCell [geo.N_seg] cell_hot(
    each ByInlet = ByInlet_hot, 
    each inlet.p.start = bc.st_hot_in.p,
    each inlet.h_outflow.start = bc.st_hot_in.h,  
    //each outlet.p.start = p_start_hot, 
    //each outlet.h_outflow.start = CP.PropsSI("H", "P", p_start_hot, "T", T_start_hot, PBMedia.mediumName),
    each T.start = bc.st_hot_in.T,     
    each Re.start = Re_hot_start,
    //each inlet.m_flow.start = mdot_start_hot,
    id = {i + 1000 for i in 1 : geo.N_seg});

  // Heat Change
  Modelica.SIunits.Heat Q[geo.N_seg];  
  
  // wall thermal conductivity - determined by material of wall and local temperature
  Modelica.SIunits.ThermalConductivity k_wall[geo.N_seg];
  
  // overall Heat transfer coefficient
  Modelica.SIunits.CoefficientOfHeatTransfer U[geo.N_seg];   
  
protected
  // internal variable, for debug or efficient purpose
  parameter Boolean debug_mode = false;   
  
  //String hot_stream_name, cold_stream_name;
  
  parameter Real testVal = 1; 

  model HXCell    
    "One cell in a HX segment"
    
    replaceable package PBMedia = Steps.Media.SCO2; 
    
    replaceable Steps.Interfaces.PBFluidPort_a inlet(redeclare package Medium = Steps.Media.SCO2) "Inlet port, previous component";
    replaceable Steps.Interfaces.PBFluidPort_b outlet(redeclare package Medium = Steps.Media.SCO2) "Outlet port, next component";    
   
    parameter Boolean debug_mode = false;
    
    parameter Boolean ByInlet = true;
    
    parameter Integer id = 0;
    
    PBMedia.ThermodynamicState state;
 
    // Heat Flux
    Modelica.SIunits.Heat Q; 
    
    // mass flux
    Modelica.SIunits.MassFlowRate G;      
    
    outer Modelica.SIunits.Diameter d_h "Hydraulic Diameter";  
    
    outer KimCorrelations kim_cor; 
        
    outer Modelica.SIunits.Area A_flow;
    
    outer Modelica.SIunits.Length length_cell;
    
    outer PCHEGeoParam geo;
    
    outer CoolPropExternalObject cp_wrapper;
    
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
    Modelica.SIunits.ReynoldsNumber Re;
    
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

    // use upstream or downstream state can increase speed of convergence
    if ByInlet then
      p = inlet.p;
      h = inlet.h_outflow;  
                                  
    else
      p = outlet.p;
      h = outlet.h_outflow; 
    end if;  
    
    // use average value is more flexible, but difficult to find solution  
/*
algorithm
      
    p := (outlet.p + inlet.p) / 2;
    h := (outlet.h_outflow + inlet.h_outflow) / 2; 
*/

equation
    
    inlet.h_outflow = inStream(inlet.h_outflow); // set up equation between inlet and previous component's outlet 

    // mass balance
    inlet.m_flow + outlet.m_flow = 0;
    
    // energy balance  
    (outlet.h_outflow - inlet.h_outflow) * inlet.m_flow = Q;
       
    inlet.p - outlet.p = dp; 
    
    state = PBMedia.setState_ph(p,h);
    
    T = PBMedia.temperature_ph(p, h);  
    mu = PBMedia.dynamicViscosity(state);
    k = PBMedia.thermalConductivity(state);
    rho = PBMedia.density(state);
    
     
    //(T, mu, k, rho) = CP.MyPropsSI(p=p, H=h, fluidName=PBMedia.mediumName);
    
    /*
    T = CP.PropsSI("V", "P", p, "H", h, PBMedia.mediumName); 
        
    mu = CP.PropsSI("V", "P", p, "T", T, PBMedia.mediumName); 

    k = CP.PropsSI("L", "P", p, "T", T, PBMedia.mediumName);   
    
    rho = CP.PropsSI("D", "P", p, "T", T, PBMedia.mediumName);     
    
    
    T = QueryProps(cp_wrapper, "HmassP_INPUTS", h, p, "T"); // ("V", "P", p, "H", h, PBMedia.mediumName); 
        
    mu = QueryProps(cp_wrapper, "HmassP_INPUTS", h, p, "VISCOSITY");

    k = QueryProps(cp_wrapper, "HmassP_INPUTS", h, p, "CONDUCTIVITY");  
    
    rho = QueryProps(cp_wrapper, "HmassP_INPUTS", h, p, "DMASS");   
    */
    
    Re = G * d_h / mu; 
    
    Nu = 4.089 + kim_cor.c * (Re ^ kim_cor.d);    
    
    u = inlet.m_flow / A_flow / rho;     
    
    hc = Nu * k / d_h;
    
    f = (15.78 + kim_cor.a * Re ^ kim_cor.b ) / Re;  
    
    //pressure drop, unit Pa  
    dp = 2 * f * length_cell * rho *  (u ^ 2) / d_h;     
    
    
/*        
algorithm 
    // ******************************************** 
    // this algorithm section is used for debug purpose
    // if a variable become invalid (zero, inf or weired)
    // 1. make a copy of its equation
    // 2. move the copy into this section, rewrite the equation into an algorithm ('=' -> ':=')
    // 3. commet the origin one
    // 4. update the calling of function MyAseert accordingly to print value in console
    //
    // once the bug fixed, DO REMEMBER to rewind the changed lines reversely. 
    // ********************************************
    //Re := G * d_h / mu; 
    
    MyUtil.myAssertNotEqual(
    debug = false, 
    val_test = Re, compared = 0,
    name_val = "Re", 
    val_ref = {id, Nu, Re, G, f, d_h, mu, T, p ,  h, rho, u, length_cell, dp}, 
    name_val_ref = {"id", "Nu", "Re", "G", "f", "d_h", "mu",  "T", "p" , "h", "rho", "u", "length_cell", "dp"}); 
    
    //k := CP.PropsSI("L", "P", p, "T", T, PBMedia.mediumName);    
      
    MyUtil.myAssert(
    debug = false, 
    val_test = k, min = 0, max = 1e5,
    name_val = "k", 
    val_ref = {id, Nu, Re, G, f, d_h, mu, T, p ,  h, rho, u, length_cell, dp}, 
    name_val_ref = {"id", "Nu", "Re", "G", "f", "d_h", "mu",  "T", "p" , "h", "rho", "u", "length_cell", "dp"});    
    
    //Nu := 4.089 + kim_cor.c * (Re ^ kim_cor.d);    
    
    MyUtil.myAssert(
    debug = false, 
    val_test = Nu, min = 0, max = 1e5,
    name_val = "Nu", 
    val_ref = {id, Nu, Re, G, f, d_h, mu, p , T, h, rho, u, length_cell, dp}, 
    name_val_ref = {"id", "Nu", "Re", "G", "f", "d_h", "mu", "p" , "T", "h", "rho", "u", "length_cell", "dp"});  
*/      
  end HXCell;
/* 
algorithm

  cold_stream_name := PBMedia.mediumName; 
  hot_stream_name := PBMedia.mediumName;
*/
equation  
  
  // connect all the segments within the heat exchanger, except for the end segment
  for i in 1 : geo.N_seg loop
  
    if i <> 1 then
      // connect current segment's cold outlet with next segment's cold inlet
      connect(cell_cold[i].outlet, cell_cold[i-1].inlet);
    end if;
    
    if i <> geo.N_seg then
      // connect current segment's hot outlet with previous segment's hot inlet
      connect(cell_hot[i].outlet, cell_hot[i+1].inlet);
    end if;
        
  end for; 
  
  // Now connect the end segement with my inlet and outlet
  connect(cell_cold[1].outlet, outlet_cold);
  connect(inlet_cold, cell_cold[geo.N_seg].inlet); 
 
  connect(cell_hot[1].inlet, inlet_hot);   
  connect(outlet_hot, cell_hot[geo.N_seg].outlet);    
  
  //inlet_cold.h_outflow = inStream(inlet_cold.h_outflow);   
  //inlet_hot.h_outflow = inStream(inlet_hot.h_outflow);    


  for i in 1 : geo.N_seg loop
  
    k_wall[i] = MyUtil.thermal_conductivity(tableID = mc.table_th_inconel_750, name = name_material, temperature = (cell_cold[i].T + cell_hot[i].T) / 2);
 
     U[i] =  1 / (1 / cell_hot[i].hc + 1 / cell_cold[i].hc + t_wall / k_wall[i]);   
    
    if cell_hot[i].T > cell_cold[i].T then
      Q[i] = U[i] * A_stack * (cell_hot[i].T - cell_cold[i].T);      
    else
      Q[i] = 0;
    end if;
        
    cell_cold[i].Q = Q[i];
    cell_hot[i].Q = -Q[i];  
    
    cell_cold[i].G = inlet_cold.m_flow / geo.N_ch / A_c;
    cell_hot[i].G = inlet_hot.m_flow / geo.N_ch / A_c;   
  
  end for;

end PCHeatExchanger;
