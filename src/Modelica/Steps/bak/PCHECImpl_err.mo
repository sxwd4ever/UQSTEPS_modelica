within Steps.Components;


model PCHECImpl
  "Printed Circuit based Heat Exchanger Model implemented in C - this is a wrapper class"
  //extends BaseExchanger; 
  
  import Steps.Utilities.CoolProp.{PropsSI, PCHE_OFFD_Simulation}; 
  import TB = Modelica.Blocks.Tables;  
  import UTIL = Modelica.Utilities;
  import MyUtil = Steps.Utilities.Util;
  import Modelica.SIunits.Conversions.{from_degC, from_bar};
  
  replaceable Steps.Interfaces.PBFluidPort_a inlet_hot(redeclare package Medium = PBMedia, p(start=bc.st_hot_in.p), h_outflow(start=bc.st_hot_in.h), m_flow(start = bc.st_hot_in.mdot)) "Inlet port, previous component";
  replaceable Steps.Interfaces.PBFluidPort_b outlet_hot(redeclare package Medium = PBMedia, p(start=bc.st_hot_out.p), h_outflow(start=bc.st_hot_out.h), m_flow(start = bc.st_hot_out.mdot)) "Outlet port, next component";
  replaceable Steps.Interfaces.PBFluidPort_a inlet_cold(redeclare package Medium = PBMedia, p(start=bc.st_cold_in.p), h_outflow(start=bc.st_cold_in.h), m_flow(start = bc.st_cold_in.mdot)) "Recuperator inlet";
  replaceable Steps.Interfaces.PBFluidPort_b outlet_cold(redeclare package Medium = PBMedia, p(start=bc.st_cold_out.p), h_outflow(start=bc.st_cold_out.h), m_flow(start = bc.st_cold_out.mdot)) "Recuperator outlet"; 

  replaceable package PBMedia = Steps.Media.SCO2; 
  
  parameter String name_material = "inconel 750";    
  
  String name = "PCHE";
  
  MaterialConductivity mc(name_material = name_material);
  // On design parameters
  // Geometry parameters
    
  parameter PCHEGeoParam geo(
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

  // Geometry determined correlation coefficients - a, b, c d
  KimCorrelations kim_cor(phi = geo.phi, pitch = geo.pitch, d_h = d_h);   
  
  KimCorrCoe cor(a=kim_cor.a, b=kim_cor.b, c=kim_cor.c, d=kim_cor.d);  
  
  // following variables are for debug purpose only
  // d_c determined variables, d_h, A_c, peri_c
  Modelica.SIunits.Diameter d_h = 4 * A_c / peri_c "Hydraulic Diameter";
  
  Modelica.SIunits.Area A_c = Modelica.Constants.pi * geo.d_c * geo.d_c / 8 "Area of semi-circular tube";   
  
  Modelica.SIunits.Length peri_c = geo.d_c * Modelica.Constants.pi / 2 + geo.d_c "perimeter of semi-circular";  
  
  Modelica.SIunits.Length t_wall = (2 - Modelica.Constants.pi  / 4) * (geo.d_c / 2) "thickness of wall between two neighboring hot and cold";   
  
  Modelica.SIunits.Area A_stack = peri_c * length_cell * geo.N_ch "surface area of all cells in a stack";
  
  Modelica.SIunits.Area A_flow = geo.N_ch * A_c "Flow area for all channels";   
  
  Modelica.SIunits.Length length_cell = geo.length / geo.N_seg;
  
  // only (p, h, mdot) of the hot/cold inlet matters
  parameter BoundaryCondition bc(
    st_hot_in(p = from_bar(80), T = from_degC(578.22), h = PropsSI("H", "P", bc.st_hot_in.p, "T", bc.st_hot_in.T, PBMedia.mediumName), mdot = 10),   
    st_cold_in(p = from_bar(200), T = from_degC(151.45), h = PropsSI("H", "P", bc.st_cold_in.p, "T", bc.st_cold_in.T, PBMedia.mediumName), mdot = 10), 
    st_hot_out(p = from_bar(80), T = from_degC(156.5), h = PropsSI("H", "P", bc.st_hot_out.p, "T", bc.st_hot_out.T, PBMedia.mediumName), mdot = 10),
    st_cold_out(p = from_bar(200), T = from_degC(533.5), h = PropsSI("H", "P", bc.st_cold_out.p, "T", bc.st_cold_out.T, PBMedia.mediumName), mdot = 10)    
  );      
    
  parameter SimParam sim_param(err=1e-2, delta_T_init = 5, N_iter = 10, step_rel=0.3, log_level = 1);
       
  // output buffer 
  Real h_hot[2];
  Real h_cold[2];
  Real p_hot[2];
  Real p_cold[2];   
  // output for debug purpose
  Modelica.SIunits.Temperature T_hot_in = PropsSI("T", "P", inlet_hot.p, "H", inlet_hot.h_outflow, PBMedia.mediumName);
  
  Modelica.SIunits.Temperature T_cold_in = PropsSI("T", "P", inlet_cold.p, "H", inlet_cold.h_outflow, PBMedia.mediumName);
  
  Modelica.SIunits.Temperature T_hot_out = PropsSI("T", "P", outlet_hot.p, "H", outlet_hot.h_outflow, PBMedia.mediumName);
  
  Modelica.SIunits.Temperature T_cold_out = PropsSI("T", "P", outlet_cold.p, "H", outlet_cold.h_outflow, PBMedia.mediumName);
/*
algorithm
  if initial() then  
    bc.st_hot_in.p := bc_hot_in.p;
    bc.st_hot_in.h := bc_hot_in.h;
    bc.st_hot_in.mdot := bc_hot_in.mdot;
    
    bc.st_cold_in.p := bc_cold_in.p;
    bc.st_cold_in.h := bc_cold_in.h;
    bc.st_cold_in.mdot := bc_cold_in.mdot;
  else
    bc.st_hot_in.p := inlet_hot.p;
    bc.st_hot_in.h := inlet_hot.h_outflow;
    bc.st_hot_in.mdot := inlet_hot.m_flow;
    
    bc.st_cold_in.p := inlet_cold.p;
    bc.st_cold_in.h := inlet_cold.h_outflow;
    bc.st_cold_in.mdot := inlet_cold.m_flow;  
  end if;
*/

equation   
 /* 
  if initial() then  
    bc.st_hot_in.p = bc_hot_in.p;
    bc.st_hot_in.h = bc_hot_in.h;
    bc.st_hot_in.mdot = bc_hot_in.mdot;
    
    bc.st_cold_in.p = bc_cold_in.p;
    bc.st_cold_in.h = bc_cold_in.h;
    bc.st_cold_in.mdot = bc_cold_in.mdot;
  else
    bc.st_hot_in.p = inlet_hot.p;
    bc.st_hot_in.h = inlet_hot.h_outflow;
    bc.st_hot_in.mdot = inlet_hot.m_flow;
    
    bc.st_cold_in.p = inlet_cold.p;
    bc.st_cold_in.h = inlet_cold.h_outflow;
    bc.st_cold_in.mdot = inlet_cold.m_flow;  
  end if;
  */
  bc.st_hot_in.p = inlet_hot.p;
  bc.st_hot_in.h = inlet_hot.h_outflow;
  bc.st_hot_in.mdot = inlet_hot.m_flow;
  
  bc.st_cold_in.p = inlet_cold.p;
  bc.st_cold_in.h = inlet_cold.h_outflow;
  bc.st_cold_in.mdot = inlet_cold.m_flow;  
   
  (h_hot, h_cold, p_hot, p_cold) = PCHE_OFFD_Simulation(name, "CO2", "CO2", geo, cor, sim_param, bc);

  inlet_hot.m_flow + outlet_hot.m_flow = 0;

  inlet_cold.m_flow + outlet_cold.m_flow = 0;
  
  inlet_hot.h_outflow = inStream(inlet_hot.h_outflow);
  
  inlet_cold.h_outflow = inStream(inlet_cold.h_outflow);  


// algorithm
  
  outlet_hot.h_outflow = h_hot[2];
  outlet_hot.p = p_hot[2];
  
  outlet_cold.h_outflow = h_cold[1];
  outlet_cold.p = p_cold[1];
  
end PCHECImpl;
