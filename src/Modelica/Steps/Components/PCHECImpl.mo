within Steps.Components;


model PCHECImpl
  "Printed Circuit based Heat Exchanger Model implemented in C - this is a wrapper class"
  //extends BaseExchanger; 
  
  import CP = Steps.Utilities.CoolProp; 
  import TB = Modelica.Blocks.Tables;  
  import UTIL = Modelica.Utilities;
  import MyUtil = Steps.Utilities.Util;
  import Modelica.SIunits.Conversions.{from_degC, from_bar};
/*  
  replaceable Steps.Interfaces.PBFluidPort_a inlet_hot(redeclare package Medium = PBMedia, p(start= p_start_hot), h_outflow(start = h_start_hot)) "Inlet port, previous component";
  replaceable Steps.Interfaces.PBFluidPort_b outlet_hot(redeclare package Medium = PBMedia, p(start= p_start_hot), h_outflow(start = h_start_hot)) "Outlet port, next component";
  replaceable Steps.Interfaces.PBFluidPort_a inlet_cold(redeclare package Medium = PBMedia, p(start= p_start_cold), h_outflow(start = h_start_cold)) "Recuperator inlet";
  replaceable Steps.Interfaces.PBFluidPort_b outlet_cold(redeclare package Medium = PBMedia, p(start= p_start_cold), h_outflow(start = h_start_cold)) "Recuperator outlet";
*/
  
  replaceable Steps.Interfaces.PBFluidPort_a inlet_hot(redeclare package Medium = PBMedia) "Inlet port, previous component";
  replaceable Steps.Interfaces.PBFluidPort_b outlet_hot(redeclare package Medium = PBMedia) "Outlet port, next component";
  replaceable Steps.Interfaces.PBFluidPort_a inlet_cold(redeclare package Medium = PBMedia) "Recuperator inlet";
  replaceable Steps.Interfaces.PBFluidPort_b outlet_cold(redeclare package Medium = PBMedia) "Recuperator outlet";  

  replaceable package PBMedia = Steps.Media.SCO2; 
  
  parameter String name_material = "inconel 750";    
  
  MaterialConductivity mc(name_material = name_material);
  // On design parameters
  // Geometry parameters  
  inner parameter Modelica.SIunits.Length pitch = 10 "pitch length of channel";
  
  inner parameter Modelica.SIunits.Angle phi = 1.0 "angle of zigzag path, unit rad";

  inner parameter Modelica.SIunits.Length length_cell = 1e-3 "length of a cell";
  
  parameter Modelica.SIunits.Diameter d_c = 0.0 "Diameter of semi-circular channel";
  
  inner parameter Integer N_ch = 10 "Number of Channels in PCHE";
  
  parameter Integer N_seg = 1 "Number of segments in a tube";

  // Geometry determined correlation coefficients - a, b, c d
  inner KimCorrelations kim_cor(phi = phi, pitch = pitch, d_h = d_h); 
  
  constant Integer MAX_SEG_NUM = Steps.Utilities.CoolProp.MAX_SEG_NUM;

  PCHEGeoParam geo(pitch = pitch, phi = phi, length_cell = length_cell, d_c = d_c,
  N_ch = N_ch, N_seg = N_seg);
  
  SimParam sim_param(err=1e-2, delta_T_init = 5, N_iter = 100, step_rel=0.3);
  
  // only (p, h, mdot) of the hot/cold inlet matters
  BoundaryCondition bc(
  st_hot_in(T = 0, p = inlet_hot.p, h = inlet_hot.h_outflow, mdot = inlet_hot.m_flow), 
  st_cold_in(T = 0, p = inlet_cold.p, h = inlet_cold.h_outflow, mdot = inlet_cold.m_flow), 
  st_hot_out(T = 0, p = 0, h = 0, mdot = inlet_hot.m_flow),  
  st_cold_out(T = 0, p = 0, h = 0, mdot = inlet_cold.m_flow));   
  
  Real h_hot[MAX_SEG_NUM];
  Real h_cold[MAX_SEG_NUM];
  Real p_hot[MAX_SEG_NUM];
  Real p_cold[MAX_SEG_NUM];    
  
  KimCorrCoe cor(a=kim_cor.a, b=kim_cor.b, c=kim_cor.c, d=kim_cor.d);  
  
  // following variables are for debug purpose only
  // d_c determined variables, d_h, A_c, peri_c
  inner Modelica.SIunits.Diameter d_h = 4 * A_c / peri_c "Hydraulic Diameter";
  
  inner Modelica.SIunits.Area A_c = Modelica.Constants.pi * d_c * d_c / 8 "Area of semi-circular tube";   
  
  inner Modelica.SIunits.Length peri_c = d_c * Modelica.Constants.pi / 2 + d_c "perimeter of semi-circular";  
  
  inner Modelica.SIunits.Length t_wall = (2 - Modelica.Constants.pi  / 4) * (d_c / 2) "thickness of wall between two neighboring hot and cold";   
  
  inner Modelica.SIunits.Area A_stack = peri_c * length_cell * N_ch "surface area of all cells in a stack";
  
  inner Modelica.SIunits.Area A_flow = N_ch * A_c "Flow area for all channels";  
 
  Modelica.SIunits.Length length_ch = length_cell * N_seg "length of one pipe in HeatExchanger unit m";   

equation   
  
  (h_hot, h_cold, p_hot, p_cold) = CP.PCHE_OFFD_Simulation("CO2", "CO2", geo, cor, sim_param, bc, geo.N_seg);

  inlet_hot.m_flow + outlet_hot.m_flow = 0;

  inlet_cold.m_flow + outlet_cold.m_flow = 0;
  
  inlet_hot.h_outflow = inStream(inlet_hot.h_outflow);
  
  inlet_cold.h_outflow = inStream(inlet_cold.h_outflow);  


// algorithm
  
  outlet_hot.h_outflow = h_hot[N_seg];
  outlet_hot.p = p_hot[N_seg];
  
  outlet_cold.h_outflow = h_cold[0];
  outlet_cold.p = p_cold[0];
  
end PCHECImpl;
