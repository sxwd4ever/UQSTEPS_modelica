within Steps.Test;

model TestPCHECImpl
  "PCHE c-based implementation"      
  
  import Modelica.SIunits.Conversions.{from_bar, from_deg};   
  import Steps.Utilities.CoolProp.PropsSI;
  
  replaceable package PBMedia = Steps.Media.SCO2;
  
  parameter Modelica.SIunits.Temp_K T_cold_in = 500 "T of cold inlet - High Temp range for zigzag in 4.1.2 Meshram [2016]";
  
  parameter Modelica.SIunits.AbsolutePressure p_cold_in = from_bar(225) "p of cold inlet - High Temp range for zigzag in 4.1.2 Meshram [2016]"; // Find THIS.

  parameter Modelica.SIunits.Temp_K T_cold_out = 639.15 "T of cold outlet - High Temp range for zigzag in 4.1.2 Meshram [2016]";
  
  parameter Modelica.SIunits.AbsolutePressure p_cold_out = from_bar(225) "p of cold outlet - High Temp range for zigzag in 4.1.2 Meshram [2016]";  

  parameter Modelica.SIunits.Temp_K T_hot_in = 730 "T of hot inlet - High Temp range for zigzag in 4.1.2 Meshram [2016]";    
  parameter Modelica.SIunits.AbsolutePressure p_hot_in = from_bar(90) "p of hot inlet - High Temp range for zigzag in 4.1.2 Meshram [2016]"; // Find THIS.  
  
  parameter Modelica.SIunits.Temp_K T_hot_out = 576.69 "T of hot outlet - High Temp range for zigzag in 4.1.2 Meshram [2016]"; 
  parameter Modelica.SIunits.AbsolutePressure p_hot_out = from_bar(90) "p of hot outlet - High Temp range for zigzag in 4.1.2 Meshram [2016]";  
  
  //parameter Modelica.SIunits.Angle phi_array[4] = {from_deg(0), from_deg(0), from_deg((180 - 108) /2), from_deg((180 - 108) /2)}; // 4.1.2 of Meshram [2016]
  
  parameter Modelica.SIunits.Angle phi = from_deg((180 - 108) /2) "phi of zigzag - High Temp range for zigzag in 4.1.2 Meshram [2016]"; // agree with Hal's steps instead of using 4.1.2 of Meshram [2016]

  parameter Integer N_ch = integer(80e3);  
  
  parameter Modelica.SIunits.ReynoldsNumber Re_hot_start = 2e3 "Hot stream's start value of Reynolds Number, used to increase convergence";
  
  parameter Modelica.SIunits.ReynoldsNumber Re_cold_start = 2e3 "Cold stream's start value of Reynolds Number, used to increase convergence";  
  
  parameter Modelica.SIunits.MassFlowRate mdot_hot = 10 "mass flow rate for hot stream";
  
  parameter Modelica.SIunits.MassFlowRate mdot_cold = 10 "mass flow rate for cold stream";

  parameter Modelica.SIunits.Length d_c = 2e-3 "diameter of the channel";

  parameter Modelica.SIunits.Length pitch = 12e-3 "pitch length betweem zigzag channels";

  parameter Modelica.SIunits.Length length_cell = 12e-3 "length of the discretized cell in a channel";

  parameter Integer N_seg = 100 "number of cells/segment for the discretization of a channel";
  
  parameter Boolean SourceFixed_hot = true;
  
  parameter Boolean SourceFixed_cold = true;
  
  Components.Source source_hot(
    p_outlet = p_hot_in,
    T_outlet = T_hot_in,
    mdot_init = mdot_hot,
    fix_state = SourceFixed_hot
  );

  Components.Source source_cold(
    p_outlet = p_cold_in,
    T_outlet = T_cold_in,
    mdot_init = mdot_cold,
    fix_state = SourceFixed_cold // True if set its state as boundary condition
  );

  Components.Sink sink_hot(
    p_inlet = p_hot_out,
    T_inlet = T_hot_out,
    mdot_init = mdot_hot,
    fix_state = not SourceFixed_hot
  );

  Components.Sink sink_cold(
    p_inlet = p_cold_out,
    T_inlet = T_cold_out,
    mdot_init = mdot_cold,
    fix_state = not SourceFixed_cold
  );
 
  Components.PCHECImpl pche(
    phi = phi, 
    d_c = d_c,
    pitch = pitch,     
    N_ch = N_ch,    
    N_seg = N_seg,
    length_cell = length_cell,    
    
    bc_hot_in(p = p_hot_in, T = T_hot_in, mdot = mdot_hot),    
    bc_cold_in(p = p_cold_in, T = T_cold_in, mdot = mdot_cold)    
    /*
    inlet_hot.p(start = p_hot_in, nominal = p_hot_in),
    inlet_hot.h_outflow(start = h_hot_in, nominal = h_hot_in),
    inlet_hot.m_flow(start = mdot_hot, nominal = mdot_hot),
    
    inlet_cold.p(start = p_cold_in, nominal = p_cold_in),
    inlet_cold.h_outflow(start = h_cold_in, nominal = h_cold_in),
    inlet_cold.m_flow(start = mdot_cold, nominal = mdot_cold) 
    */
  );
  

  
equation
  
  connect(source_hot.outlet, pche.inlet_hot);
  connect(pche.outlet_hot, sink_hot.inlet);
  connect(source_cold.outlet, pche.inlet_cold);
  connect(pche.outlet_cold, sink_cold.inlet);  
  
end TestPCHECImpl;
