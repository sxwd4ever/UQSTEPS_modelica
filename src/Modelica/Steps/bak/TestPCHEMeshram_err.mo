within Steps.Test;

model TestPCHEMeshram
  "PCHE Test against Meshram [2016]"      
  
  import Modelica.SIunits.Conversions.{from_bar, from_degC};   
  import Steps.Utilities.CoolProp.PropsSI;
  
  replaceable package PBMedia = Steps.Media.SCO2;
  
  // **** cold_inlet **** 
  parameter Modelica.SIunits.Temp_K T_cold_in = from_degC(151.4) "T of cold inlet - High Temp range for zigzag in 4.1.2 Meshram [2016]";  
  parameter Modelica.SIunits.AbsolutePressure p_cold_in = from_bar(200) "p of cold inlet - High Temp range for zigzag in 4.1.2 Meshram [2016]"; // Find THIS.
  parameter Modelica.SIunits.SpecificEnthalpy h_cold_in = PropsSI("H", "T", T_cold_in, "P", p_cold_in, PBMedia.mediumName);
  
  // **** cold_outlet ****
  parameter Modelica.SIunits.Temp_K T_cold_out = from_degC(533.4) "T of cold outlet - High Temp range for zigzag in 4.1.2 Meshram [2016]";  
  parameter Modelica.SIunits.AbsolutePressure p_cold_out = p_cold_in "p of cold outlet - High Temp range for zigzag in 4.1.2 Meshram [2016]";  
  parameter Modelica.SIunits.SpecificEnthalpy h_cold_out = PropsSI("H", "T", T_cold_out, "P", p_cold_out, PBMedia.mediumName);

  // **** hot_inlet ****
  parameter Modelica.SIunits.Temp_K T_hot_in = from_degC(578.2) "T of hot inlet - High Temp range for zigzag in 4.1.2 Meshram [2016]";    
  parameter Modelica.SIunits.AbsolutePressure p_hot_in = from_bar(80) "p of hot inlet - High Temp range for zigzag in 4.1.2 Meshram [2016]"; // Find THIS.  
  parameter Modelica.SIunits.SpecificEnthalpy h_hot_in = PropsSI("H", "T", T_hot_in, "P", p_hot_in, PBMedia.mediumName);
  
  
  // **** hot_outlet ****
  parameter Modelica.SIunits.Temp_K T_hot_out = from_degC(156.4) "T of hot outlet - High Temp range for zigzag in 4.1.2 Meshram [2016]"; 
  parameter Modelica.SIunits.AbsolutePressure p_hot_out = p_hot_in "p of hot outlet - High Temp range for zigzag in 4.1.2 Meshram [2016]"; 
  parameter Modelica.SIunits.SpecificEnthalpy h_hot_out = PropsSI("H", "T", T_hot_out, "P", p_hot_out, PBMedia.mediumName); 
  
  //parameter Modelica.SIunits.Angle phi_array[4] = {from_deg(0), from_deg(0), from_deg((180 - 108) /2), from_deg((180 - 108) /2)}; // 4.1.2 of Meshram [2016]
  
  parameter Modelica.SIunits.Angle phi = from_deg((180 - 108) /2) "phi of zigzag - High Temp range for zigzag in 4.1.2 Meshram [2016]"; // agree with Hal's steps instead of using 4.1.2 of Meshram [2016]

  parameter Integer N_ch = integer(94e3);  
  
  parameter Modelica.SIunits.ReynoldsNumber Re_hot_start = 14.5e3 "Hot stream's start value of Reynolds Number, used to increase convergence";
  
  parameter Modelica.SIunits.ReynoldsNumber Re_cold_start = 14.5e3 "Cold stream's start value of Reynolds Number, used to increase convergence";  
  
  parameter Modelica.SIunits.MassFlowRate mdot_hot = 51.91 "mass flow rate for hot stream";
  
  parameter Modelica.SIunits.MassFlowRate mdot_cold = 51.91 "mass flow rate for cold stream";

  parameter Modelica.SIunits.Length d_c = 2e-3 "diameter of the channel";

  parameter Modelica.SIunits.Length pitch = 12e-3 "pitch length of zigzag channel";

  parameter Modelica.SIunits.Length length = 2860 "length of the discretized cell in a channel";

  parameter Integer N_seg = 50 "number of cells/segment for the discretization of a channel";
  
  parameter Boolean SourceFixed_hot = true;
  
  parameter Boolean SourceFixed_cold = true;
  
  Components.Source source_hot(
    p_outlet = p_hot_in,
    T_outlet = T_hot_in,
    mdot_init = mdot_hot,
    fix_state = SourceFixed_hot,
    outlet.p(start = p_hot_in),
    outlet.h_outflow(start = h_hot_in)
  );

  Components.Source source_cold(
    p_outlet = p_cold_in,
    T_outlet = T_cold_in,
    mdot_init = mdot_cold,
    fix_state = SourceFixed_cold, // True if set its state as boundary condition
    outlet.p(start = p_cold_in),
    outlet.h_outflow(start = h_cold_in)
  );

  Components.Sink sink_hot(
    p_inlet = p_hot_out,
    T_inlet = T_hot_out,
    mdot_init = mdot_hot,
    fix_state = not SourceFixed_hot,
    inlet.p(start = p_hot_out),
    inlet.h_outflow(start = h_hot_out)
  );

  Components.Sink sink_cold(
    p_inlet = p_cold_out,
    T_inlet = T_cold_out,
    mdot_init = mdot_cold,
    fix_state = not SourceFixed_cold,
    inlet.p(start = p_cold_out),
    inlet.h_outflow(start = h_cold_out)
  );
 
  Components.PCHeatExchanger pche(
    phi = phi, 
    d_c = d_c,
    pitch = pitch,     
    N_ch = N_ch,    
    Re_cold_start = Re_cold_start,
    Re_hot_start = Re_hot_start,
    
    h_start_hot = h_hot_out,
    h_start_cold = h_cold_out,
    p_start_hot = p_hot_in,
    p_start_cold = p_cold_out,
    
    
    inlet_hot.p(start = p_hot_in),
    inlet_hot.h_outflow(start = h_hot_in),    
    
    inlet_cold.p(start = p_cold_in),
    inlet_cold.h_outflow(start = h_cold_in),    
    
    //outlet_cold.p(start = p_cold_out),
    //outlet_hot.p(start = p_hot_out),
    
    ByInlet_hot = SourceFixed_hot,
    ByInlet_cold = SourceFixed_cold,
    N_seg = N_seg,
    length = length,
    mdot_start_hot = mdot_hot,
    mdot_start_cold = mdot_cold 
  );
  
equation
  
  connect(source_hot.outlet, pche.inlet_hot);
  connect(pche.outlet_hot, sink_hot.inlet);
  connect(source_cold.outlet, pche.inlet_cold);
  connect(pche.outlet_cold, sink_cold.inlet);  
  
annotation(
    experiment(StartTime = 0, StopTime = 1, Interval = 1, Tolerance = 1e-6),
    __OpenModelica_commandLineOptions = "--matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection -d=initialization,NLSanalyticJacobian,aliasConflicts");  
end TestPCHEMeshram;
