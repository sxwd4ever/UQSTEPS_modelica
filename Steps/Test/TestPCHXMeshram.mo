within Steps.Test;

model TestPCHXMeshram
  "PCHE Test against Meshram [2016]"      
  
  import Modelica.SIunits.Conversions.{from_bar, from_deg};   
  /*
  type TestScenario = enumeration(straight_low_T, straight_high_T, zigzag_low_T, zigzag_high_T);
  // selection index of test scenario in Meshram [2016]
  parameter Integer scenario = Integer(TestScenario.zigzag_high_T);
  
  // configuration of different test scenario in meshram [2016] - table 3
  // arranged in same order as above enumeration values
  parameter Modelica.SIunits.Temp_K T_cold_in[4] = {400, 500, 400, 500};
  parameter Modelica.SIunits.AbsolutePressure p_cold_in[4] = {from_bar(225), from_bar(225), from_bar(225), from_bar(225)}; // Find THIS.

  parameter Modelica.SIunits.Temp_K T_cold_out[4] = {498.45, 615.48, 522.23, 639.15}; // 4.1.2 of Meshram [2016]
  parameter Modelica.SIunits.AbsolutePressure p_cold_out[4] = {from_bar(225), from_bar(225), from_bar(225), from_bar(225)};  

  parameter Modelica.SIunits.Temp_K T_hot_in[4] = {630, 730, 630, 730};    
  parameter Modelica.SIunits.AbsolutePressure p_hot_in[4] = {from_bar(90), from_bar(90), from_bar(90), from_bar(90)}; // Find THIS.  
  
  parameter Modelica.SIunits.Temp_K T_hot_out[4] = {494.37, 601.83, 466.69, 576.69}; // 4.1.2 of Meshram [2016]
  parameter Modelica.SIunits.AbsolutePressure p_hot_out[4] = {from_bar(90),from_bar(90),from_bar(90),from_bar(90)};  
  
  //parameter Modelica.SIunits.Angle phi_array[4] = {from_deg(0), from_deg(0), from_deg((180 - 108) /2), from_deg((180 - 108) /2)}; // 4.1.2 of Meshram [2016]
  
  parameter Modelica.SIunits.Angle phi_array[4] = {from_deg(0), from_deg(0), from_deg((180 - 108) /2), from_deg((180 - 108) /2)}; // agree with Hal's steps instead of using 4.1.2 of Meshram [2016]
  */

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

  parameter Modelica.SIunits.Length pitch = 12e-3 "pitch length of zigzag channel";

  parameter Modelica.SIunits.Length length_cell = 12e-3 "length of the discretized cell in a channel";

  parameter Integer N_seg = 10 "number of cells/segment for the discretization of a channel";
  
  Components.Source source_hot(
    p_outlet = p_hot_in,
    T_outlet = T_hot_in,
    mdot_init = mdot_hot,
    fix_state = true
  );

  Components.Source source_cold(
    p_outlet = p_cold_in,
    T_outlet = T_cold_in,
    mdot_init = mdot_cold,
    fix_state = true
  );

  Components.Sink sink_hot(
    p_inlet = p_hot_out,
    T_inlet = T_hot_out,
    mdot_init = mdot_hot,
    fix_state = false
  );

  Components.Sink sink_cold(
    p_inlet = p_cold_out,
    T_inlet = T_cold_out,
    mdot_init = mdot_cold,
    fix_state = false
  );
 
  Components.PCHeatExchanger pchx(
    phi = phi, 
    d_c = d_c,
    pitch = pitch,     
    N_ch = N_ch,
    Re_cold_start = Re_cold_start,
    Re_hot_start = Re_hot_start,
    N_seg = N_seg,
    length_cell = length_cell   
  );
  
equation
  
  connect(source_hot.outlet, pchx.inlet_hot);
  connect(pchx.outlet_hot, sink_hot.inlet);
  connect(source_cold.outlet, pchx.inlet_cool);
  connect(pchx.outlet_cool, sink_cold.inlet);  
  
end TestPCHXMeshram;
