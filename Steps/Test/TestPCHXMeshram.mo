within Steps.Test;

model TestPCHXMeshram
  "PCHE Test against Meshram [2016]"      
  
  import Modelica.SIunits.Conversions.{from_bar, from_deg};
  
  type TestScenario = enumeration(straight_low_T, straight_high_T, zigzag_low_T, zigzag_high_T);
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
  
  parameter Modelica.SIunits.Angle phi_array[4] = {from_deg(0), from_deg(0), from_deg((180 - 110) /2), from_deg((180 - 110) /2)}; // agree with Hal's steps instead of using 4.1.2 of Meshram [2016]
  
  // select one test scenario in Meshram [2016]
  parameter Integer scenario = Integer(TestScenario.zigzag_high_T);
  
  parameter Modelica.SIunits.MassFlowRate mdot_hot_init = 99.4621884201685;
  
  parameter Modelica.SIunits.MassFlowRate mdot_cold_init = 100.891159281523;
  
  Components.Source source_hot(
    p_outlet = p_hot_in[scenario],
    T_outlet = T_hot_in[scenario],
    mdot_init = mdot_hot_init,
    fix_state = true
  );

  Components.Source source_cold(
    p_outlet = p_cold_in[scenario],
    T_outlet = T_cold_in[scenario],
    mdot_init = mdot_cold_init,
    fix_state = true
  );

  Components.Sink sink_hot(
    p_inlet = p_hot_out[scenario],
    T_inlet = T_hot_out[scenario],
    mdot_init = mdot_hot_init,
    fix_state = false
  );

  Components.Sink sink_cold(
    p_inlet = p_cold_out[scenario],
    T_inlet = T_cold_out[scenario],
    mdot_init = mdot_cold_init,
    fix_state = false
  );
 
  Components.PCHeatExchanger pchx(
    phi = phi_array[scenario], 
    Re_design = 2000,
    d_c = 2e-3,
    T_hot_in = T_hot_in[scenario],
    T_cool_in = T_cold_in[scenario],
    p_hot = p_hot_in[scenario],
    p_cool = p_cold_in[scenario],
    m_dot_hot = mdot_hot_init,
    m_dot_cool = mdot_cold_init,
    pitch = 12.3e-3,
    length_cell = 12e-3,
    N_seg = 10
  );
  
equation
  
  connect(source_hot.outlet, pchx.inlet_hot);
  connect(pchx.outlet_hot, sink_hot.inlet);
  connect(source_cold.outlet, pchx.inlet_cool);
  connect(pchx.outlet_cool, sink_cold.inlet);  
  
end TestPCHXMeshram;
