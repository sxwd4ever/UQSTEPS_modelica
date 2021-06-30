within Steps.Test.Components;

model TestPCHE
  "Simple cycle to test PC Heat Exchange"     
  
  /*
  Components.EndPoint source(
    p_outlet = 9 * 1e6,
    T_outlet = Modelica.SIunits.Conversions.from_degC(451),
    m_dot_flow = 8.3
  );
  
  Components.EndPoint sink(
    p_outlet = 20 * 1e6,
    T_outlet = Modelica.SIunits.Conversions.from_degC(332.9),
    m_dot_flow = 8.3
  );
  */
  
  Components.Source source_hot(
    p_outlet = 9 * 1e6,
    T_outlet = Modelica.SIunits.Conversions.from_degC(451),
    mdot_init = 8.3,
    fix_state = false
  );

  Components.Source source_cold(
    p_outlet = 20 * 1e6,
    T_outlet = Modelica.SIunits.Conversions.from_degC(41),
    mdot_init = 8.3,
    fix_state = true
  );

  Components.Sink sink_hot(
    p_inlet = 8.88 * 1e6,
    T_inlet = Modelica.SIunits.Conversions.from_degC(51),
    mdot_init = 8.3,
    fix_state = true
  );

  Components.Sink sink_cold(
    p_inlet = 20 * 1e6,
    T_inlet = Modelica.SIunits.Conversions.from_degC(332.9),
    mdot_init = 8.3,
    fix_state = false
  );
 
  Components.PCHeatExchanger pche(
    
phi = Modelica.SIunits.Conversions.from_deg(45), 
    Re_design = 5000,
    d_c = 1.51 * 1e-3,
    T_hot_in = Modelica.SIunits.Conversions.from_degC(451),
    T_cold_in = Modelica.SIunits.Conversions.from_degC(41),
    p_hot = 9 * 1e6,
    p_cold = 20 * 1e6,
    m_dot_hot = 8.3,
    m_dot_cold = 8.3,
    pitch = 24.6 * 1e-3,
    length_cell = 8e-2,
    N_seg = 20
  );
  
  /*
  Components.PCHeatExchangerV2 pche(
    phi = Modelica.SIunits.Conversions.from_deg(45), 
    Re_design = 5000,
    d_c = 1.51 * 1e-3,
    T_hot_in = Modelica.SIunits.Conversions.from_degC(451),
    T_cold_in = Modelica.SIunits.Conversions.from_degC(41),
    p_hot = 9 * 1e6,
    p_cold = 20 * 1e6,
    m_dot_hot = 8.3,
    m_dot_cold = 8.3,
    pitch = 24.6 * 1e-3,
    length_cell = 3e-3,
    N_seg = 2
  );
  */
  
  /*
  Components.PCHeatExchanger pche(
    phi = Modelica.Constants.pi / 4, 
    Re_design = 5000,
    d_c = 1.51 * 1e-3,
    T_hot_in = Modelica.SIunits.Conversions.from_degC(451),
    T_cold_in = Modelica.SIunits.Conversions.from_degC(41),
    p_hot = 9 * 1e6,
    p_cold = 20 * 1e6
  );
  */
  
  /*
  Components.Recuperator pche(
    eta = 0.99  
  );
  */
  
equation
  
  connect(source_hot.outlet, pche.inlet_hot);
  connect(pche.outlet_hot, sink_hot.inlet);
  connect(source_cold.outlet, pche.inlet_cold);
  connect(pche.outlet_cold, sink_cold.inlet);  
  
end TestPCHE;
