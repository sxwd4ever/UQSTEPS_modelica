within Steps;

model TestPCHX
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
    m_dot_flow = 8.3
  );

  Components.Source source_cool(
    p_outlet = 20 * 1e6,
    T_outlet = Modelica.SIunits.Conversions.from_degC(41),
    m_dot_flow = 8.3
  );

  Components.Sink sink_hot(
    p_inlet = 8.88 * 1e6,
    T_inlet = Modelica.SIunits.Conversions.from_degC(51),
    fix_sink = false
  );

  Components.Sink sink_cool(
    p_inlet = 20 * 1e6,
    T_inlet = Modelica.SIunits.Conversions.from_degC(332.9),
    fix_sink = false
  );
 
  Components.MockPCHeatExchanger pchx(
    phi = Modelica.SIunits.Conversions.from_deg(45), 
    Re_design = 5000,
    d_c = 1.51 * 1e-3,
    T_hot_in = Modelica.SIunits.Conversions.from_degC(451),
    T_cool_in = Modelica.SIunits.Conversions.from_degC(41),
    p_hot = 9 * 1e6,
    p_cool = 20 * 1e6,
    pitch = 24.6 * 1e-3
  );
  
  /*
  Components.PCHeatExchanger pchx(
    phi = Modelica.Constants.pi / 4, 
    Re_design = 5000,
    d_c = 1.51 * 1e-3,
    T_hot_in = Modelica.SIunits.Conversions.from_degC(451),
    T_cool_in = Modelica.SIunits.Conversions.from_degC(41),
    p_hot = 9 * 1e6,
    p_cool = 20 * 1e6
  );
  */
  
  /*
  Components.Recuperator pchx(
    eta = 0.99  
  );
  */
  
equation
  /*
  connect(source.outlet, pchx.inlet_hot);
  connect(pchx.outlet_hot, sink.inlet);
  connect(sink.outlet, pchx.inlet_cool);
  connect(pchx.outlet_cool, source.inlet);
  */
  
  connect(source_hot.outlet, pchx.inlet_hot);
  connect(pchx.outlet_hot, sink_hot.inlet);
  connect(source_cool.outlet, pchx.inlet_cool);
  connect(pchx.outlet_cool, sink_cool.inlet);
  
  //connect(source_hot.outlet, sink_hot.inlet);
  
end TestPCHX;
