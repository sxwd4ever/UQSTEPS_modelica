within Steps;

model TestPCHX
  "Simple cycle to test PC Heat Exchange"     
  
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
  
  Components.PCHeatExchanger pchx(
    phi = Modelica.Constants.pi / 4, 
    Re_design = 5000,
    d_c = 0.8493,
    T_hot_in = Modelica.SIunits.Conversions.from_degC(451),
    T_cool_in = Modelica.SIunits.Conversions.from_degC(41),
    p_hot = 9 * 1e6,
    p_cool = 20 * 1e6
  );  
  
equation
  
  connect(source.outlet, pchx.inlet_hot);
  connect(pchx.outlet_hot, sink.inlet);
  connect(sink.outlet, pchx.inlet_cool);
  connect(pchx.outlet_cool, source.inlet);
 
end TestPCHX;
