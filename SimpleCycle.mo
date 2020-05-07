within Steps;

model SimpleCycle
  "Simplest Cycle with pump - oilheater - valve - Water cooler - pump"   
  constant Real CONV_DEG_TO_K = 273.15 "Temperature conversion constant degC -> K"; 
  
  parameter Modelica.SIunits.Pressure P_ATM = 101325; // Pa
  String FLUID_NAME = "CO2";
  
  parameter Modelica.SIunits.MassFlowRate M_CO2 = 0.25;  //# Co2 flow rate, kg/s    
  
  parameter Modelica.SIunits.Pressure P_PUMP_I = 7.5 * 1e6 ; //# Pa
  parameter Modelica.SIunits.AbsolutePressure P_PUMP_E = 16.0 * 1e6;  //# Pump exit pressure, Pa  
  parameter Modelica.SIunits.Temperature T_PUMP_I = 10 + CONV_DEG_TO_K; // K, "Pump inlet temperature - trial value" ;   
    
  //Source  
  //Source source(outlet(T = T_PUMP_I, p = P_PUMP_I, m_flow = M_CO2));
  //Sink sink(T = T_PUMP_I, p = P_PUMP_E);
  
  //Pump
  //parameter Real ETA_PUMP = 0.80;      
  
  Steps.Components.Regulator regulator(   
    p_init = P_PUMP_I,
    T_init = T_PUMP_I,
    m_flow_init = M_CO2
  );
  
  Steps.Components.Pump pump(
    //inlet(T = T_I, p = P_I, m_flow = M_CO2),   
    //fluid(name = FLUID_NAME, cp = PropsSI("C", "T", pump.inlet.T, "P", pump.inlet.p, pump.fluid.name)),
    p_outlet = P_PUMP_E, 
    eta = 0.80
  );   
  
  parameter Real M_OIL = 1.2;
  
  Steps.Components.OilHeater oil_heater(
    T_hot_in = 250 + CONV_DEG_TO_K, 
    p_hot_in = P_ATM, 
    m_dot_hot = M_OIL,
    //hot_out(m_flow = M_OIL),
    //hot_fluid(name = "Oil", cp = 3.0 * 1e3), // kJ/kg-K -> J/kg-K
    cp_hot = 3.0 * 1e3,
    //fluid(name = FLUID_NAME, cp = PropsSI("C", "T", oil_heater.inlet.T, "P", oil_heater.inlet.p, FLUID_NAME)),
    eta = 0.52
  ); 
  
  Steps.Components.Valve valve(
    //fluid(name = FLUID_NAME, cp = PropsSI("C", "T", valve.inlet.T, "P", valve.inlet.p, valve.fluid.name)),
    //outlet(p = P_PUMP_I)
    p_outlet = P_PUMP_I
  );
  
  parameter Real M_WATER = 2; //kg/s
  
  Steps.Components.WaterCooler water_cooler(
    //fluid(name = FLUID_NAME, cp = PropsSI("C", "T", water_cooler.inlet.T, "P", water_cooler.inlet.p, water_cooler.fluid.name)),
    T_cool_in = 5.0 + CONV_DEG_TO_K,
    //cool_in(T = 5.0 + CONV_DEG_TO_K, p = P_ATM, m_flow = M_WATER),
    //cool_out(m_flow = M_WATER),
    //cool_fluid(name = "WATER", cp = 4.186 * 1e3), // kJ/kg-K -> J/kg-K
    pinch = 10  
  );
  
equation

  connect(regulator.outlet, pump.inlet);
  connect(pump.outlet, oil_heater.inlet);
  connect(oil_heater.outlet, valve.inlet);
  connect(valve.outlet, water_cooler.inlet);
  connect(water_cooler.outlet, regulator.inlet);
 
end SimpleCycle;
