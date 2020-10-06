within Steps.Cycle;

model SimpleCycle
  "Simplest Cycle with pump - oilheater - valve - Water colder - pump"   
  constant Real CONV_DEG_TO_K = 273.15 "Temperature conversion constant degC -> K"; 
  
  parameter Modelica.SIunits.Pressure P_ATM = 101325; // Pa
  String FLUID_NAME = "CO2";
  
  parameter Modelica.SIunits.MassFlowRate M_CO2 = 0.25;  //# Co2 flow rate, kg/s    
  
  parameter Modelica.SIunits.Pressure P_PUMP_I = 7.5 * 1e6 ; //# Pa
  parameter Modelica.SIunits.AbsolutePressure P_PUMP_E = 16.0 * 1e6;  //# Pump exit pressure, Pa  
  parameter Modelica.SIunits.Temperature T_PUMP_I = 10 + CONV_DEG_TO_K; // K, "Pump inlet temperature - trial value" ;      
   
  
  Steps.Components.Regulator regulator(   
    p_init = P_PUMP_I,
    T_init = T_PUMP_I,
    m_flow_init = M_CO2
  );
  
  Steps.Components.Pump pump(
    p_outlet = P_PUMP_E, 
    eta = 0.80
  );   
  
  parameter Real M_OIL = 1.2;
  
  Steps.Components.OilHeater oil_heater(
    T_hot_in = 250 + CONV_DEG_TO_K, 
    p_hot_in = P_ATM, 
    m_dot_hot = M_OIL,
    cp_hot = 3.0 * 1e3,
    eta = 0.52
  ); 
  
  Steps.Components.Valve valve(
    p_outlet = P_PUMP_I
  );
  
  parameter Real M_WATER = 2; //kg/s
  
  Steps.Components.WaterCooler water_colder(    
    T_cold_in = 5.0 + CONV_DEG_TO_K,
    pinch = 10  
  );
  
equation

  connect(regulator.outlet, pump.inlet);
  connect(pump.outlet, oil_heater.inlet);
  connect(oil_heater.outlet, valve.inlet);
  connect(valve.outlet, water_colder.inlet);
  connect(water_colder.outlet, regulator.inlet);
 
end SimpleCycle;
