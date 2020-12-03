within Steps.Model;

model PBConfig_5MW
  "Preset param set for Off-design 5MW power block test"
  extends PBConfiguration(
    p_pump_in = 8e6,
    p_pump_out = 20e6,
    p_ATM = 1.01e6,
    p_heater = 20e6,
    
    mdot_main = 51.91,
    mdot_pump = 31.31,
    mdot_heater = 42.73,
    mdot_cooler = 335.2,
    
    T_amb = from_degC(15),
    T_HTR_hot_in = from_degC(578.22),
    T_HTR_cold_out = from_degC(533.5),
    
    T_HTR_hot_out = from_degC(156.45),
    T_HTR_cold_in = from_degC(151.45), 
    
    T_LTR_cold_in = from_degC(62.229),
    T_LTR_hot_out = from_degC(67.229),  
    
    T_heater_hot_out = from_degC(600),
    T_heater_hot_in = from_degC(800),
    // fixed pre-defined condition
    T_heater_cold_out = from_degC(700),
    
    T_cooler_cold_out = from_degC(30),
    T_cooler_cold_in = pbcfg_5MW.T_amb, 
    // fixed pre-defined condition
    T_cooler_hot_out = from_degC(33)
  );


end PBConfig_5MW;
