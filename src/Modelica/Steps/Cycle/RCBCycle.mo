within Steps.Cycle;

model RCBCycle
  "Brayton Cycle with Recuperator + Recompression"  
// Adjustable parameters
  import Util = Utilities.Util;
  parameter Modelica.SIunits.Pressure P_ATM = 101325; // Pa
  parameter Modelica.SIunits.Temperature T_AMB = Utilities.Util.ToK(15) "Ambinent temperature";
  
  parameter Modelica.SIunits.Time stop_time = 1.0 "time length of the experiment";
  
  parameter Real M_CO2 = 51.91;    //# Co2 flow rate, kg/s
  parameter Modelica.SIunits.Pressure P_PUMP_I = 8.0 * 1e6 ; //# Pa
  parameter Modelica.SIunits.Pressure P_PUMP_E = 20.0 * 1e6;    //# Pump exit pressure, Pa
  //Modelica.SIunits.Temperature T_PUMP_I = T_AMB + DT_COOLER; // K, "Pump inlet temperature - trial value" ;
  // efficiency of main compressor, bypass_compressor and turbine
  parameter Real eta_main_compressor = 0.89;
  
  parameter Real eta_bypass_compressor = 0.89;
  
  parameter Real eta_turbine = 0.9;
  
  // effectiveness for two recuperators
  // these two parameters will effect the difference between crec_in.T and hot_in.T of recuperator
  parameter Real eta_recuperator_high = 0.99;
  
  parameter Real eta_recuperator_low = 0.99;
  
  // mass split ratio of splitter
  parameter Real splitter_split_ratio = 0.4; 

  Steps.Components.Regulator regulator(
    T_init = Modelica.SIunits.Conversions.from_degC(10),
    p_init = P_PUMP_E,   
    m_flow_init = M_CO2
  );
  
  parameter Modelica.SIunits.TemperatureDifference DT_COOLER = 18.0;
  
  Steps.Components.FanCooler fan_colder(
    T_amb = T_AMB,
    delta_T = DT_COOLER
  );
  
  // use pump as compressor
  Steps.Components.Pump pump(
    p_outlet = P_PUMP_E, 
    eta = eta_main_compressor
  );   
  
  Steps.Components.Pump recom_pump(
    p_outlet = P_PUMP_E, 
    eta = eta_bypass_compressor
  );
  
  /*
  //high temperature recuperator
  Steps.Components.Recuperator recup_high(
    eta = eta_recuperator_high  
  );
  
  //low temperature recuperator
  Steps.Components.Recuperator recup_low(
    eta = eta_recuperator_low   
  );
  */
  
  Components.MockPCHeatExchanger recup_high(
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
  
  Components.MockPCHeatExchanger recup_low(
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
 
  parameter Real dT_pcm = 10.0;

  Steps.Components.TemperatureOutput temp_out(
    T_start = Utilities.Util.ToK(700),
    T_stop = Utilities.Util.ToK(700),
    dT_step = 50,
    t_sim_duration = stop_time
  );
  
  Steps.Components.PCMHeater pcm_heater(
    delta_T = dT_pcm
  );  
  
  parameter Modelica.SIunits.Pressure P_TURBINE_E = P_PUMP_I;
  
  Steps.Components.Turbine turbine(    
    p_out = P_TURBINE_E,
    eta = eta_turbine
  ); 
  
  Steps.Components.Splitter splitter(
    split_ratio = splitter_split_ratio
  );
  
  Steps.Components.Merger merger(
    T_init = T_AMB,
    p_init = P_ATM
  );
  
  //total efficiency
  Real eta_total;

equation
  
  connect(regulator.outlet, turbine.inlet);
  
  connect(turbine.outlet, recup_high.inlet_hot);
  
  connect(recup_high.outlet_hot, recup_low.inlet_hot); 
 
  connect(recup_low.outlet_hot, splitter.inlet);
  //recompression loop
  connect(splitter.outlet_split, recom_pump.inlet);
  connect(recom_pump.outlet, merger.inlet_merge);  
  
  connect(splitter.outlet, fan_colder.inlet);
  
  connect(fan_colder.outlet, pump.inlet);
  
  connect(pump.outlet, recup_low.inlet_cold);  
   
  connect(recup_low.outlet_cold, merger.inlet);
  
  connect(merger.outlet, recup_high.inlet_cold);
  
  connect(recup_high.outlet_cold, pcm_heater.inlet);    
  
  connect(temp_out.y, pcm_heater.T_input);
  connect(pcm_heater.outlet, regulator.inlet);

algorithm
  eta_total := if initial() then 0 else (turbine.W_turbine - pump.W_comp - recom_pump.W_comp) / pcm_heater.Q * 100;
annotation(
    experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-6, Interval = 0.002),
    __OpenModelica_commandLineOptions = "--matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection -d=initialization,NLSanalyticJacobian");
end RCBCycle;
