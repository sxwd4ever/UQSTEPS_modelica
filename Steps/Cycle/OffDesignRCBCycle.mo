within Steps.Cycle;

model OffDesignRCBCycle
  "OffDesign Brayton Cycle with Recuperator + Recompression"  
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
  
  parameter Modelica.SIunits.Angle phi = from_deg((180 - 108) /2) "phi of zigzag - High Temp range for zigzag in 4.1.2 Meshram [2016]"; // agree with Hal's steps instead of using 4.1.2 of Meshram [2016]

  parameter Integer N_ch = integer(80e3);  
  
  parameter Modelica.SIunits.Length d_c = 2e-3 "diameter of the channel";

  parameter Modelica.SIunits.Length pitch = 12e-3 "pitch length of zigzag channel";

  parameter Modelica.SIunits.Length length_cell = 12e-3 "length of the discretized cell in a channel";

  parameter Integer N_seg = 10 "number of cells/segment for the discretization of a channel";  
  
  parameter Modelica.SIunits.ReynoldsNumber Re_hot_start = 2e3 "Hot stream's start value of Reynolds Number, used to increase convergence";
  
  parameter Modelica.SIunits.ReynoldsNumber Re_cold_start = 2e3 "Cold stream's start value of Reynolds Number, used to increase convergence";    

  parameter Modelica.SIunits.TemperatureDifference DT_COOLER = 18.0;

  Steps.Components.Regulator regulator(
    T_init = Modelica.SIunits.Conversions.from_degC(10),
    p_init = P_PUMP_E,   
    m_flow_init = M_CO2
  );
  
  Steps.Components.FanCooler fan_cooler(
    T_amb = T_AMB,
    delta_T = DT_COOLER
  );
  
  // use pump as compressor
  Steps.Components.Pump pump(
    p_outlet = P_PUMP_E, 
    eta = eta_main_compressor,
    medium_in.T.start = 500, 
    medium_in.T.nominal = 500
  );   
  
  Steps.Components.Pump recom_pump(
    p_outlet = P_PUMP_E, 
    eta = eta_bypass_compressor
  );
  
  
  //high temperature recuperator
  Steps.Components.Recuperator recup_high(
    eta = eta_recuperator_high      
  );
  
  //low temperature recuperator
  Steps.Components.Recuperator recup_low(
    eta = eta_recuperator_low,
    medium_hot_max.T.start = 500,
    medium_hot_max.T.nominal=500
  );
  
  /*
  Components.PCHeatExchanger recup_high(
    phi = phi, 
    d_c = d_c,
    pitch = pitch,     
    N_ch = N_ch,
    Re_cold_start = Re_cold_start,
    Re_hot_start = Re_hot_start,
    N_seg = N_seg,
    length_cell = length_cell  
  );
  
  Components.PCHeatExchanger recup_low(
    phi = phi, 
    d_c = d_c,
    pitch = pitch,     
    N_ch = N_ch,
    Re_cold_start = Re_cold_start,
    Re_hot_start = Re_hot_start,
    N_seg = N_seg,
    length_cell = length_cell  
  );
 */
  parameter Real dT_pcm = 10.0;

  Steps.Components.TemperatureOutput temp_out(
    T_start = Modelica.SIunits.Conversions.from_degC(700),
    T_stop = Modelica.SIunits.Conversions.from_degC(700),
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
  
  connect(splitter.outlet, fan_cooler.inlet);
  
  connect(fan_cooler.outlet, pump.inlet);
  
  connect(pump.outlet, recup_low.inlet_cool);  
   
  connect(recup_low.outlet_cool, merger.inlet);
  
  connect(merger.outlet, recup_high.inlet_cool);
  
  connect(recup_high.outlet_cool, pcm_heater.inlet);    
  
  connect(temp_out.y, pcm_heater.T_input);
  connect(pcm_heater.outlet, regulator.inlet);

algorithm
  eta_total := if initial() then 0 else (turbine.W_turbine - pump.W_comp - recom_pump.W_comp) / pcm_heater.Q * 100;
annotation(
    experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-6, Interval = 0.002),
    __OpenModelica_commandLineOptions = "--matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection -d=initialization,NLSanalyticJacobian");
end OffDesignRCBCycle;
