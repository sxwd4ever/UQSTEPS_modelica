within Steps;

model RCBCCycle
  "Brayton Cycle with Recuperator + Recompression"  
  
  // Adjustable parameters  
  import Util = Utilities.Util;
  parameter Modelica.SIunits.Pressure P_ATM = 101325; // Pa
  parameter Modelica.SIunits.Temperature T_AMB = Utilities.Util.ToK(15) "Ambinent temperature";
  //String FLUID_NAME = "CO2";  
  
  parameter Modelica.SIunits.Time stop_time = 1.0 "time length of the experiment";
  
  parameter Real M_CO2 = 51.91;  //# Co2 flow rate, kg/s    
  
  parameter Modelica.SIunits.Pressure P_PUMP_I = 8.0 * 1e6 ; //# Pa
  parameter Modelica.SIunits.Pressure P_PUMP_E = 20.0 * 1e6;  //# Pump exit pressure, Pa 
  
  Modelica.SIunits.Temperature T_PUMP_I = T_AMB + DT_COOLER; // K, "Pump inlet temperature - trial value" ; 
   
  // efficiency of main compressor, bypass_compressor and turbine
  parameter Real eta_main_compressor = 0.89;
  
  parameter Real eta_bypass_compressor = 0.89;
  
  parameter Real eta_turbine = 0.9;
  
  // efeectiveness for two recuperators
  // these two parameters will effect the difference between crec_in.T and hot_in.T of recuperator  
  parameter Real eta_recuperator_high = 0.99;
  
  parameter Real eta_recuperator_low = 0.99;
  
  // mass split ratio of splitter
  parameter Real splitter_split_ratio = 0.4; 
  
  	// to-be swiped parameters - now only split_ration contained 
  //parameter String data_file = Modelica.Utilities.Files.loadResource("D:/sxwd/Projects/UQ/2019.10.01 Thermal Cycle Simulation/dev/steps_modelica/Steps/data.dat");
  
  //DataTable data(file = data_file);

  Steps.Components.Regulator regulator(
    //fluid(name = FLUID_NAME, cp = 1 * 1e3),
    p_init = P_PUMP_E,   
    m_flow_init = M_CO2
  );
  
  parameter Modelica.SIunits.TemperatureDifference DT_COOLER = 18.0;
  
  Steps.Components.FanCooler fan_cooler(
    //fluid(name = FLUID_NAME, cp = 1),
    T_amb = T_AMB,
    delta_T = DT_COOLER
  );
  
  // use pump as compressor
  Steps.Components.Pump pump(
    //fluid(name = FLUID_NAME, cp = PropsSI("C", "T", pump.inlet.T, "P", pump.inlet.p, pump.fluid.name)),
    p_outlet = P_PUMP_E, 
    eta = eta_main_compressor
  );   
  
  Steps.Components.Pump recom_pump(
    //fluid(name = FLUID_NAME, cp = PropsSI("C", "T", recom_pump.inlet.T, "P", recom_pump.inlet.p, recom_pump.fluid.name)),
     p_outlet = P_PUMP_E, 
    eta = eta_bypass_compressor
  );
  
  //high temperature recuperator
  
  Steps.Components.Recuperator recup_high(
    //fluid(name = FLUID_NAME, cp = PropsSI("C", "T", recup_high.inlet.T, "P", recup_high.inlet.p, recup_high.fluid.name)),
    //crec_fluid(name = FLUID_NAME, cp = PropsSI("C", "T", recup_high.crec_in.T, "P", recup_high.crec_in.p, recup_high.crec_fluid.name)),     
    eta = eta_recuperator_high  
  );
  
  Steps.Components.Recuperator recup_low(
    //fluid(name = FLUID_NAME, cp = PropsSI("C", "T", recup_low.inlet.T, "P", recup_low.inlet.p, recup_low.fluid.name)),
    //crec_fluid(name = FLUID_NAME, cp = PropsSI("C", "T", recup_low.crec_in.T, "P", recup_low.crec_in.p, recup_low.crec_fluid.name)),     
    eta = eta_recuperator_low   
  );
  
  /*
  Steps.Recuperator_2 recup_high(
    fluid(name = FLUID_NAME, cp = PropsSI("C", "T", recup_high.inlet.T, "P", recup_high.inlet.p, recup_high.fluid.name)),
    crec_fluid(name = FLUID_NAME, cp = PropsSI("C", "T", recup_high.crec_in.T, "P", recup_high.crec_in.p, recup_high.crec_fluid.name)),     
    pinch = 5   
  );
   
  Steps.Recuperator_2 recup_low(
    fluid(name = FLUID_NAME, cp = PropsSI("C", "T", recup_low.inlet.T, "P", recup_low.inlet.p, recup_low.fluid.name)),
    crec_fluid(name = FLUID_NAME, cp = PropsSI("C", "T", recup_low.crec_in.T, "P", recup_low.crec_in.p, recup_low.crec_fluid.name)),     
    pinch = 5    
  );
  */
  
  parameter Real dT_pcm = 10.0;

  Steps.Components.TemperatureOutput temp_out(
    T_start = Utilities.Util.ToK(700),
    T_stop = Utilities.Util.ToK(700),
    dT_step = 50,
    t_sim_duration = stop_time
  );
  
  Steps.Components.PCMHeater pcm_heater(
    //fluid(name = FLUID_NAME, cp = 1 * 1e3),    
    delta_T = dT_pcm
  );  
  
  parameter Modelica.SIunits.Pressure P_TURBINE_E = P_PUMP_I;
  
  Steps.Components.Turbine turbine(
    //fluid(name = FLUID_NAME, cp = PropsSI("C", "T", turbine.inlet.T, "P", turbine.inlet.p, turbine.fluid.name)),
    p_out = P_TURBINE_E,
    eta = eta_turbine
  ); 
  
  Steps.Components.Splitter splitter(
    //fluid(name = FLUID_NAME, cp = 1 * 1e3),
    //split_ratio = data.split_ratio
    split_ratio = splitter_split_ratio
  );
  
  Steps.Components.Merger merger(
    //fluid(name = FLUID_NAME, cp = PropsSI("C", "T", merger.inlet.T, "P", merger.inlet.p, merger.fluid.name)),
    T_init = T_AMB,
    p_init = P_ATM
    //m_init = 0
  );
  
  //total efficiency
  //Real eta_total;
  
  //Real W_comp, W_turbine, Q_H;    
  
equation
  
  connect(regulator.outlet, turbine.inlet);
  
  connect(turbine.outlet, recup_high.inlet);
  
  connect(recup_high.outlet, recup_low.inlet); 
 
  connect(recup_low.outlet, splitter.inlet);
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
  /*
  W_comp := pump.inlet.m_flow * (pump.h_ea - pump.h_i) + recom_pump.inlet.m_flow * (recom_pump.h_ea - recom_pump.h_i);
  
  W_turbine := turbine.inlet.m_flow * (turbine.h_i - turbine.h_ea);
  Q_H := pcm_heater.inlet.m_flow * (pcm_heater.h_e - pcm_heater.h_i);
  
  eta_total := if initial() then 0 else (W_turbine - W_comp) / Q_H * 100;
  */
end RCBCCycle;
