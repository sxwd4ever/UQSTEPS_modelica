within Steps.Cycle;

model OffDesignRCBCycle_v2
  "OffDesign Brayton Cycle with Recuperator + Recompression"  
// Adjustable parameters
  import Util = Utilities.Util;
  import Modelica.SIunits.Conversions.{from_degC, from_deg};
  import Steps.Utilities.CoolProp.PropsSI;
  import Modelica.SIunits.{Temperature, Pressure, SpecificEnthalpy};
  
  replaceable package PBMedia = Steps.Media.SCO2;
  
  parameter Modelica.SIunits.Time stop_time = 1.0 "time length of the experiment";
  
  parameter Modelica.SIunits.Pressure P_ATM = 101325; // Pa
  
  //parameter Modelica.SIunits.Temperature T_AMB = from_degC(33) "Ambinent temperature";
  parameter Real M_CO2 = 10;    //# Co2 flow rate, kg/s  
  
  // Start values for all components
  // thermo states for outlet of each component are defined
  parameter Modelica.SIunits.Pressure P_COOLER_E = 8.0e6;
  parameter Modelica.SIunits.Temperature T_COOLER_E = from_degC(33);
  parameter Modelica.SIunits.SpecificEnthalpy H_COOLER_E = PropsSI("H", "T", T_COOLER_E, "P", P_COOLER_E, PBMedia.mediumName);   
  
  parameter Modelica.SIunits.Pressure P_PUMP_I = P_COOLER_E ; //# Pa
  parameter Modelica.SIunits.Temperature T_PUMP_I = T_COOLER_E; //#K
  parameter Modelica.SIunits.SpecificEnthalpy H_PUMP_I = PropsSI("H", "T", T_PUMP_I, "P", P_PUMP_I, PBMedia.mediumName);   
  
  parameter Modelica.SIunits.Pressure P_PUMP_E = 20.0e6;    //# Pump exit pressure, Pa
  parameter Modelica.SIunits.Temperature T_PUMP_E = from_degC(63); //#K
  parameter Modelica.SIunits.SpecificEnthalpy H_PUMP_E = PropsSI("H", "T", T_PUMP_E, "P", P_PUMP_E, PBMedia.mediumName);   
  
  parameter Modelica.SIunits.Pressure P_HEATER_E = P_PUMP_E;
  parameter Modelica.SIunits.Temperature T_HEATER_E = from_degC(700);  
  parameter Modelica.SIunits.SpecificEnthalpy H_HEATER_E = PropsSI("H", "T", T_HEATER_E, "P", P_HEATER_E, PBMedia.mediumName);  
  
  parameter Modelica.SIunits.Pressure P_TURBINE_E = P_PUMP_I;  
  parameter Modelica.SIunits.Temperature T_TURBINE_E = from_degC(578) "Ambinent temperature";
  parameter Modelica.SIunits.SpecificEnthalpy H_TURBINE_E = PropsSI("H", "T", T_TURBINE_E, "P", P_TURBINE_E, PBMedia.mediumName); 
  
  parameter Modelica.SIunits.Pressure P_HTR_H_E = P_TURBINE_E;
  parameter Modelica.SIunits.Temperature T_HTR_H_E = from_degC(157); 
  parameter Modelica.SIunits.SpecificEnthalpy H_HTR_H_E = PropsSI("H", "T", T_HTR_H_E, "P", P_HTR_H_E, PBMedia.mediumName);  
  
  parameter Modelica.SIunits.Pressure P_LTR_C_E = P_PUMP_E;
  parameter Modelica.SIunits.Temperature T_LTR_C_E = from_degC(151); 
  parameter Modelica.SIunits.SpecificEnthalpy H_LTR_C_E = PropsSI("H", "T", T_LTR_C_E, "P", P_LTR_C_E, PBMedia.mediumName);
  
  parameter Modelica.SIunits.Pressure P_LTR_H_E = P_TURBINE_E;
  parameter Modelica.SIunits.Temperature T_LTR_H_E = from_degC(68);
  parameter Modelica.SIunits.SpecificEnthalpy H_LTR_H_E = PropsSI("H", "T", T_LTR_H_E, "P", P_LTR_H_E, PBMedia.mediumName);
  
  parameter Real dT_pcm = 10.0;  
  
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
  parameter Real splitter_split_ratio = 0.3; 
  
  parameter Modelica.SIunits.Angle phi = from_deg((180 - 108) /2) "phi of zigzag - High Temp range for zigzag in 4.1.2 Meshram [2016]"; // agree with Hal's steps instead of using 4.1.2 of Meshram [2016]

  parameter Integer N_ch = integer(80e3);  
  
  parameter Integer N_seg = 100 "number of cells/segment for the discretization of a channel";
  
  parameter Modelica.SIunits.Length d_c = 2e-3 "diameter of the channel";

  parameter Modelica.SIunits.Length pitch = 12e-3 "pitch length of zigzag channel";

  parameter Modelica.SIunits.Length length_cell = 12e-3 "length of the discretized cell in a channel";
 
  parameter Modelica.SIunits.ReynoldsNumber Re_hot_start = 14.5e3 "Hot stream's start value of Reynolds Number, used to increase convergence";
  
  parameter Modelica.SIunits.ReynoldsNumber Re_cold_start = 14.5e3 "Cold stream's start value of Reynolds Number, used to increase convergence";    
  
  /*
  Steps.Components.Regulator regulator_HTR_hot(
    T_init = from_degC(578),
    p_init = P_TURBINE_E,   
    m_flow_init = M_CO2
  );
  
  Steps.Components.Regulator regulator_HTR_cold(
    T_init = T_AMB,
    p_init = P_PUMP_E,   
    m_flow_init = M_CO2
  );
  
  Steps.Components.Regulator regulator_LTR_hot(
    T_init = from_degC(578),
    p_init = P_TURBINE_E,   
    m_flow_init = M_CO2
  );
  
  Steps.Components.Regulator regulator_LTR_cold(
    T_init = T_AMB,
    p_init = P_PUMP_E,   
    m_flow_init = M_CO2
  );*/
  
  Steps.Components.FanCooler fan_cooler(
    T_output = T_COOLER_E,      
    inlet.p(start = P_LTR_H_E),
    inlet.h_outflow(start = H_LTR_H_E),
    outlet.p(start = P_COOLER_E),
    outlet.h_outflow(start = H_COOLER_E)
  );
  
  // use pump as compressor
  Steps.Components.Pump pump(
    p_outlet = P_PUMP_E, 
    eta = eta_main_compressor,/*
    outlet.p(start = P_PUMP_E),
    outlet.h_outflow(start = H_PUMP_E),*/
    inlet.p(start = P_PUMP_I),
    inlet.h_outflow(start = H_PUMP_I)
  );   
  
  Steps.Components.Pump recom_pump(
    p_outlet = P_PUMP_E, 
    eta = eta_bypass_compressor/*,
    inlet.p(start = P_PUMP_I),
    inlet.h_outflow(start = H_PUMP_I),
    outlet.p(start = P_PUMP_E),
    outlet.h_outflow(start = H_PUMP_E)*/
  );
  
  Components.PCHECImpl HTR(
    name = "HTR", 
    phi = phi, 
    d_c = d_c,
    pitch = pitch,     
    N_ch = N_ch,    
    N_seg = N_seg,
    length_cell = length_cell,
    sim_param(log_level = 4),
    
    inlet_hot.p(start = P_TURBINE_E),
    inlet_hot.h_outflow(start = H_TURBINE_E),
    inlet_hot.m_flow(start = M_CO2),
    
    inlet_cold.p(start = P_LTR_C_E),
    inlet_cold.h_outflow(start = H_LTR_C_E),
    inlet_cold.m_flow(start = M_CO2),
    
    outlet_hot.p(start = P_HTR_H_E),
    outlet_hot.h_outflow(start = H_HTR_H_E),
    outlet_hot.m_flow(start = M_CO2)   
    /*
    
    bc_hot_in(p = P_TURBINE_E, T = from_degC(578), h = PropsSI("H", "T", from_degC(578), "P", P_TURBINE_E, PBMedia.mediumName), mdot = M_CO2),    
    bc_cold_in(p = P_PUMP_E, T = from_degC(151), h = PropsSI("H", "T", from_degC(151), "P", P_PUMP_E, PBMedia.mediumName), mdot = M_CO2)*/   
        
  );
  
  Components.PCHECImpl LTR(
    name = "LTR", 
    phi = phi, 
    d_c = d_c,
    pitch = pitch,     
    N_ch = N_ch,    
    N_seg = N_seg,
    length_cell = length_cell,
    sim_param(log_level = 4),
    
    inlet_hot.p(start = P_HTR_H_E),
    inlet_hot.h_outflow(start = H_HTR_H_E),
    inlet_cold.p(start = P_PUMP_E),
    inlet_cold.h_outflow(start = H_PUMP_E)/*,   
    
    bc_hot_in(p = P_TURBINE_E, T = from_degC(156), h = PropsSI("H", "T", from_degC(156), "P", P_TURBINE_E, PBMedia.mediumName), mdot = M_CO2),    
    bc_cold_in(p = P_PUMP_E, T = from_degC(62), h = PropsSI("H", "T", from_degC(62), "P", P_PUMP_E, PBMedia.mediumName), mdot = M_CO2)*/
  );

  Steps.Components.TemperatureOutput temp_out(
    T_start = T_HEATER_E,
    T_stop = T_HEATER_E,
    dT_step = 50,
    t_sim_duration = stop_time
  );
  
  Steps.Components.PCMHeater pcm_heater(
    mdot_init = M_CO2,
    T_input = T_HEATER_E,
    /*
    init_outlet = true,
    p_init_outlet = P_PUMP_E,
    T_init_outlet = T_AMB
    
    
    inlet.p(start = P_PUMP_E),
    inlet.h_outflow(start = H_PUMP_E),*/
    outlet.p(start = P_HEATER_E),
    outlet.h_outflow(start = H_HEATER_E)    
    
  );  
   
  Steps.Components.Turbine turbine(    
    p_out = P_TURBINE_E,
    eta = eta_turbine, 
    inlet.p(start = P_HEATER_E),
    inlet.h_outflow(start = H_HEATER_E),
    outlet.p(start = P_TURBINE_E),
    outlet.h_outflow(start = H_TURBINE_E)    
  ); 
  
  Steps.Components.Splitter splitter(
    split_ratio = splitter_split_ratio/*,
    outlet.p(start = P_PUMP_I),
    outlet.h_outflow(start = H_PUMP_I),
    inlet.p(start = P_PUMP_I),
    inlet.h_outflow(start = H_PUMP_I),
    outlet_split.p(start = P_PUMP_I),
    outlet_split.h_outflow(start = H_PUMP_I)*/
  );
  
  Steps.Components.Merger merger(
    //T_init = T_AMB,
    //p_init = P_ATM,
    outlet.p(start = P_PUMP_E),
    outlet.h_outflow(start = H_PUMP_E)/*,
    inlet.p(start = P_PUMP_E),
    inlet.h_outflow(start = H_PUMP_E),    
    inlet_merge.p(start = P_PUMP_E),
    inlet_merge.h_outflow(start = H_PUMP_E)*/
  );
  
  //total efficiency
  Real eta_total;


equation
  
  // connect(regulator.outlet, turbine.inlet);
  
  connect(turbine.outlet, HTR.inlet_hot);
  
  // connect(turbine.outlet, regulator_HTR_hot.inlet);  
  // connect(regulator_HTR_hot.outlet, HTR.inlet_hot);
  
  connect(HTR.outlet_hot,LTR.inlet_hot);
  // connect(HTR.outlet_hot,regulator_LTR_hot.inlet);
  // connect(regulator_LTR_hot.outlet, LTR.inlet_hot); 
 
  connect(LTR.outlet_hot, splitter.inlet);
  //recompression loop
  connect(splitter.outlet_split, recom_pump.inlet);
  connect(recom_pump.outlet, merger.inlet_merge);  
  
  connect(splitter.outlet, fan_cooler.inlet);
  
  connect(fan_cooler.outlet, pump.inlet);
  
  connect(pump.outlet, LTR.inlet_cold);  
  
  // connect(pump.outlet, regulator_LTR_cold.inlet);  
  
  // connect(regulator_LTR_cold.outlet, LTR.inlet_cold);  
   
  connect(LTR.outlet_cold, merger.inlet);
  
  connect(merger.outlet, HTR.inlet_cold);
  
  // connect(merger.outlet, regulator_HTR_cold.inlet);
  
  // connect(regulator_HTR_cold.outlet, HTR.inlet_cold);
  
  connect(HTR.outlet_cold, pcm_heater.inlet);    
  
  //connect(temp_out.y, pcm_heater.T_input);
  // connect(pcm_heater.outlet, regulator.inlet);
  connect(pcm_heater.outlet, turbine.inlet);

// algorithm
  eta_total = if initial() then 0 else (turbine.W_turbine - pump.W_comp - recom_pump.W_comp) / pcm_heater.Q * 100;
annotation(
    experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-6, Interval = 0.002),
    __OpenModelica_commandLineOptions = "--matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection -d=initialization,NLSanalyticJacobian");
end OffDesignRCBCycle_v2;
