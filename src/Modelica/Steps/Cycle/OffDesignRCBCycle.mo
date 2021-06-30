within Steps.Cycle;

model OffDesignRCBCycle
  "OffDesign Brayton Cycle with Recuperator + Recompression"  
// Adjustable parameters
  import Util = Utilities.Util;
  import Modelica.SIunits.Conversions.{from_degC, from_deg};
  import Steps.Utilities.CoolProp.PropsSI;
  
  replaceable package PBMedia = Steps.Media.SCO2;
  parameter Modelica.SIunits.Time stop_time = 1.0 "time length of the experiment";
  
  parameter Real M_CO2 = 51.91;    //# Co2 flow rate, kg/s
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
  parameter Modelica.SIunits.Temperature T_TURBINE_E = from_degC(578);
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
  parameter Real splitter_split_ratio = 0.39; 
  
  parameter Modelica.SIunits.Angle phi = from_deg((180 - 108) /2) "phi of zigzag - High Temp range for zigzag in 4.1.2 Meshram [2016]"; // agree with Hal's steps instead of using 4.1.2 of Meshram [2016]

  parameter Integer N_ch_LTR = integer(80e3);  
  
  parameter Integer N_ch_HTR = integer(200e3);  
  
  parameter Modelica.SIunits.Length d_c = 2e-3 "diameter of the channel";

  parameter Modelica.SIunits.Length pitch = 12e-3 "pitch length of zigzag channel";
  
  parameter Modelica.SIunits.Length length_LTR = 360e-3;
  
  parameter Modelica.SIunits.Length length_HTR = 360e-3;

  //parameter Modelica.SIunits.Length length_cell = length_pche / N_seg "length of the discretized cell in a channel";

  parameter Integer N_seg = 50 "number of cells/segment for the discretization of a channel";  
  
  parameter Modelica.SIunits.ReynoldsNumber Re_hot_start = 14.5e3 "Hot stream's start value of Reynolds Number, used to increase convergence";
  
  parameter Modelica.SIunits.ReynoldsNumber Re_cold_start = 14.5e3 "Cold stream's start value of Reynolds Number, used to increase convergence";    
  
  parameter Boolean SourceFixed_hot = true;
  
  parameter Boolean SourceFixed_cold = true;  
  /*
  Steps.Components.Regulator regulator(
    T_init = T_AMB,
    p_init = P_PUMP_E,   
    m_flow_init = M_CO2,
    inlet.p(start = P_PUMP_E),
    inlet.h_outflow(start = H_TURBINE_I),
    outlet.p(start = P_PUMP_E),
    outlet.h_outflow(start = H_TURBINE_I)
  );
  */
  Steps.Components.FanCooler fan_cooler(
    T_output = T_COOLER_E,
    outlet.p(start = P_PUMP_I),
    outlet.h_outflow(start = H_PUMP_I),
    inlet.p(start = P_LTR_H_E),
    inlet.h_outflow(start = H_LTR_H_E)
  );
  
  // use pump as compressor
  Steps.Components.Pump pump(
    p_outlet = P_PUMP_E, 
    eta = eta_main_compressor,
    inlet.p(start = P_PUMP_I),
    inlet.h_outflow(start = H_PUMP_I),
    outlet.p(start = P_PUMP_E),
    outlet.h_outflow(start = H_PUMP_E)    
  );   
  
  Steps.Components.Pump recom_pump(
    p_outlet = P_PUMP_E, 
    eta = eta_bypass_compressor,
    inlet.p(start = P_PUMP_I),
    inlet.h_outflow(start = H_PUMP_I),
    outlet.p(start = P_PUMP_E),
    outlet.h_outflow(start = H_PUMP_E)      
  );
  
  Components.PCHeatExchanger HTR(
    phi = phi, 
    d_c = d_c,
    pitch = pitch,     
    N_ch = N_ch_HTR,
    Re_cold_start = Re_cold_start,
    Re_hot_start = Re_hot_start,
    N_seg = N_seg,
    length = length_HTR,
    
    inlet_hot.p(start = P_TURBINE_E),
    inlet_hot.h_outflow(start = PropsSI("H", "T", T_TURBINE_E, "P", P_TURBINE_E, PBMedia.mediumName)),
    inlet_hot.m_flow(start = M_CO2),
    
    inlet_cold.p(start = P_LTR_C_E),
    inlet_cold.h_outflow(start = PropsSI("H", "T", T_LTR_C_E, "P", P_LTR_C_E, PBMedia.mediumName)),
    inlet_cold.m_flow(start = M_CO2),
    
    outlet_cold.p(start = P_LTR_C_E),
    outlet_hot.p(start = P_HTR_H_E),
    //outlet_cold.h_outflow(start = PropsSI("H", "T", T_LTR_C_E, "P", P_LTR_C_E, PBMedia.mediumName)),
    //outlet_cold.m_flow(start = M_CO2),
    
    p_start_hot = P_TURBINE_E,
    T_start_hot = T_TURBINE_E,
    //h_start_hot = H_PUMP_I,
    
    p_start_cold = P_LTR_C_E,
    T_start_cold = T_LTR_C_E,    
    //h_start_cold = H_PUMP_E,
    
    ByInlet_hot = SourceFixed_hot,
    ByInlet_cold = SourceFixed_cold     
  );
  
  Components.PCHeatExchanger LTR(
    phi = phi, 
    d_c = d_c,
    pitch = pitch,     
    N_ch = N_ch_LTR,
    Re_cold_start = Re_cold_start,
    Re_hot_start = Re_hot_start,
    N_seg = N_seg,
    length = length_HTR, 
     
    inlet_hot.p(start = P_HTR_H_E),
    inlet_hot.h_outflow(start = H_HTR_H_E),
    inlet_hot.m_flow(start = M_CO2),
    
    inlet_cold.p(start = P_PUMP_E),
    inlet_cold.h_outflow(start = H_PUMP_E),
    inlet_cold.m_flow(start = M_CO2 * (1 - splitter_split_ratio)),  
    
    outlet_hot.p(start = P_LTR_H_E),
    outlet_hot.h_outflow(start = H_LTR_H_E),
    
    outlet_cold.p(start = P_PUMP_E),
    
    p_start_hot = P_HTR_H_E,   
    T_start_hot = T_HTR_H_E,  
    //h_start_hot = H_PUMP_I,
    
    p_start_cold = P_PUMP_E,
    T_start_cold = T_PUMP_E,    
    //h_start_cold = H_PUMP_E,
    
    ByInlet_hot = SourceFixed_hot,
    ByInlet_cold = SourceFixed_cold    
  );

  Steps.Components.TemperatureOutput temp_out(
    T_start = T_HEATER_E,
    T_stop = T_HEATER_E,
    dT_step = 50,
    t_sim_duration = stop_time
  );
  
  Steps.Components.PCMHeater pcm_heater(
    mdot_init = M_CO2,
    T_output_set = T_HEATER_E,
    inlet.p(start = P_HEATER_E),
    //inlet.h_outflow(start = H_PUMP_E),
    outlet.p(start = P_HEATER_E),
    outlet.h_outflow(start = H_HEATER_E)    
  );  
   
  Steps.Components.Turbine turbine(    
    p_out = P_TURBINE_E,
    eta = eta_turbine,
    inlet.p(start = P_HEATER_E),
    inlet.h_outflow(start = H_HEATER_E)
    
    //outlet.p(start = P_TURBINE_E)
    //outlet.h_outflow(start = H_PUMP_I)    
  ); 
  
  Steps.Components.Splitter splitter(
    split_ratio = splitter_split_ratio,
    outlet.p(start = P_LTR_H_E),
    outlet.h_outflow(start = H_LTR_H_E),
    inlet.p(start = P_LTR_H_E),
    //inlet.h_outflow(start = H_PUMP_I),
    outlet_split.p(start = P_LTR_H_E)
    //outlet_split.h_outflow(start = H_PUMP_I)
  );
  
  Steps.Components.Mixer mixer(
    //T_init = T_AMB,
    //p_init = P_ATM
    inlet.p(start = P_LTR_C_E),
    inlet.h_outflow(start = H_PUMP_E),    
    outlet.p(start = P_LTR_C_E),
    outlet.h_outflow(start = H_PUMP_E),
    inlet_mix.p(start = P_PUMP_E),
    inlet_mix.h_outflow(start = H_PUMP_E)
  );
  
  //total efficiency
  Real eta_total;

equation
  
  //connect(regulator.outlet, turbine.inlet);
  
  connect(turbine.outlet, HTR.inlet_hot);
  
  connect(HTR.outlet_hot, LTR.inlet_hot); 
 
  connect(LTR.outlet_hot, splitter.inlet);
  //recompression loop
  connect(splitter.outlet_split, recom_pump.inlet);
  connect(recom_pump.outlet, mixer.inlet_mix);  
  
  connect(splitter.outlet, fan_cooler.inlet);
  
  connect(fan_cooler.outlet, pump.inlet);
  
  connect(pump.outlet, LTR.inlet_cold);  
   
  connect(LTR.outlet_cold, mixer.inlet);
  
  connect(mixer.outlet, HTR.inlet_cold);
  
  connect(HTR.outlet_cold, pcm_heater.inlet);    
  
  //connect(temp_out.y, pcm_heater.T_input);
  connect(pcm_heater.outlet, turbine.inlet);

// algorithm
  eta_total = if initial() then 0 else (turbine.W_turbine - pump.W_comp - recom_pump.W_comp) / pcm_heater.Q * 100;
annotation(
    experiment(StartTime = 0, StopTime = 1, Interval = 1, Tolerance = 1e-6),
    __OpenModelica_commandLineOptions = "--matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection -d=initialization,NLSanalyticJacobian");
end OffDesignRCBCycle;
