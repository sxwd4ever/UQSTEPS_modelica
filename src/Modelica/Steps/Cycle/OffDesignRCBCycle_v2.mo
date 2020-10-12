within Steps.Cycle;

model OffDesignRCBCycle_v2
  "OffDesign Brayton Cycle with Recuperator + Recompression"  
// Adjustable parameters
  import Util = Utilities.Util;
  import Modelica.SIunits.Conversions.{from_degC, from_deg};
  import Steps.Utilities.CoolProp.PropsSI;
  
  replaceable package PBMedia = Steps.Media.SCO2;
  
  parameter Modelica.SIunits.Pressure P_ATM = 101325; // Pa
  parameter Modelica.SIunits.Temperature T_AMB = from_degC(15) "Ambinent temperature";
  
  parameter Modelica.SIunits.Time stop_time = 1.0 "time length of the experiment";
  
  parameter Real M_CO2 = 100;    //# Co2 flow rate, kg/s
  parameter Modelica.SIunits.Pressure P_PUMP_I = 8.0 * 1e6 ; //# Pa
  parameter Modelica.SIunits.Pressure P_PUMP_E = 20.0 * 1e6;    //# Pump exit pressure, Pa
  parameter Modelica.SIunits.Pressure P_TURBINE_E = P_PUMP_I;
  
  parameter Modelica.SIunits.SpecificEnthalpy H_PUMP_I = PropsSI("H", "T", T_AMB, "P", P_PUMP_I, PBMedia.mediumName); 
  parameter Modelica.SIunits.SpecificEnthalpy H_PUMP_E = PropsSI("H", "T", T_AMB, "P", P_PUMP_E, PBMedia.mediumName); 
  parameter Modelica.SIunits.SpecificEnthalpy H_TURBINE_I = PropsSI("H", "T", T_HEATER_E, "P", P_PUMP_E, PBMedia.mediumName); 
  
  parameter Modelica.SIunits.Temp_K T_HEATER_E = from_degC(700);
  
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
  parameter Real splitter_split_ratio = 0.4; 
  
  parameter Modelica.SIunits.Angle phi = from_deg((180 - 108) /2) "phi of zigzag - High Temp range for zigzag in 4.1.2 Meshram [2016]"; // agree with Hal's steps instead of using 4.1.2 of Meshram [2016]

  parameter Integer N_ch = integer(80e3);  
  
  parameter Modelica.SIunits.Length d_c = 2e-3 "diameter of the channel";

  parameter Modelica.SIunits.Length pitch = 12e-3 "pitch length of zigzag channel";

  parameter Modelica.SIunits.Length length_cell = 12e-3 "length of the discretized cell in a channel";

  parameter Integer N_seg = Utilities.CoolProp.MAX_SEG_NUM "number of cells/segment for the discretization of a channel";  
  
  parameter Modelica.SIunits.ReynoldsNumber Re_hot_start = 14.5e3 "Hot stream's start value of Reynolds Number, used to increase convergence";
  
  parameter Modelica.SIunits.ReynoldsNumber Re_cold_start = 14.5e3 "Cold stream's start value of Reynolds Number, used to increase convergence";    

  parameter Modelica.SIunits.TemperatureDifference DT_COOLER = 18.0;
  
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
  Steps.Components.FanCooler fan_colder(
    T_amb = T_AMB,
    delta_T = DT_COOLER,
    outlet.p(start = P_PUMP_I)
    //outlet.h_outflow(start = H_PUMP_I),
    //inlet.p(start = P_PUMP_I),
    //inlet.h_outflow(start = H_PUMP_I),
  );
  
  // use pump as compressor
  Steps.Components.Pump pump(
    p_outlet = P_PUMP_E, 
    eta = eta_main_compressor
    //inlet.p(start = P_PUMP_I),
    //inlet.h_outflow(start = H_PUMP_I),
    //outlet.p(start = P_PUMP_E)
    //outlet.h_outflow(start = H_PUMP_E)    
  );   
  
  Steps.Components.Pump recom_pump(
    p_outlet = P_PUMP_E, 
    eta = eta_bypass_compressor
    //inlet.p(start = P_PUMP_I),
    //inlet.h_outflow(start = H_PUMP_I),
    //outlet.p(start = P_PUMP_E)
    //outlet.h_outflow(start = H_PUMP_E)      
  );
  
  Components.PCHECImpl recup_high(
    phi = phi, 
    d_c = d_c,
    pitch = pitch,     
    N_ch = N_ch,    
    N_seg = N_seg,
    length_cell = length_cell,  
    
    inlet_hot.p(start = P_PUMP_I, nominal = P_PUMP_I),
    inlet_hot.h_outflow(start = h_HTR_hot_in, nominal = h_HTR_hot_in),
    inlet_hot.m_flow(start = M_CO2, nominal = M_CO2),
    
    inlet_cold.p(start = P_PUMP_E, nominal = P_PUMP_E),
    inlet_cold.h_outflow(start = h_HTR_cold_in, nominal = h_HTR_cold_in),
    inlet_cold.m_flow(start = M_CO2, nominal = M_CO2) 
        
  );
  
  Components.PCHECImpl recup_low(
    phi = phi, 
    d_c = d_c,
    pitch = pitch,     
    N_ch = N_ch,    
    N_seg = N_seg,
    length_cell = length_cell,    
    
    inlet_hot.p(start = P_TURBINE_E, nominal = P_TURBINE_E),
    inlet_hot.h_outflow(start = h_LTR_hot_in, nominal = h_LTR_hot_in),
    inlet_hot.m_flow(start = M_CO2, nominal = M_CO2),
    
    inlet_cold.p(start = P_PUMP_E, nominal = P_PUMP_E),
    inlet_cold.h_outflow(start = h_LTR_cold_in, nominal = h_LTR_cold_in),
    inlet_cold.m_flow(start = M_CO2, nominal = M_CO2)   

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
    //inlet.p(start = P_PUMP_E),
    //inlet.h_outflow(start = H_PUMP_E),
    outlet.p(start = P_PUMP_E),
    outlet.h_outflow(start = H_TURBINE_I)    
  );  
   
  Steps.Components.Turbine turbine(    
    p_out = P_TURBINE_E,
    eta = eta_turbine
    //outlet.p(start = P_TURBINE_E)
    //outlet.h_outflow(start = H_PUMP_I)    
  ); 
  
  Steps.Components.Splitter splitter(
    split_ratio = splitter_split_ratio
    //outlet.p(start = P_PUMP_I),
    //outlet.h_outflow(start = H_PUMP_I),
    //inlet.p(start = P_PUMP_I),
    //inlet.h_outflow(start = H_PUMP_I),
    //outlet_split.p(start = P_PUMP_I)
    //outlet_split.h_outflow(start = H_PUMP_I)
  );
  
  Steps.Components.Merger merger(
    T_init = T_AMB,
    p_init = P_ATM
    //inlet.p(start = P_PUMP_E),
    //inlet.h_outflow(start = H_PUMP_E),    
    //outlet.p(start = P_PUMP_E),
    //outlet.h_outflow(start = H_PUMP_E),
    //inlet_merge.p(start = P_PUMP_E)
    //inlet_merge.h_outflow(start = H_PUMP_E)
  );
  
  //total efficiency
  Real eta_total;

protected
  parameter Modelica.SIunits.SpecificEnthalpy h_HTR_hot_in = PropsSI("H", "T", from_degC(156), "P", P_PUMP_I, PBMedia.mediumName);
  parameter Modelica.SIunits.SpecificEnthalpy h_HTR_cold_in = PropsSI("H", "T", from_degC(63), "P", P_PUMP_E, PBMedia.mediumName);

  parameter Modelica.SIunits.SpecificEnthalpy h_LTR_hot_in = PropsSI("H", "T", from_degC(578), "P", P_TURBINE_E, PBMedia.mediumName);
  parameter Modelica.SIunits.SpecificEnthalpy h_LTR_cold_in = PropsSI("H", "T", from_degC(151), "P", P_PUMP_E, PBMedia.mediumName);
equation
  
  //connect(regulator.outlet, turbine.inlet);
  
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
  
  //connect(temp_out.y, pcm_heater.T_input);
  connect(pcm_heater.outlet, turbine.inlet);

// algorithm
  eta_total = if initial() then 0 else (turbine.W_turbine - pump.W_comp - recom_pump.W_comp) / pcm_heater.Q * 100;
annotation(
    experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-6, Interval = 0.002),
    __OpenModelica_commandLineOptions = "--matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection -d=initialization,NLSanalyticJacobian");
end OffDesignRCBCycle_v2;
