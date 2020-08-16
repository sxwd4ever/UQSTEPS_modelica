within Steps;

model RECUPBCycle
  "Brayton Cycle with Recuperator"  
// Adjustable parameters
  import Util = Utilities.Util;
  import Steps.Utilities.Util.checkPortType;
  import Steps.Interfaces.PortType;
  
  parameter Modelica.SIunits.Pressure P_ATM = 101325; // Pa
  parameter Modelica.SIunits.Temperature T_AMB = Modelica.SIunits.Conversions.from_degC(15) "Ambinent temperature";
  
  parameter Modelica.SIunits.Time stop_time = 1.0 "time length of the experiment";
  
  parameter Real M_CO2 = 51.91;    //# Co2 flow rate, kg/s
  parameter Modelica.SIunits.Pressure P_PUMP_I = 8.0 * 1e6 ; //# Pa
  parameter Modelica.SIunits.Pressure P_PUMP_E = 20.0 * 1e6;    //# Pump exit pressure, Pa

  // efficiency of compressor, bypass_compressor and turbine
  parameter Real eta_compressor = 0.89;
  
  parameter Real eta_turbine = 0.9;
  
  // effectiveness for two recuperators
  // these two parameters will effect the difference between crec_in.T and hot_in.T of recuperator
  parameter Real eta_recuperator = 0.99;   

  Steps.Components.Regulator regulator(
    T_init = Modelica.SIunits.Conversions.from_degC(10),
    p_init = P_PUMP_E,   
    m_flow_init = M_CO2,
    outlet.PT = PortType.pT_fixed
  );
  
  parameter Modelica.SIunits.TemperatureDifference DT_COOLER = 18.0;
  
  Steps.Components.FanCooler fan_cooler(
    T_amb = T_AMB,
    delta_T = DT_COOLER,
    outlet.PT = PortType.T_fixed
  );
  
  // use pump as compressor
  Steps.Components.Pump pump(
    p_outlet = P_PUMP_E, 
    eta = eta_compressor,
    outlet.PT = PortType.p_fixed
  ); 
    
  /*
  // temperature recuperator
  Steps.Components.Recuperator recup(
    eta = eta_recuperator  
  );  
  */ 
   
  Components.MockPCHeatExchanger recup(
    phi = Modelica.SIunits.Conversions.from_deg(45), 
    Re_design = 5000,
    d_c = 1.51 * 1e-3,
    T_hot_in = Modelica.SIunits.Conversions.from_degC(451),
    T_cool_in = Modelica.SIunits.Conversions.from_degC(41),
    p_hot = 9 * 1e6,
    p_cool = 20 * 1e6,
    m_dot_hot = 8.3,
    m_dot_cool = 8.3,
    pitch = 24.6 * 1e-3,
    length_cell = 3e-3,
    N_seg = 1
    //inlet_hot.PT = PortType.pT_fixed,
    //outlet_cool.PT = PortType.pT_fixed
  );  
   
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
    eta = eta_turbine,
    outlet.PT = PortType.p_fixed
  ); 
  

  //total efficiency
  Real eta_total;

initial algorithm
  
  //check the compatibility between two ports
  // required to make the equations solveable  

  checkPortType(regulator.outlet.PT, turbine.inlet.PT, name_conn = "regulator->turbine");
  
  checkPortType(turbine.outlet.PT, recup.inlet_hot.PT, name_conn = "turbine -> recup");
  
  checkPortType(recup.outlet_hot.PT, fan_cooler.inlet.PT, name_conn = "recup -> fan_cooler");
  
  checkPortType(fan_cooler.outlet.PT, pump.inlet.PT, name_conn = "fan_cooler -> pump");
  
  checkPortType(pump.outlet.PT, recup.inlet_cool.PT, name_conn = "pump -> recup");
  
  checkPortType(recup.outlet_cool.PT, pcm_heater.inlet.PT, name_conn = "recup -> pcm_heater");    

  checkPortType(pcm_heater.outlet.PT, regulator.inlet.PT, name_conn = "pcm_heater -> regulator");

equation
  
  connect(regulator.outlet, turbine.inlet);
  
  connect(turbine.outlet, recup.inlet_hot);
  
  connect(recup.outlet_hot, fan_cooler.inlet);
  
  connect(fan_cooler.outlet, pump.inlet);
  
  connect(pump.outlet, recup.inlet_cool);
  
  connect(recup.outlet_cool, pcm_heater.inlet);    
  
  connect(temp_out.y, pcm_heater.T_input);
  
  connect(pcm_heater.outlet, regulator.inlet);

algorithm
  //eta_total := if initial() then 0 else (turbine.W_turbine - pump.W_comp) / pcm_heater.Q * 100;
  eta_total := 1.0;
annotation(
    experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-6, Interval = 0.002),
    __OpenModelica_commandLineOptions = "--matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection -d=initialization,NLSanalyticJacobian");
end RECUPBCycle;
