within Steps.Test.TPComponents;

model TestTP_PCHE_PHE
  "Test for ThermoPower based PCHE model against Pinjarra Hill Experimental Data - Off design"  
    
  import Modelica.SIunits.Conversions.{from_degC, from_deg};
  import Modelica.SIunits.{Temperature, Pressure, SpecificEnthalpy};
  import Util = Utilities.Util;
  import Steps.Utilities.CoolProp.PropsSI;  
  import Steps.Components.{PCHEGeoParam};
  import Steps.Model.{PBConfiguration, PBConfigs, SimParam, EntityConfig, EntityGeoParam, EntityThermoParam, ThermoState, HEBoundaryCondition} ;
  import Model.PBConfiguration;
  import ThermoPower.Choices.Init.Options;
  import ThermoPower.System;
  import ThermoPower.Gas;
  import Steps.Components.KimCorrelations;
  import Steps.Components.MaterialConductivity;    

  // package medium_hot = Modelica.Media.IdealGases.SingleGases.CO2;
  package medium_hot = Steps.Media.ThermiaOilD;
  // package medium_cold = Modelica.Media.IdealGases.SingleGases.CO2;
  package medium_cold = Steps.Media.SCO2;  
  
  // geometry parameters
  constant  Real pi        = Modelica.Constants.pi;
  parameter Integer N_ch   = integer(2400) "channel number";
  parameter Integer N_seg  = 10 "number of segments in one tube";
  parameter SI.Length D_ch = 1.72e-3  "channel diameter, semi circular tube";
  parameter SI.Length r_ch = D_ch / 2 "channel radius";
  // parameter SI.Length L_fp = 270e-3 "channel flow path length";  
  // parameter SI.Length L_fp = (200 + 186) * 1e-3 "equivalent valid channel flow path length";  
  parameter SI.Length L_fp    = 250 * 1e-3 "equivalent valid channel flow path length";
  parameter SI.Length L_pitch = 12e-3 "pitch length";
  parameter Real a_phi        = 36 "pitch angle degree";
  parameter SI.Length H_ch    = 4.17e-3 "Height of the solid domain, containing one cold tube and one hot tube";
  parameter SI.Length W_ch    = 2.3e-3 "Width of the solid domain";
  parameter SI.Length T_wall  = 0.51e-3 "Wall thickness";
  parameter SI.Length L_wall  = 420e-3 "Length of wall, not necessarily equals to length of flow path";
  parameter SI.Area A         = pi * r_ch ^2 / 2 "Area of cross section of semi circular tube";

  // boundary conditon
  
  // zigzag higher T
  // parameter SI.Velocity u_hot_in = 7.564 "hot inlet velocity m/s";
  // parameter SI.Velocity u_cold_in = 1.876 "cold inlet velocity m/s";
  parameter SI.Pressure p_hot_in      = from_bar(10) "hot inlet pressure - Irrelevant for Incompressible Thermia Oil" ;
  parameter SI.Pressure p_cold_in     = 12.567478e6 "cold inlet pressure";
  parameter SI.Temperature T_hot_in   = from_degC(103.222748) "hot inlet temperature, K";
  parameter SI.Temperature T_hot_out  = from_degC(96.145935) "cold outlet temperature, K";
  parameter SI.Temperature T_cold_in  = from_degC(28.910231) "cold inlet temperature, K";
  parameter SI.Temperature T_cold_out = from_degC(99.666342) "cold outlet temperature, K";

  parameter Real Cf_C1 = 1, Cf_C2 = 1, Cf_C3 = 1;
  parameter Real use_rho_bar  = -1.0;  
  parameter Real rho_bar_hot  = 1.0;
  parameter Real rho_bar_cold = 1.0;  
  
  // meshram's cp and rho for alloy Inconel 617
  // parameter Modelica.SIunits.Density rho_wall = 8360 "density of wall, kg/m3";
  // parameter Modelica.SIunits.SpecificHeatCapacity cp_wall = 417 "cp of wall, J/kg-K";  
   
  // Stainless 316, 316L, 317, 317L
  parameter Modelica.SIunits.Density rho_wall = 8030 "density of wall, kg/m3";
  parameter Modelica.SIunits.SpecificHeatCapacity cp_wall = 485 "cp of wall, J/kg-K";  
  // thermal conductivity (T in K) https://www.theworldmaterial.com/aisi-316-ss316-stainless-steel-properties-composition/
  // parameter Real table_k_metalwall[:,:] = [20, 12.1; 100, 16.3; 500, 21.5];
  parameter Real table_k_metalwall[:,:] = [293.15, 12.1; 373.15, 16.3; 773.15, 21.5];

  parameter SI.MassFlowRate mdot_hot_in  = 1.218944;
  parameter SI.MassFlowRate mdot_cold_in = 0.0854299999999999;

  parameter Model.PCHEConfig cfg(
    redeclare package medium_main   = medium_cold,
    redeclare package medium_heater = medium_hot,
    N_ch = 1200
  );

  // set the parameters accordingly
  parameter Model.HeatExchangerConfig cfg_HE     = cfg.cfg_PCHE;
  parameter Model.ThermoState st_source_hot      = cfg_HE.cfg_hot.st_in;
  parameter Model.ThermoState st_sink_hot        = cfg_HE.cfg_hot.st_out;
  parameter Model.ThermoState st_source_cold     = cfg_HE.cfg_cold.st_in;
  parameter Model.ThermoState st_sink_cold       = cfg_HE.cfg_cold.st_out;
  parameter Integer N_seg_HE                     = cfg.cfg_PCHE.cfg_hot.geo_path.N_seg;

  //Components
  // for transient simulation, set initOpt = steadyState, which can be used for steadyState simulation too. 
  inner ThermoPower.System system(allowFlowReversal = false, initOpt = ThermoPower.Choices.Init.Options.steadyState) annotation(
    Placement(transformation(extent = {{80, 80}, {100, 100}})));  
  
  parameter Boolean SSInit = false "Steady-state initialization";
  
  ThermoPower.Gas.SourceMassFlow source_cold(
    redeclare package Medium = medium_cold, 
    T  = st_source_cold.T,
    p0 = st_source_cold.p,
    //use_in_T = false, 
    w0 = st_source_cold.mdot,
    gas(
      p(nominal = st_source_cold.p), 
      T(nominal = st_source_cold.T))) 
  annotation(
    Placement(transformation(origin = {0, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 270)));
  
  ThermoPower.Gas.SinkPressure sink_cold(
    redeclare package Medium = medium_cold, 
    p0 = st_sink_cold.p, 
    T = st_sink_cold.T) 
  annotation(
    Placement(transformation(origin = {0, -80}, extent = {{-10, -10}, {10, 10}}, rotation = 270)));
 
  ThermoPower.Gas.SinkPressure sink_hot(
    redeclare package Medium = medium_hot,
    T  = st_sink_hot.T,
    p0 = st_sink_hot.p)
  annotation(
    Placement(transformation(extent = {{60, -10}, {80, 10}}, rotation = 0)));
  
  ThermoPower.Gas.SourceMassFlow source_hot(
    redeclare package Medium = medium_hot, 
    T  = st_source_hot.T,
    p0 = st_source_hot.p,
    w0 = st_source_hot.mdot,
    //use_in_T = false,
    gas(
      p(nominal = st_source_hot.p), 
      T(nominal = st_source_hot.T))) 
  annotation(
    Placement(transformation(extent = {{-70, -10}, {-50, 10}}, rotation = 0))); 

  Steps.TPComponents.GasStateReader r_HE_hin(redeclare package Medium = medium_hot) annotation(
    Placement(visible = true, transformation(origin = {66, -56}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
  Steps.TPComponents.GasStateReader r_HE_hout(redeclare package Medium = medium_hot) annotation(
    Placement(visible = true, transformation(origin = {96, 44}, extent = {{-6, -6}, {6, 6}}, rotation = 180)));
  Steps.TPComponents.GasStateReader r_HE_cin(redeclare package Medium = medium_cold) annotation(
    Placement(visible = true, transformation(origin = {-22, 50}, extent = {{-6, -6}, {6, 6}}, rotation = -90)));
  Steps.TPComponents.GasStateReader r_HE_cout(redeclare package Medium = medium_cold) annotation(
    Placement(visible = true, transformation(origin = {-32, 20}, extent = {{-6, -6}, {6, 6}}, rotation = 180)));

  Steps.TPComponents.PCHE HE(
    redeclare package FluidMedium = medium_cold, 
    redeclare package FlueGasMedium = medium_hot, 
     
  // with GnielinskiHeatTransferFV Correlation    
    redeclare replaceable model HeatTransfer_F = Steps.TPComponents.GnielinskiHeatTransferFV(),
    redeclare replaceable model HeatTransfer_G = Steps.TPComponents.GnielinskiHeatTransferFV(),
    gasFlow(
      heatTransfer(
        pitch       = cfg_HE.cfg_hot.l_pitch,
        phi         = cfg_HE.cfg_hot.a_phi,         
        Cf_C1       = Cf_C1, 
        Cf_C2       = Cf_C2, 
        Cf_C3       = Cf_C3, 
        use_rho_bar = use_rho_bar,
        rho_bar     = rho_bar_hot,
        gamma_min   = 100,
        useAverageTemperature = false)),
    fluidFlow(
      heatTransfer(
        pitch       = cfg_HE.cfg_cold.l_pitch, 
        phi         = cfg_HE.cfg_cold.a_phi,         
        Cf_C1       = Cf_C1, 
        Cf_C2       = Cf_C2, 
        Cf_C3       = Cf_C3, 
        use_rho_bar = use_rho_bar, 
        rho_bar     = rho_bar_cold,
        gamma_min   = 100,
        useAverageTemperature = false)
    ),  
    redeclare model HeatExchangerTopology = ThermoPower.Thermal.HeatExchangerTopologies.CounterCurrentFlow, 
    cfg               = cfg_HE,
    N_G               = N_seg_HE,
    N_F               = N_seg_HE,
    SSInit            = SSInit,
    gasQuasiStatic    = true,
    fluidQuasiStatic  = true,
    metalWall(
      WallRes=false),
    table_k_metalwall = table_k_metalwall
  )
  annotation(
    Placement(transformation(extent = {{-20, -20}, {20, 20}}, rotation = 0)));

  // variable for validation
  Modelica.SIunits.Power Q_out = (HE.gasIn.h_outflow - HE.gasOut.h_outflow) * HE.gasIn.m_flow; 
  Modelica.SIunits.Power Q_in = (HE.waterOut.h_outflow - HE.waterIn.h_outflow) * HE.waterIn.m_flow;
  Boolean isQMatch = abs(Q_out - Q_in) < 1e-3;   
  
equation
 
  connect(source_cold.flange, r_HE_cin.inlet);
  connect(r_HE_cin.outlet, HE.waterIn) annotation(
    Line(points = {{-1.83697e-015, 50}, {-1.83697e-015, 20}, {0, 20}}, color = {0, 0, 255}, thickness = 0.5, smooth = Smooth.None));
  connect(HE.waterOut, r_HE_cout.inlet) annotation(
    Line(points = {{8.88178e-016, -44}, {8.88178e-016, -20}, {0, -20}}, thickness = 0.5, color = {0, 0, 255}));      
  connect(r_HE_cout.outlet, sink_cold.flange) annotation(
    Line(points = {{1.83697e-015, -70}, {1.83697e-015, -56}, {-8.88178e-016, -56}}, thickness = 0.5, color = {0, 0, 255}));
  
  connect(source_hot.flange, r_HE_hin.inlet);
  connect(r_HE_hin.outlet, HE.gasIn) annotation(
        Line(points = {{-50, 0}, {-20, 0}}, color = {159, 159, 223}, thickness = 0.5, smooth = Smooth.None)); 
  connect(HE.gasOut, r_HE_hout.inlet) annotation(
    Line(points = {{34, 0}, {34, 0}, {20, 0}}, color = {159, 159, 223}, thickness = 0.5));
  connect(r_HE_hout.outlet, sink_hot.flange) annotation(
    Line(points = {{46, 0}, {46, 0}, {60, 0}}, color = {159, 159, 223}, thickness = 0.5));      
    
annotation(
    Diagram(graphics),
    // for steady-state simulation - value check
    experiment(StartTime = 0, StopTime = 2, Tolerance = 1e-3, Interval = 1),
    // for complete transient simulation
    // experiment(StartTime = 0, StopTime = 600, Tolerance = 1e-3, Interval = 10),
    __OpenModelica_commandLineOptions = "--matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection -d=initialization,NLSanalyticJacobian,aliasConflicts",        
    // remove the option flag --matchingAlgorithm=PFPlusExt, which may lead to 'Internal error - IndexReduction.dynamicStateSelectionWork failed!' during Translation
    // __OpenModelica_commandLineOptions = "--matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection -d=initialization,NLSanalyticJacobian,aliasConflicts,bltdump",    
    __OpenModelica_simulationFlags(lv = "LOG_DEBUG,LOG_NLS,LOG_NLS_V,LOG_STATS,LOG_INIT,LOG_STDOUT, -w", outputFormat = "mat", s = "dassl", nls = "homotopy")
    );
end TestTP_PCHE_PHE;
