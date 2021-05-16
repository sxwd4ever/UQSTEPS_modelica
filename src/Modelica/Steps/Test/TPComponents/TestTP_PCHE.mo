within Steps.Test.TPComponents;

model TestTP_PCHE
  "Test for ThermoPower based PCHE model"  
    
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
  
  // package medium_hot = Steps.Media.CO2;
  // package medium_cold = Steps.Media.CO2;
  // package medium_hot = ExternalMedia.Examples.CO2CoolProp;
  // package medium_cold = ExternalMedia.Examples.CO2CoolProp;    
  package medium_hot = Steps.Media.SCO2;
  package medium_cold = Steps.Media.SCO2; 
  // package medium_hot = Modelica.Media.IdealGases.SingleGases.CO2;
  // package medium_cold = Modelica.Media.IdealGases.SingleGases.CO2;    
  parameter Real table_k_metalwall[:,:] = [293.15, 12.1; 373.15, 16.3; 773.15, 21.5];    
  
  parameter Real Cf_C1 = 1.626, Cf_C2 = 1, Cf_C3 = 1;
  // parameter Real Cf_C1_cold = 1, Cf_C2_cold = 1, Cf_C3_cold = 1;
  parameter Real use_rho_bar = -1.0;  
  parameter Real rho_bar_hot = 1.0;
  parameter Real rho_bar_cold = 1.0;    
  
  parameter Model.RCBCycleConfig cfg(
    N_ch_HTR = 30000,
    L_HTR    = 2.5, 
    r_i_HTR  = 1.5e-3,   
    r_o_HTR  = 1.5e-3,    
    N_ch_LTR = 30000,
    L_LTR    = 2.5,
    r_i_LTR  = 1.5e-3,   
    r_o_LTR  = 1.5e-3         
  );
  /*
  (
    mdot_heater = 90,
    table_k_LTR_wall = table_k_metalwall,
    table_k_HTR_wall = table_k_metalwall
  );
  */
  
  // set the values of parameters accordingly - For HTR test.
  // modify it to cfg_HE = cfg.cfg_LTR for LTR test
  parameter Model.HeatExchangerConfig cfg_HE     = cfg.cfg_HTR;
  parameter Model.ThermoState st_source_hot      = cfg_HE.cfg_hot.st_in;
  parameter Model.ThermoState st_sink_hot        = cfg_HE.cfg_hot.st_out;
  parameter Model.ThermoState st_source_cold     = cfg_HE.cfg_cold.st_in;
  parameter Model.ThermoState st_sink_cold       = cfg_HE.cfg_cold.st_out;
  parameter Integer N_seg_HE                     = cfg.cfg_HTR.cfg_hot.geo_path.N_seg;

  //Components
  inner ThermoPower.System system(allowFlowReversal = false, initOpt=ThermoPower.Choices.Init.Options.fixedState) annotation(
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
    /*
    // with Kim Correlation
    redeclare replaceable model HeatTransfer_F = Steps.TPComponents.KimPCHEHeatTransferFV(),    
    // ThermoPower.Thermal.HeatTransferFV.ConstantHeatTransferCoefficient(gamma = thermo_LTR_cold.gamma_he),
    // redeclare replaceable model HeatTransfer_G = ThermoPower.Thermal.HeatTransferFV.IdealHeatTransfer, 
    redeclare replaceable model HeatTransfer_G = Steps.TPComponents.KimPCHEHeatTransferFV(),
    fluidFlow(
      heatTransfer(
        pitch = cfg_HE.cfg_cold.l_pitch,
        phi   = cfg_HE.cfg_cold.a_phi
      )
    ),
    gasFlow(
      heatTransfer(
        pitch = cfg_HE.cfg_hot.l_pitch,
        phi   = cfg_HE.cfg_hot.a_phi
      )
    ),
    */
    
    // with Marchionni Correlation    
    redeclare replaceable model HeatTransfer_F = Steps.TPComponents.MarchionniPCHEHeatTransferFV(),
    redeclare replaceable model HeatTransfer_G = Steps.TPComponents.MarchionniPCHEHeatTransferFV(),
    gasFlow(
      heatTransfer(
        pitch       = cfg_HE.cfg_hot.l_pitch,
        phi         = cfg_HE.cfg_hot.a_phi,         
        Cf_C1       = Cf_C1, 
        Cf_C2       = Cf_C2, 
        Cf_C3       = Cf_C3, 
        use_rho_bar = use_rho_bar,
        rho_bar     = rho_bar_hot,
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
        useAverageTemperature = false)
    ),  
    redeclare model HeatExchangerTopology = ThermoPower.Thermal.HeatExchangerTopologies.CounterCurrentFlow, 
    cfg               = cfg_HE,
    N_G               = N_seg_HE,
    N_F               = N_seg_HE,
    SSInit            = SSInit,
    gasQuasiStatic    = true,
    fluidQuasiStatic  = true
    //metalWall(WallRes=false) 
    //table_k_metalwall = table_k_metalwall
    // override the values of Am and L of metaltubeFV
    // to make them agree with semi-circular tube of PCHE
    // ('final' modifier of Am in metalTubeFv was removed as well)
    //metalTube(WallRes=false, L = 1, rhomcm=200, Am = HE.metalVol / 1) 
  )
  annotation(
    Placement(transformation(extent = {{-20, -20}, {20, 20}}, rotation = 0)));

  // variable for validation
  Modelica.SIunits.Power Q_out    = (HE.gasIn.h_outflow - HE.gasOut.h_outflow) * HE.gasIn.m_flow;
  Modelica.SIunits.Power Q_in     = (HE.waterOut.h_outflow - HE.waterIn.h_outflow) * HE.waterIn.m_flow;
  Boolean                isQMatch = abs(Q_out -Q_in) < 1e-3;
equation

/*
  // HE alone
  connect(source_cold.flange, HE.waterIn) annotation(
    Line(points = {{-1.83697e-015, 50}, {-1.83697e-015, 20}, {0, 20}}, color = {0, 0, 255}, thickness = 0.5, smooth = Smooth.None));
  connect(HE.waterOut, sink_cold.flange) annotation(
    Line(points = {{1.83697e-015, -70}, {1.83697e-015, -56}, {-8.88178e-016, -56}}, thickness = 0.5, color = {0, 0, 255}));
  
  connect(source_hot.flange, HE.gasIn) annotation(
        Line(points = {{-50, 0}, {-20, 0}}, color = {159, 159, 223}, thickness = 0.5, smooth = Smooth.None)); 
  connect(HE.gasOut, sink_hot.flange) annotation(
    Line(points = {{46, 0}, {46, 0}, {60, 0}}, color = {159, 159, 223}, thickness = 0.5));    
*/   
  
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
    experiment(StartTime = 0, StopTime = 10, Tolerance = 1e-3, Interval = 2),
    __OpenModelica_commandLineOptions = "--matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection -d=initialization,NLSanalyticJacobian,aliasConflicts",    
    __OpenModelica_simulationFlags(lv = "LOG_DEBUG,LOG_NLS,LOG_NLS_V,LOG_STATS,LOG_INIT,LOG_STDOUT, -w", outputFormat = "mat", s = "dassl", nls = "homotopy"));
end TestTP_PCHE;
