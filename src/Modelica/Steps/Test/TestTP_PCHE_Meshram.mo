within Steps.Test;

model TestTP_PCHE_Meshram
  "Test for ThermoPower based PCHE model against Meshram [2016]"  
    
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

  package medium_hot = Steps.Media.SCO2;
  package medium_cold = Steps.Media.SCO2; 
  
  // geometry parameters
  constant Real pi = Modelica.Constants.pi;
  parameter Integer N_ch = integer(1e4) "channel number";
  parameter Integer N_seg = 10 "number of segments in one tube";
  parameter SI.Length D_ch = 2e-3 "channel diameter, semi circular tube";
  parameter SI.Length r_ch = D_ch / 2 "channel radiaus";
  parameter SI.Length L_fp = 200e-3 "channel flow path length";  
  parameter SI.Length L_pitch = 12e-3 "pitch length";
  parameter Real a_phi = 36 "pitch angle degree";
  parameter SI.Length H_ch = 3.2e-3 "Height of the solid domain, containing one cold tube and one hot tube";
  parameter SI.Length W_ch = 2.5e-3 "Width of the solid domain";
  parameter SI.Area A = pi * r_ch ^2 / 2 "Area of cross section of semi circular tube";

  // boundary conditon
  
  // zigzag higher T
  parameter SI.Velocity u_hot_in = 7.564 "hot inlet velocity m/s";
  parameter SI.Velocity u_cold_in = 1.876 "cold inlet velocity m/s";
  parameter SI.Pressure p_hot_in =  from_bar(90) "hot inlet pressure";
  parameter SI.Pressure p_cold_in = from_bar(225) "cold inlet pressure";
  parameter SI.Temperature T_hot_in = 730 "hot inlet temperature, K";
  parameter SI.Temperature T_hot_out = 576.69 "cold outlet temperature, K";
  parameter SI.Temperature T_cold_in = 500 "cold inlet temperature, K";
  parameter SI.Temperature T_cold_out = 639.15 "cold outlet temperature, K";
  
  // pressure drop correction coefficient 
  parameter Real kc_dp = 1.0;
  parameter Real C1 = 1.0;
  parameter Real C2 = 1.0;
  
  // meshram's cp and rho for alloy Inconel 617
  parameter Modelica.SIunits.Density rho_wall = 8360 "density of wall, kg/m3";
  parameter Modelica.SIunits.SpecificHeatCapacity cp_wall = 417 "cp of wall, J/kg-K";  

/*  
  // zigzag lower T
  parameter SI.Velocity u_hot_in = 1.345 "hot inlet velocity m/s";
  parameter SI.Velocity u_cold_in = 0.806 "cold inlet velocity m/s";
  parameter SI.Pressure p_hot_in =  from_bar(90) "hot inlet pressure";
  parameter SI.Pressure p_cold_in = from_bar(225) "cold inlet pressure";
  parameter SI.Temperature T_hot_in = 630 "hot inlet temperature, K";
  parameter SI.Temperature T_hot_out = 466.69 "cold outlet temperature, K";
  parameter SI.Temperature T_cold_in = 400 "cold inlet temperature, K";
  parameter SI.Temperature T_cold_out = 522.23 "cold outlet temperature, K";  
*/  
  parameter SI.Density rho_hot_in = medium_hot.density_pT(p_hot_in, T_hot_in);
  parameter SI.Density rho_cold_in = medium_cold.density_pT(p_cold_in, T_cold_in);
  parameter SI.MassFlowRate mdot_hot_in = rho_hot_in * A * u_hot_in * N_ch;
  parameter SI.MassFlowRate mdot_cold_in = rho_cold_in * A * u_cold_in * N_ch;

  // use configuration of LTR for this test since the mdot are different for hot and cold side
  parameter Model.PBConfig_PCHE cfg(
    p_pump_in = p_hot_in,
    p_pump_out = p_cold_in,
    mdot_main = mdot_hot_in,
    mdot_pump = mdot_cold_in, 
    T_LTR_cold_in = T_cold_in, 
    T_LTR_cold_out = T_cold_out,
    T_HTR_hot_out = T_hot_in, // T_LTR_hot_in = T_HTR_hot_out,
    T_LTR_hot_out = T_hot_out,
    r_LTR = r_ch,
    L_LTR = L_fp,
    N_ch_LTR = N_ch,
    N_seg = N_seg,
    pitch = L_pitch,
    phi = a_phi,
    rho_wall = rho_wall,
    cp_wall = cp_wall
  );
  
  // set the values of parameters accordingly
  parameter HEBoundaryCondition bc_HE = cfg.bc_LTR; 

  //Components
  inner ThermoPower.System system(allowFlowReversal = false, initOpt=ThermoPower.Choices.Init.Options.noInit) annotation(
    Placement(transformation(extent = {{80, 80}, {100, 100}})));  
  
  ThermoPower.Gas.SourceMassFlow source_cold(
    redeclare package Medium = medium_cold, 
    T = bc_HE.st_cold_in.T, 
    p0 = bc_HE.st_cold_in.p, 
    //use_in_T = false, 
    w0 = bc_HE.st_cold_in.mdot,
    gas(p(nominal = bc_HE.st_cold_in.p), 
    T(nominal=bc_HE.st_cold_in.T))) 
  annotation(
    Placement(transformation(origin = {0, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 270)));
  
  ThermoPower.Gas.SinkPressure sink_cold(
    redeclare package Medium = medium_cold, 
    p0 = bc_HE.st_cold_out.p, 
    T = bc_HE.st_cold_out.T) 
  annotation(
    Placement(transformation(origin = {0, -80}, extent = {{-10, -10}, {10, 10}}, rotation = 270)));
 
  ThermoPower.Gas.SinkPressure sink_hot(
    redeclare package Medium = medium_hot,
    T = bc_HE.st_hot_out.T, 
    p0 = bc_HE.st_hot_out.p) 
  annotation(
    Placement(transformation(extent = {{60, -10}, {80, 10}}, rotation = 0)));
  
  ThermoPower.Gas.SourceMassFlow source_hot(
    redeclare package Medium = medium_hot, 
    T = bc_HE.st_hot_in.T, 
    p0 = bc_HE.st_hot_in.p, 
    w0 = bc_HE.st_hot_in.mdot,
    //use_in_T = false,
    gas(p(nominal = bc_HE.st_hot_in.p), 
    T(nominal=bc_HE.st_hot_in.T))) 
  annotation(
    Placement(transformation(extent = {{-70, -10}, {-50, 10}}, rotation = 0))); 
  
  ThermoPower.Gas.SensT T_waterIn(redeclare package Medium = medium_cold) annotation(
    Placement(transformation(origin = {4, -50}, extent = {{-10, -10}, {10, 10}}, rotation = 270)));
  ThermoPower.Gas.SensT T_waterOut(redeclare package Medium = medium_cold) annotation(
    Placement(transformation(origin = {4, -50}, extent = {{-10, -10}, {10, 10}}, rotation = 270)));
  
  ThermoPower.Gas.SensT T_gasIn(redeclare package Medium = medium_hot);
  ThermoPower.Gas.SensT T_gasOut(redeclare package Medium = medium_hot);
    
  TPComponents.PCHE HE(
    redeclare package FluidMedium = medium_cold, 
    redeclare package FlueGasMedium = medium_hot,     
    redeclare replaceable model HeatTransfer_F = TPComponents.MarchionniPCHEHeatTransferFV(),
    //redeclare replaceable model HeatTransfer_F = TPComponents.KimPCHEHeatTransferFV(), 
    // ThermoPower.Thermal.HeatTransferFV.ConstantHeatTransferCoefficient(gamma = thermo_LTR_cold.gamma_he),
    // redeclare replaceable model HeatTransfer_G = ThermoPower.Thermal.HeatTransferFV.IdealHeatTransfer, 
    redeclare replaceable model HeatTransfer_G = TPComponents.MarchionniPCHEHeatTransferFV(),
    // redeclare replaceable model HeatTransfer_G = ThermoPower.Thermal.HeatTransferFV.IdealHeatTransfer, 
    //redeclare replaceable model HeatTransfer_G = Steps.TPComponents.KimPCHEHeatTransferFV(),
    redeclare model HeatExchangerTopology = ThermoPower.Thermal.HeatExchangerTopologies.CounterCurrentFlow, 
    bc = bc_HE, 
    geo_hot = cfg.cfg_LTR_hot.geo,
    geo_cold = cfg.cfg_LTR_cold.geo,
    geo_tube = cfg.cfg_LTR_tube.geo,  
    thermo_hot = cfg.cfg_LTR_hot.thermo,
    thermo_cold = cfg.cfg_LTR_cold.thermo,
    thermo_tube = cfg.cfg_LTR_tube.thermo, 
    L = L_fp,
    SSInit = true,
    gasQuasiStatic = true,
    fluidQuasiStatic = true,
    gasFlow(heatTransfer(pitch = cfg.pitch, phi = cfg.phi, kc_dp = kc_dp, C1 = C1, C2 = C2)),
    fluidFlow(heatTransfer(pitch = cfg.pitch, phi = cfg.phi, kc_dp = kc_dp, C1 = C1, C2 = C2))
    // override the values of Am and L of metaltubeFV
    // to make them agree with semi-circular tube of PCHE
    // ('final' modifier of Am in metalTubeFv was removed as well)
    //metalTube(WallRes=false, L = 1, rhomcm=200, Am = HE.metalVol / 1) 
  )
  annotation(
    Placement(transformation(extent = {{-20, -20}, {20, 20}}, rotation = 0)));

  // variable for validation
  Modelica.SIunits.Power Q_out = (HE.gasIn.h_outflow - HE.gasOut.h_outflow) * HE.gasIn.m_flow; 
  Modelica.SIunits.Power Q_in = (HE.waterOut.h_outflow - HE.waterIn.h_outflow) * HE.waterIn.m_flow;
  Boolean isQMatch = abs(Q_out -Q_in) < 1e-3;  
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
  
  connect(source_cold.flange, T_waterIn.inlet);
  connect(T_waterIn.outlet, HE.waterIn) annotation(
    Line(points = {{-1.83697e-015, 50}, {-1.83697e-015, 20}, {0, 20}}, color = {0, 0, 255}, thickness = 0.5, smooth = Smooth.None));
  connect(HE.waterOut, T_waterOut.inlet) annotation(
    Line(points = {{8.88178e-016, -44}, {8.88178e-016, -20}, {0, -20}}, thickness = 0.5, color = {0, 0, 255}));      
  connect(T_waterOut.outlet, sink_cold.flange) annotation(
    Line(points = {{1.83697e-015, -70}, {1.83697e-015, -56}, {-8.88178e-016, -56}}, thickness = 0.5, color = {0, 0, 255}));
  
  connect(source_hot.flange, T_gasIn.inlet);
  connect(T_gasIn.outlet, HE.gasIn) annotation(
        Line(points = {{-50, 0}, {-20, 0}}, color = {159, 159, 223}, thickness = 0.5, smooth = Smooth.None)); 
  connect(HE.gasOut, T_gasOut.inlet) annotation(
    Line(points = {{34, 0}, {34, 0}, {20, 0}}, color = {159, 159, 223}, thickness = 0.5));
  connect(T_gasOut.outlet, sink_hot.flange) annotation(
    Line(points = {{46, 0}, {46, 0}, {60, 0}}, color = {159, 159, 223}, thickness = 0.5));    
  
annotation(
    Diagram(graphics),
    experiment(StartTime = 0, StopTime = 10, Tolerance = 1e-3, Interval = 2),
    __OpenModelica_commandLineOptions = "--matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection -d=initialization,NLSanalyticJacobian,aliasConflicts,bltdump",    
    __OpenModelica_simulationFlags(lv = "LOG_DEBUG,LOG_NLS,LOG_NLS_V,LOG_STATS,LOG_INIT,LOG_STDOUT, -w", outputFormat = "mat", s = "dassl", nls = "homotopy"));
end TestTP_PCHE_Meshram;
