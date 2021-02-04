within Steps.Test;

model TestTP_PCHE_Marchionni
  "Test for ThermoPower based PCHE model against Marchionni [2019]"  
    
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
  parameter Integer N_ch = 10000 "channel number";
  parameter Integer N_seg = 20 "number of segments in one tube";
  parameter SI.Length D_ch = 2e-3 "channel diameter, semi circular tube";
  parameter SI.Length r_ch = D_ch / 2 "channel radiaus";
  parameter SI.Length L_fp = 272e-3 "channel flow path length";  
  parameter SI.Length L_pitch = 12.3e-3 "pitch length"; 
  parameter Real a_phi "pitch angle degree";
  parameter SI.Length H_ch = 3.26e-3 "Height of the solid domain, containing one cold tube and one hot tube";
  parameter SI.Length W_ch = 1.27e-3 * 2"Width of the solid domain";
  parameter SI.Area A = pi * r_ch ^2 / 2 "Area of cross section of semi circular tube";

  // boundary conditon
  
  // zigzag higher T
  parameter Real G_hot_in = 509.3 "hot inlet mass flux kg/(m^2 s";
  parameter Real G_cold_in = 509.3 "cold inlet mass flux kg/(m^2 s";
  parameter SI.Pressure p_hot_in =  from_bar(75) "hot inlet pressure";
  parameter SI.Pressure p_cold_in = from_bar(150) "cold inlet pressure";
  parameter SI.Temperature T_hot_in = from_degC(400) "hot inlet temperature, K";
  parameter SI.Temperature T_hot_out = from_degC(140) "cold outlet temperature, K";
  parameter SI.Temperature T_cold_in = from_degC(100) "cold inlet temperature, K";
  parameter SI.Temperature T_cold_out = from_degC(300) "cold outlet temperature, K";
  
  // pressure drop correction coefficient 
  // parameter Real kc_dp = 1.0;
  
  parameter Real Cf_C1 = 1, Cf_C2 = 1, Cf_C3 = 1;
  parameter Real use_rho_bar;  
  parameter Real rho_bar_hot;
  parameter Real rho_bar_cold;
  
  // meshram's cp and rho for alloy Inconel 617
  parameter Modelica.SIunits.Density rho_wall = 7990 "density of wall, kg/m3";
  parameter Modelica.SIunits.SpecificHeatCapacity cp_wall = 500 "cp of wall, J/kg-K";  

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
  parameter SI.MassFlowRate mdot_hot_in = G_hot_in * A * N_ch;
  parameter SI.MassFlowRate mdot_cold_in = G_cold_in * A * N_ch;

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
    Placement(visible = true, transformation(origin = {0, 80}, extent = {{-10, -10}, {10, 10}}, rotation = 270)));
   
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
    Placement(visible = true, transformation(extent = {{70, -10}, {90, 10}}, rotation = 0)));
  
  ThermoPower.Gas.SourceMassFlow source_hot(
    redeclare package Medium = medium_hot, 
    T = bc_HE.st_hot_in.T, 
    p0 = bc_HE.st_hot_in.p, 
    w0 = bc_HE.st_hot_in.mdot,
    //use_in_T = false,
    gas(p(nominal = bc_HE.st_hot_in.p), 
    T(nominal=bc_HE.st_hot_in.T))) 
  annotation(
    Placement(visible = true, transformation(extent = {{-92, -10}, {-72, 10}}, rotation = 0)));
  
  TPComponents.GasStateReader sr_water_in(redeclare package Medium = medium_cold) annotation(
    Placement(visible = true, transformation(origin = {0, 48}, extent = {{-10, -10}, {10, 10}}, rotation = 270)));
 
  TPComponents.GasStateReader sr_water_out(redeclare package Medium = medium_cold) annotation(
    Placement(visible = true, transformation(origin = {0, -48}, extent = {{-10, -10}, {10, 10}}, rotation = 270)));
 
  TPComponents.GasStateReader sr_gas_in(redeclare package Medium = medium_hot) annotation(
    Placement(visible = true, transformation(origin = {-50, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
 
  TPComponents.GasStateReader sr_gas_out(redeclare package Medium = medium_hot) annotation(
    Placement(visible = true, transformation(origin = {44, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
     
  TPComponents.PCHE HE(
    redeclare package FluidMedium = medium_cold, 
    redeclare package FlueGasMedium = medium_hot,     
    // use Marchionni PCHE HeatTransfer
    redeclare replaceable model HeatTransfer_F = TPComponents.MarchionniPCHEHeatTransferFV(),
    redeclare replaceable model HeatTransfer_G = TPComponents.MarchionniPCHEHeatTransferFV(),
    gasFlow(heatTransfer(
      pitch = cfg.pitch, 
      phi = cfg.phi, 
      // kc_dp = kc_dp, 
      Cf_C1 = Cf_C1, 
      Cf_C2 = Cf_C2, 
      Cf_C3 = Cf_C3, 
      use_rho_bar = use_rho_bar,
      rho_bar = rho_bar_hot)),
    fluidFlow(heatTransfer(
      pitch = cfg.pitch, 
      phi = cfg.phi, 
      //kc_dp = kc_dp, 
      Cf_C1 = Cf_C1, 
      Cf_C2 = Cf_C2, 
      Cf_C3 = Cf_C3, 
      use_rho_bar = use_rho_bar, 
      rho_bar = rho_bar_cold)),    
    /*    
    // use Kim PCHE HeatTransfer
    redeclare replaceable model HeatTransfer_F = TPComponents.KimPCHEHeatTransferFV(), 
    redeclare replaceable model HeatTransfer_G = Steps.TPComponents.KimPCHEHeatTransferFV(),    
    gasFlow(heatTransfer(pitch = cfg.pitch, phi = cfg.phi, kc_dp = kc_dp)),
    fluidFlow(heatTransfer(pitch = cfg.pitch, phi = cfg.phi, kc_dp = kc_dp)), 
    */
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
    fluidQuasiStatic = true
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
  
  Real T_hot_out_act = sr_gas_out.T;
  Real T_cold_out_act = sr_water_out.T;
  Real dp_hot_act = sum(HE.gasFlow.heatTransfer.dp);
  Real dp_cold_act = sum(HE.fluidFlow.heatTransfer.dp);
  // Real dp_hot_act_m = sum(HE.gasFlow.heatTransfer.dp_m) * 10 "actual dp calculated by Eq. 1 [Marchionni 2019]"; 
  // Real dp_cold_act_m = sum(HE.fluidFlow.heatTransfer.dp_m) * 10 "actual dp calculated by Eq. 1 [Marchionni 2019]"; 
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
  
  connect(source_cold.flange, sr_water_in.inlet) annotation(
    Line(points = {{0, 54}, {0, 54}, {0, 70}, {0, 70}}, color = {0, 0, 255}, thickness = 0.5));
  connect(sr_water_in.outlet, HE.waterIn) annotation(
    Line(points = {{0, 42}, {0, 42}, {0, 20}, {0, 20}}, color = {0, 0, 255}, thickness = 0.5));
  connect(HE.waterOut, sr_water_out.inlet) annotation(
    Line(points = {{0, -42}, {0, -42}, {0, -20}, {0, -20}}, color = {0, 0, 255}, thickness = 0.5));
  connect(sr_water_out.outlet, sink_cold.flange) annotation(
    Line(points = {{0, -54}, {0, -54}, {0, -70}, {0, -70}}, color = {0, 0, 255}, thickness = 0.5));

  connect(source_hot.flange, sr_gas_in.inlet) annotation(
    Line(points = {{-72, 0}, {-56, 0}, {-56, 0}, {-56, 0}}, color = {159, 159, 223}));
  connect(sr_gas_in.outlet, HE.gasIn) annotation(
    Line(points = {{-44, 0}, {-22, 0}, {-22, 0}, {-20, 0}}, color = {159, 159, 223}));
  connect(HE.gasOut, sr_gas_out.inlet) annotation(
    Line(points = {{20, 0}, {38, 0}, {38, 0}, {38, 0}}, color = {159, 159, 223}));
  connect(sr_gas_out.outlet, sink_hot.flange) annotation(
    Line(points = {{50, 0}, {70, 0}, {70, 0}, {70, 0}}, color = {159, 159, 223}));
 
annotation(
    Diagram(graphics),
    experiment(StartTime = 0, StopTime = 10, Tolerance = 1e-3, Interval = 2),
    __OpenModelica_commandLineOptions = "--matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection -d=initialization,NLSanalyticJacobian,aliasConflicts,bltdump",    
    __OpenModelica_simulationFlags(lv = "LOG_DEBUG,LOG_NLS,LOG_NLS_V,LOG_STATS,LOG_INIT,LOG_STDOUT, -w", outputFormat = "mat", s = "dassl", nls = "homotopy"));
end TestTP_PCHE_Marchionni;
