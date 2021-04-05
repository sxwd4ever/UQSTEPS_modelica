within Steps.Test.TPComponents.Transient;

model TestDyn_Comps_PCHE_ramp
  "Test for ThermoPower based components (heater + turbine + HTR + LTR, mainly PCHE) with ramp change in mass flow rate"  
  import Modelica.SIunits.Conversions.{from_degC, from_deg};
  import Modelica.SIunits.{Temperature, Pressure, SpecificEnthalpy};
  import Util = Utilities.Util;
  import Steps.Utilities.CoolProp.PropsSI;  
  import Steps.Components.{PCHEGeoParam};
  import Steps.Model.{PBConfiguration, SimParam, EntityConfig, EntityGeoParam, EntityThermoParam, ThermoState, HEBoundaryCondition} ;
  import Model.PBConfiguration;
  import ThermoPower.Choices.Init.Options;
  import ThermoPower.System;
  import ThermoPower.Gas;
  
  // package medium_hot = Modelica.Media.IdealGases.SingleGases.CO2; //Steps.Media.CO2;
  // package medium_cold = Modelica.Media.IdealGases.SingleGases.CO2; //Steps.Media.CO2;
  package medium_hot = Steps.Media.SCO2;
  package medium_cold = Steps.Media.SCO2;  
  
  // package medium_hot = ExternalMedia.Examples.CO2CoolProp;
  // package medium_cold = ExternalMedia.Examples.CO2CoolProp;
  
  package medium_heater = Steps.Media.ThermiaOilD;
  
  //package medium_heater = ThermoPower.Water.StandardWater;// Modelica.Media.IdealGases.SingleGases.CO2;

  parameter Real Cf_C1 = 1, Cf_C2 = 1, Cf_C3 = 1;
  
  // geometry parameters
  constant Real pi = Modelica.Constants.pi;
  parameter Integer N_ch = integer(2400) "channel number";
  parameter Integer N_seg = 10 "number of segments in one tube";
  parameter SI.Length D_ch = 1.72e-3  "channel diameter, semi circular tube";
  parameter SI.Length r_ch = D_ch / 2 "channel radiaus";
  parameter SI.Length L_fp = 270 * 1e-3 "equivalent valid channel flow path length"; 
  parameter SI.Length L_pitch = 12e-3 "pitch length";
  parameter Real a_phi = 36 "pitch angle °";
  parameter SI.Length H_ch = 4.17e-3 "Height of the solid domain, containing one cold tube and one hot tube";
  parameter SI.Length W_ch = 2.3e-3 "Width of the solid domain";
  parameter SI.Length T_wall = 0.51e-3 "Wall thinckness";
  parameter SI.Length L_wall = 420e-3 "Length of wall, not necessarily equals to length of flow path";
  parameter SI.Area A = pi * r_ch ^2 / 2 "Area of cross section of semi circular tube";
    
  
  // Stainless 316, 316L, 317, 317L
  parameter Modelica.SIunits.Density rho_wall = 8030 "density of wall, kg/m3";
  parameter Modelica.SIunits.SpecificHeatCapacity cp_wall = 485 "cp of wall, J/kg-K";  
  // thermal conductivity (T in K) https://www.theworldmaterial.com/aisi-316-ss316-stainless-steel-properties-composition/
  // parameter Real table_k_metalwall[:,:] = [20, 12.1; 100, 16.3; 500, 21.5];
  parameter Real table_k_metalwall[:,:] = [293.15, 12.1; 373.15, 16.3; 773.15, 21.5];  

  // input parameters of the power block
  parameter Modelica.SIunits.MassFlowRate mdot_main = 125 "kg/s, mass flow in the main path of PB, which follows the power demand";
  parameter Modelica.SIunits.MassFlowRate mdot_heater_hot = 90 "kg/s, mass flow rate of heater's hot fluid";
  parameter Real gamma = 0.4 "split ratio, mdot_bypass/mdot_main";
  
  parameter Modelica.SIunits.Temperature T_heater_hot = from_degC(800) "K, Temperature of heater's hot fluid";  
  parameter Modelica.SIunits.Temperature T_cooler_cold = from_degC(45) "K, Temperature of cooler's cold fluid";
 
  // results based on sscar's simulation for 10 Mw power block
  parameter Model.PBConfig_PCHE cfg(
    mdot_main = mdot_main, 
    mdot_pump = mdot_main * (1 - gamma),
    // mdot_pump = mdot_main,
    N_ch_LTR = N_ch,
    N_ch_HTR = N_ch,
    // N_seg = N_seg,
    mdot_heater = mdot_heater_hot,
    T_heater_hot_in = T_heater_hot,
    T_heater_hot_out = T_heater_hot - 200,
    T_heater_cold_out = T_heater_hot - 90,
    T_cooler_cold_in = T_cooler_cold,
    pitch = L_pitch
  );
 
  // set the values of parameters accordingly
  parameter HEBoundaryCondition bc_HTR = cfg.bc_HTR;  
  parameter HEBoundaryCondition bc_LTR = cfg.bc_LTR;
  parameter HEBoundaryCondition bc_heater = cfg.bc_heater;
  
  parameter ThermoState st_bypass = cfg.st_bypass;
  parameter EntityThermoParam thermo_mixer = cfg.cfg_mixer.thermo;
  
  // load design parameters such as geometry values.
  parameter EntityGeoParam geo_HTR_hot = cfg.cfg_HTR_hot.geo;
  parameter EntityGeoParam geo_HTR_cold = cfg.cfg_HTR_cold.geo;
  
  parameter EntityGeoParam geo_LTR_hot = cfg.cfg_LTR_hot.geo;
  parameter EntityGeoParam geo_LTR_cold = cfg.cfg_LTR_cold.geo; 

  // load thermodynamic parameters. 
  parameter EntityThermoParam thermo_HTR_hot = cfg.cfg_HTR_hot.thermo;
  parameter EntityThermoParam thermo_HTR_cold = cfg.cfg_HTR_cold.thermo;
  
  parameter EntityThermoParam thermo_LTR_hot = cfg.cfg_LTR_hot.thermo;
  parameter EntityThermoParam thermo_LTR_cold = cfg.cfg_LTR_cold.thermo;

  //Components
  inner ThermoPower.System system(allowFlowReversal = false, initOpt=ThermoPower.Choices.Init.Options.steadyState) annotation(
    Placement(transformation(extent = {{80, 80}, {100, 100}})));
  
  // global init opition (system.initOpt) leads to order reduction error
  // use this flag to control the initialization of all components instead. 
  // parameter Boolean SSInit = false "Steady-state initialization";

  ThermoPower.Gas.SourceMassFlow source_cold(
    redeclare package Medium = medium_cold, 
    T = bc_LTR.st_cold_in.T, 
    p0 = bc_LTR.st_cold_in.p, 
    use_in_T = false, 
    w0 = bc_LTR.st_cold_in.mdot) 
  annotation(
    Placement(transformation(origin = {0, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 270)));
  

  ThermoPower.Gas.SourceMassFlow source_mixer_in(
    redeclare package Medium = medium_cold,
    T = st_bypass.T,
    p0 = st_bypass.p,
    use_in_T = false,
    w0 = st_bypass.mdot    
  );  
  
  ThermoPower.Gas.SinkPressure sink_cold(
    redeclare package Medium = medium_cold, 
    // p0 = bc_HTR.st_cold_out.p, 
    // p0 = bc_heater.st_cold_out.p,
    p0 = bc_HTR.st_cold_out.p,
    // T = bc_HTR.st_cold_out.T)
    T = bc_HTR.st_cold_out.T)
  annotation(
    Placement(transformation(origin = {0, -80}, extent = {{-10, -10}, {10, 10}}, rotation = 270)));

  ThermoPower.Gas.SinkPressure sink_hot(
    redeclare package Medium = medium_hot,
    T = bc_LTR.st_hot_out.T, 
    p0 = bc_LTR.st_hot_out.p) 
  annotation(
    Placement(transformation(extent = {{60, -10}, {80, 10}}, rotation = 0)));

  ThermoPower.Gas.SourceMassFlow source_hot(
    redeclare package Medium = medium_hot, 
    // T = from_degC(730.32), 
    //T = bc_heater.st_cold_out.T, 
    T = bc_HTR.st_hot_in.T, 
    //p0 = bc_heater.st_cold_out.p, 
    p0 = bc_HTR.st_hot_in.p, 
    //w0 = bc_heater.st_cold_out.mdot,
    w0 = bc_HTR.st_hot_in.mdot,
    use_in_T = false) 
  annotation(
    Placement(transformation(extent = {{-70, -10}, {-50, 10}}, rotation = 0))); 

  //gas
  parameter Modelica.SIunits.MassFlowRate gasNomFlowRate = 125 "Nominal mass flowrate";
  parameter Modelica.SIunits.Pressure gasNomPressure = 9e6 "Nominal pressure in the gas side inlet";
  parameter Modelica.SIunits.Temperature Tstart_G_In = 883 "Inlet gas temperature start value";
  parameter Modelica.SIunits.Temperature Tstart_G_Out = 643 "Outlet gas temperature start value";
  //fluid
  parameter Modelica.SIunits.MassFlowRate fluidNomFlowRate = 125 "Nominal flow rate through the fluid side";
  parameter Modelica.SIunits.Pressure fluidNomPressure = 9e6 "Nominal pressure in the fluid side inlet";
  parameter Modelica.SIunits.CoefficientOfHeatTransfer gamma_G = 200 "Constant heat transfer coefficient in the gas side";
  parameter Modelica.SIunits.CoefficientOfHeatTransfer gamma_F = 200 "Constant heat transfer coefficient in the fluid side";
  parameter Modelica.SIunits.Temperature Tstart_M_In = Tstart_F_In "Inlet metal wall temperature start value";
  parameter Modelica.SIunits.Temperature Tstart_M_Out = Tstart_F_Out "Outlet metal wall temperature start value";
  parameter Modelica.SIunits.Temperature Tstart_F_In = 633 "Inlet fluid temperature start value";
  parameter Modelica.SIunits.Temperature Tstart_F_Out = 843 "Outlet fluid temperature start value";

  // use FlowJoin to mix flow
  Gas.FlowJoin mixer(redeclare package Medium = medium_cold);  

  /*
  ThermoPower.Gas.SensT T_waterOut(
    redeclare package Medium = medium_cold) 
  annotation(
    Placement(transformation(origin = {4, -50}, extent = {{-10, -10}, {10, 10}}, rotation = 270)));

  ThermoPower.Gas.SensT T_gasOut(
    redeclare package Medium = medium_hot) 
  annotation(
    Placement(transformation(extent = {{30, -6}, {50, 14}}, rotation = 0)));
  */
/*
  // ThermoPower.PowerPlants.HRSG.Components.HEG2G HTR(
  // ThermoPower.PowerPlants.HRSG.Components.HE HTR(
  // Steps.TPComponents.HE HTR(
  Steps.TPComponents.HEG2G HTR(
    redeclare package FluidMedium = medium_cold, 
    redeclare package FlueGasMedium = medium_hot, 
    redeclare replaceable model HeatTransfer_F = ThermoPower.Thermal.HeatTransferFV.ConstantHeatTransferCoefficient(gamma = gamma_F), 
    redeclare replaceable model HeatTransfer_G = ThermoPower.Thermal.HeatTransferFV.ConstantHeatTransferCoefficientTwoGrids(gamma = gamma_G),
    redeclare model HeatExchangerTopology = ThermoPower.Thermal.HeatExchangerTopologies.CounterCurrentFlow, 
    N_F = 7, 
    N_G = 7, 
    Nw_G = 6, 
    SSInit = true, 
    Tstartbar_G = 823.15, 
    Tstartbar_M = 773.15, 
    exchSurface_F = 225.073, 
    exchSurface_G = 1708.2, 
    extSurfaceTub = 252.286, 
    fluidNomFlowRate = fluidNomFlowRate, 
    fluidNomPressure = fluidNomPressure, 
    fluidVol = 2.234, 
    gasNomFlowRate = gasNomFlowRate, 
    gasNomPressure = gasNomPressure, 
    // fluidQuasiStatic = true,
    gasVol = 10, 
    lambda = 20, 
    metalVol = 0.573, 
    pstart_F = fluidNomPressure, 
    rhomcm = 7900 * 578.05) 
    annotation(
      Placement(transformation(extent = {{-20, -20}, {20, 20}}, rotation = 0)));
 
  // ThermoPower.PowerPlants.HRSG.Components.HEG2G LTR(
  // ThermoPower.PowerPlants.HRSG.Components.HE LTR(
  // Steps.TPComponents.HE LTR(
  Steps.TPComponents.HEG2G LTR(
    redeclare package FluidMedium = medium_cold, 
    redeclare package FlueGasMedium = medium_hot, 
    redeclare replaceable model HeatTransfer_F = ThermoPower.Thermal.HeatTransferFV.ConstantHeatTransferCoefficient(gamma = gamma_F), 
    redeclare replaceable model HeatTransfer_G = ThermoPower.Thermal.HeatTransferFV.ConstantHeatTransferCoefficientTwoGrids(gamma = gamma_G),
    redeclare model HeatExchangerTopology = ThermoPower.Thermal.HeatExchangerTopologies.CounterCurrentFlow, 
    N_F = 7, 
    N_G = 7, 
    Nw_G = 6, 
    SSInit = true, 
    Tstartbar_G = 823.15, 
    Tstartbar_M = 773.15, 
    exchSurface_F = 225.073, 
    exchSurface_G = 1708.2, 
    extSurfaceTub = 252.286, 
    fluidNomFlowRate = fluidNomFlowRate, 
    fluidNomPressure = fluidNomPressure, 
    // fluidQuasiStatic = true,
    fluidVol = 2.234, 
    gasNomFlowRate = gasNomFlowRate, 
    gasNomPressure = gasNomPressure, 
    gasVol = 10, 
    lambda = 20, 
    metalVol = 0.573, 
    pstart_F = fluidNomPressure, 
    rhomcm = 7900 * 578.05) 
    annotation(
      Placement(transformation(extent = {{-20, -20}, {20, 20}}, rotation = 0)));
*/

  Steps.TPComponents.PCHE LTR(
    redeclare package FluidMedium = medium_cold, 
    redeclare package FlueGasMedium = medium_hot, 
     
    // use Marchionni PCHE HeatTransfer
    // slow but can have a result - set a_phi = 0 to use Gnielinski's correlation 
    // redeclare replaceable model HeatTransfer_F = Steps.TPComponents.MarchionniPCHEHeatTransferFV(),
    // redeclare replaceable model HeatTransfer_G = Steps.TPComponents.MarchionniPCHEHeatTransferFV(),
    // gasFlow(heatTransfer(pitch = cfg.pitch, phi = cfg.phi, Cf_C1 = Cf_C1, Cf_C2 = Cf_C2, Cf_C3 = Cf_C3)),
    // fluidFlow(heatTransfer(pitch = cfg.pitch, phi = cfg.phi, Cf_C1 = Cf_C1, Cf_C2 = Cf_C2, Cf_C3 = Cf_C3)),    
    
    // fast and works fine for now. Error occurs when mass flow rate is zero, i.e. one flow is shut down. 
    redeclare replaceable model HeatTransfer_F = ThermoPower.Thermal.HeatTransferFV.IdealHeatTransfer,      
    redeclare replaceable model HeatTransfer_G = ThermoPower.Thermal.HeatTransferFV.IdealHeatTransfer, 
    
    bc = bc_LTR, 
    geo_hot = cfg.cfg_LTR_hot.geo,
    geo_cold = cfg.cfg_LTR_cold.geo,
    geo_tube = cfg.cfg_LTR_tube.geo,  
    thermo_hot = cfg.cfg_LTR_hot.thermo,
    thermo_cold = cfg.cfg_LTR_cold.thermo,
    thermo_tube = cfg.cfg_LTR_tube.thermo,   
    L = L_fp,
    SSInit = true,
    gasQuasiStatic = false,
    fluidQuasiStatic = false,
    metalWall(L = L_wall, w_ch = W_ch, h_ch = H_ch, dx = T_wall),
    table_k_metalwall =   table_k_metalwall
    // metalQuasiStatic = true
    // override the values of Am and L of metaltubeFV
    // to make them agree with semi-circular tube of PCHE
    // ('final' modifier of Am in metalTubeFv was removed as well)
    //metalTube(WallRes=false, L = 1, rhomcm=200, Am = HE.metalVol / 1) 
  )
  annotation(
    Placement(transformation(extent = {{-20, -20}, {20, 20}}, rotation = 0)));

   
  Steps.TPComponents.PCHE HTR(
    redeclare package FluidMedium = medium_cold, 
    redeclare package FlueGasMedium = medium_hot, 
     
    // use Marchionni PCHE HeatTransfer
    // slow but can have a result - set a_phi = 0 to use Gnielinski's correlation 
    // redeclare replaceable model HeatTransfer_F = Steps.TPComponents.MarchionniPCHEHeatTransferFV(),
    // redeclare replaceable model HeatTransfer_G = Steps.TPComponents.MarchionniPCHEHeatTransferFV(),
    // gasFlow(heatTransfer(pitch = cfg.pitch, phi = cfg.phi, Cf_C1 = Cf_C1, Cf_C2 = Cf_C2, Cf_C3 = Cf_C3)),
    // fluidFlow(heatTransfer(pitch = cfg.pitch, phi = cfg.phi, Cf_C1 = Cf_C1, Cf_C2 = Cf_C2, Cf_C3 = Cf_C3)),    
    
    // fast and works fine for now. Error occurs when mass flow rate is zero, i.e. one flow is shut down. 
    redeclare replaceable model HeatTransfer_F = ThermoPower.Thermal.HeatTransferFV.IdealHeatTransfer,      
    redeclare replaceable model HeatTransfer_G = ThermoPower.Thermal.HeatTransferFV.IdealHeatTransfer,  
    bc = bc_HTR, 
    geo_hot = cfg.cfg_HTR_hot.geo,
    geo_cold = cfg.cfg_HTR_cold.geo,
    geo_tube = cfg.cfg_HTR_tube.geo,  
    thermo_hot = cfg.cfg_HTR_hot.thermo,
    thermo_cold = cfg.cfg_HTR_cold.thermo,
    thermo_tube = cfg.cfg_HTR_tube.thermo,   
    L = L_fp,
    SSInit = true,
    gasQuasiStatic = false,
    fluidQuasiStatic = false,
    metalWall(L = L_wall, w_ch = W_ch, h_ch = H_ch, dx = T_wall),
    table_k_metalwall = table_k_metalwall
    // metalQuasiStatic = true
    // override the values of Am and L of metaltubeFV
    // to make them agree with semi-circular tube of PCHE
    // ('final' modifier of Am in metalTubeFv was removed as well)
    //metalTube(WallRes=false, L = 1, rhomcm=200, Am = HE.metalVol / 1) 
  )
  annotation(
    Placement(transformation(extent = {{-20, -20}, {20, 20}}, rotation = 0)));

/* 
  ThermoPower.Gas.Turbine Turbine1(
  redeclare package Medium = medium_hot, 
  fileName = Modelica.Utilities.Files.loadResource("modelica://Steps/Resources/Data/turbine_map.txt"),   
  tablePhic = fill(0.0, 14, 12), //tablePhic, 
  tableEta = fill(0.0, 14, 12), //tableEta, 
  pstart_in = bc_heater.st_cold_out.p, 
  pstart_out = bc_HTR.st_hot_in.p, 
  Tstart_in = bc_heater.st_cold_out.T, 
  Tstart_out = bc_HTR.st_hot_in.T, 
  Ndesign = 60000.0, 
  Tdes_in = bc_heater.st_cold_out.T,  
  Table = ThermoPower.Choices.TurboMachinery.TableTypes.file,
  explicitIsentropicEnthalpy = false,
  gas_in(
    p(nominal = Turbine1.pstart_in), 
    T(nominal = Turbine1.Tstart_in)),
  gas_iso(
    p(nominal = Turbine1.pstart_out), 
    T(nominal = Turbine1.Tstart_out))) 
  annotation(
    Placement(transformation(extent = {{-40, -20}, {0, 20}}, rotation = 0)));

  Modelica.Mechanics.Rotational.Sources.Speed speed1 annotation(
    Placement(visible = true, transformation(origin = {84, 0}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.Constant const1(k = 60000) annotation(
    Placement(visible = true, transformation(origin = {130, 0}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));  
 
  ThermoPower.Gas.SensT sens_turbine(redeclare package Medium = medium_hot) annotation(
    Placement(visible = true, transformation(origin = {20, 44}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
*/
/*
  // Input signals for transient simulation
  // hot/gas side
  Modelica.Blocks.Sources.IntegerConstant const_T_step_h(k = 20) annotation(
    Placement(visible = true, transformation(origin = {-230, 10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Interaction.Show.IntegerValue disp_T_h annotation(
    Placement(visible = true, transformation(origin = {-108, 44}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.IntegerConstant const_T_offset_h(k = 630) annotation(
    Placement(visible = true, transformation(origin = {-180, -54}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.MathInteger.TriggeredAdd triadd_T_h annotation(
    Placement(visible = true, transformation(origin = {-182, 8}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
  Modelica.Blocks.MathInteger.Sum sum_T_h(nu = 2) annotation(
    Placement(visible = true, transformation(origin = {-144, 8}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.BooleanPulse en_triadd_T_h(period = 10, width = 10) annotation(
    Placement(visible = true, transformation(origin = {-224, -24}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Interaction.Show.IntegerValue disp_T_step_h annotation(
    Placement(visible = true, transformation(origin = {-140, 42}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Math.IntegerToReal I2R_T_h annotation(
    Placement(visible = true, transformation(origin = {-102, 8}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));

  // cold/fluid side
  Modelica.Blocks.Sources.IntegerConstant const_T_offset_c(k = 400) annotation(
    Placement(visible = true, transformation(origin = {-128, 72}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Interaction.Show.IntegerValue disp_T_c annotation(
    Placement(visible = true, transformation(origin = {-144, 158}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.MathInteger.TriggeredAdd triadd_T_c annotation(
    Placement(visible = true, transformation(origin = {-128, 126}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
  Modelica.Blocks.Sources.BooleanPulse en_triadd_T_c(period = 10, startTime = 5, width = 10) annotation(
    Placement(visible = true, transformation(origin = {-156, 102}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.MathInteger.Sum sum_T_c(nu = 2) annotation(
    Placement(visible = true, transformation(origin = {-94, 126}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Interaction.Show.IntegerValue disp_T_step_c annotation(
    Placement(visible = true, transformation(origin = {-56, 158}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Math.IntegerToReal I2R_T_c annotation(
    Placement(visible = true, transformation(origin = {-48, 124}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.IntegerConstant const_T_step_c(k = 20) annotation(
    Placement(visible = true, transformation(origin = {-194, 126}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
*/
protected
  parameter Real tablePhic[5, 4] = [1, 37, 80, 100; 1.5, 7.10E-05, 7.10E-05, 7.10E-05; 2, 8.40E-05, 8.40E-05, 8.40E-05; 2.5, 8.70E-05, 8.70E-05, 8.70E-05; 3, 1.04E-04, 1.04E-04, 1.04E-04];
  parameter Real tableEta[5, 4] = [1, 37, 80, 100; 1.5, 0.57, 0.89, 0.81; 2, 0.46, 0.82, 0.88; 2.5, 0.41, 0.76, 0.85; 3, 0.38, 0.72, 0.82];
  
equation
 
/*
  //HTR + mixer

  connect(source_mixer_in.flange, mixer.inlet1);  
  connect(source_cold.flange, mixer.inlet2);  
  connect(mixer.outlet, HTR.waterIn);

  connect(HTR.waterOut, T_waterOut.inlet) annotation(
    Line(points = {{8.88178e-016, -44}, {8.88178e-016, -20}, {0, -20}}, thickness = 0.5, color = {0, 0, 255}));
  connect(sink_cold.flange, T_waterOut.outlet) annotation(
    Line(points = {{1.83697e-015, -70}, {1.83697e-015, -56}, {-8.88178e-016, -56}}, thickness = 0.5, color = {0, 0, 255}));

  // connect(HTR.waterOut, sink_cold.flange) 
   //annotation(Line(points = {{1.83697e-015, -70}, {1.83697e-015, -56}, {-8.88178e-016, -56}}, thickness = 0.5, color = {0, 0, 255}));   
  
  connect(source_hot.flange, HTR.gasIn) annotation(
    Line(points = {{-50, 0}, {-20, 0}}, color = {159, 159, 223}, thickness = 0.5, smooth = Smooth.None));   
    
  connect(HTR.gasOut, T_gasOut.inlet ) annotation(
    Line(points = {{34, 0}, {34, 0}, {20, 0}}, color = {159, 159, 223}, thickness = 0.5));
  connect(T_gasOut.outlet, sink_hot.flange) annotation(
    Line(points = {{46, 0}, {46, 0}, {60, 0}}, color = {159, 159, 223}, thickness = 0.5));
  
  // connect(HTR.gasOut, sink_hot.flange) annotation(
  //   Line(points = {{46, 0}, {46, 0}, {60, 0}}, color = {159, 159, 223}, thickness = 0.5));
*/  
/*
  // mixer + LTR
  // water/cold side  
  connect(source_mixer_in.flange, mixer.inlet1);
  
  connect(source_cold.flange, LTR.waterIn);
  
  connect(LTR.waterOut, mixer.inlet2);
  
  connect(mixer.outlet, T_waterOut.inlet);
  
  connect(T_waterOut.outlet, sink_cold.flange) annotation(
    Line(points = {{1.83697e-015, -70}, {1.83697e-015, -56}, {-8.88178e-016, -56}}, thickness = 0.5, color = {0, 0, 255}));
  
  // gas/hot side
  connect(source_hot.flange, LTR.gasIn);
  
  connect(LTR.gasOut, T_gasOut.inlet);
  connect(T_gasOut.outlet, sink_hot.flange) annotation(
    Line(points = {{46, 0}, {46, 0}, {60, 0}}, color = {159, 159, 223}, thickness = 0.5));
*/


  //HTR + mixer + LTR 
  // water/cold side    
  connect(source_cold.flange, LTR.waterIn);
  
  connect(source_mixer_in.flange, mixer.inlet1);  
  connect(LTR.waterOut, mixer.inlet2);  
  connect(mixer.outlet, HTR.waterIn);
  
  // connect(LTR.waterOut, HTR.waterIn);

  connect(HTR.waterOut, sink_cold.flange) annotation(
    Line(points = {{1.83697e-015, -70}, {1.83697e-015, -56}, {-8.88178e-016, -56}}, thickness = 0.5, color = {0, 0, 255}));
  
  // gas/hot side
  connect(source_hot.flange, HTR.gasIn) annotation(
   Line(points = {{-50, 0}, {-20, 0}}, color = {159, 159, 223}, thickness = 0.5, smooth = Smooth.None));

  connect(HTR.gasOut, LTR.gasIn);
  
  connect(LTR.gasOut, sink_hot.flange) annotation(
    Line(points = {{46, 0}, {46, 0}, {60, 0}}, color = {159, 159, 223}, thickness = 0.5));

/*
  //LTR + HTR + mixer
  // water/cold side    
  connect(source_mixer_in.flange, mixer.inlet1);  
  connect(source_cold.flange, mixer.inlet2);  
  connect(mixer.outlet, HTR.waterIn);
  
  // connect(LTR.waterOut, HTR.waterIn);

  connect(HTR.waterOut, LTR.waterIn);
  connect(LTR.waterOut, sink_cold.flange) annotation(
    Line(points = {{1.83697e-015, -70}, {1.83697e-015, -56}, {-8.88178e-016, -56}}, thickness = 0.5, color = {0, 0, 255}));
  
  // gas/hot side
  connect(source_hot.flange, LTR.gasIn) annotation(
   Line(points = {{-50, 0}, {-20, 0}}, color = {159, 159, 223}, thickness = 0.5, smooth = Smooth.None));

  connect(LTR.gasOut, HTR.gasIn);
  
  connect(HTR.gasOut, sink_hot.flange) annotation(
    Line(points = {{46, 0}, {46, 0}, {60, 0}}, color = {159, 159, 223}, thickness = 0.5));
*/

/*
  // temperature input
  // hot / gas side
  connect(const_T_step_h.y, triadd_T_h.u) annotation(
    Line(points = {{-219, 10}, {-202.5, 10}, {-202.5, 8}, {-190, 8}}, color = {255, 127, 0}));
  connect(const_T_offset_h.y, sum_T_h.u[2]) annotation(
    Line(points = {{-169, -54}, {-154, -54}, {-154, 8}}, color = {255, 127, 0}));
  connect(sum_T_h.y, disp_T_h.numberPort) annotation(
    Line(points = {{-132.5, 8}, {-121, 8}, {-121, 44}, {-119.5, 44}}, color = {255, 127, 0}));
  connect(triadd_T_h.y, disp_T_step_h.numberPort) annotation(
    Line(points = {{-175, 8}, {-160.5, 8}, {-160.5, 42}, {-151.5, 42}}, color = {255, 127, 0}));
  connect(triadd_T_h.y, sum_T_h.u[1]) annotation(
    Line(points = {{-175, 8}, {-154, 8}}, color = {255, 127, 0}));
  connect(en_triadd_T_h.y, triadd_T_h.trigger) annotation(
    Line(points = {{-213, -24}, {-186, -24}, {-186, 1}}, color = {255, 0, 255}));
  connect(sum_T_h.y, I2R_T_h.u) annotation(
    Line(points = {{-132.5, 8}, {-114, 8}}, color = {255, 127, 0}));
  connect(I2R_T_h.y, source_hot.in_T) annotation(
    Line(points = {{-91, 8}, {-90, 8}, {-90, 9}, {-80, 9}, {-80, 7.75}, {-78, 7.75}, {-78, 6}}, color = {0, 0, 127}));

  // cold / fluid side 
  connect(en_triadd_T_c.y, triadd_T_c.trigger) annotation(
    Line(points = {{-145, 102}, {-132, 102}, {-132, 119}}, color = {255, 0, 255}));
  connect(triadd_T_c.y, sum_T_c.u[1]) annotation(
    Line(points = {{-121, 126}, {-104, 126}}, color = {255, 127, 0}));
  connect(triadd_T_c.y, disp_T_c.numberPort) annotation(
    Line(points = {{-121, 126}, {-160.5, 126}, {-160.5, 158}, {-155.5, 158}}, color = {255, 127, 0}));
  connect(sum_T_c.y, I2R_T_c.u) annotation(
    Line(points = {{-82.5, 126}, {-71.25, 126}, {-71.25, 124}, {-60, 124}}, color = {255, 127, 0}));
  connect(sum_T_c.y, disp_T_step_c.numberPort) annotation(
    Line(points = {{-82.5, 126}, {-71, 126}, {-71, 158}, {-67.5, 158}}, color = {255, 127, 0}));
  connect(const_T_step_c.y, triadd_T_c.u) annotation(
    Line(points = {{-183, 126}, {-136, 126}}, color = {255, 127, 0}));
  connect(const_T_offset_c.y, sum_T_c.u[2]) annotation(
    Line(points = {{-117, 72}, {-104, 72}, {-104, 126}}, color = {255, 127, 0}));
  connect(I2R_T_c.y, source_cold.in_T) annotation(
    Line(points = {{-37, 124}, {6, 124}, {6, 64}}, color = {0, 0, 127}));
*/

annotation(
    Diagram(graphics),
    experiment(StartTime = 0, StopTime = 20, Tolerance = 1e-3, Interval = 1),    
    // options = "-showErrorMessages -demoMode",
    __OpenModelica_commandLineOptions = "--matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection -d=initialization,NLSanalyticJacobian,aliasConflicts,dumpCSE",    
    __OpenModelica_simulationFlags(lv = "LOG_DEBUG,LOG_NLS,LOG_NLS_V,LOG_STATS,LOG_INIT,LOG_STDOUT, -w", outputFormat = "mat", s = "dassl", nls = "homotopy"));
end TestDyn_Comps_PCHE_ramp;
