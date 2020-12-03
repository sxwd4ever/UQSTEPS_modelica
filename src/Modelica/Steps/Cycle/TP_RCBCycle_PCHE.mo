within Steps.Cycle;

model TP_RCBCycle_PCHE
  "ThermoPower based RCBCycle with PCHE as recuperator"
import Modelica.SIunits.Conversions.{from_degC,from_deg};
  import Modelica.SIunits.{Temperature,Pressure,SpecificEnthalpy};
  import Util = Utilities.Util;
  import Steps.Utilities.CoolProp.PropsSI;
  import Steps.Components.{PCHEGeoParam};
  import Steps.Model.{PBConfiguration,SimParam,EntityConfig,EntityGeoParam,EntityThermoParam,ThermoState,HEBoundaryCondition};
  import Model.PBConfiguration;
  import ThermoPower.Choices.Init.Options;
  import ThermoPower.System;
  import ThermoPower.Gas;
  //package medium_main = Modelica.Media.IdealGases.SingleGases.CO2; //Steps.Media.CO2; //
  //package medium_main = Modelica.Media.IdealGases.SingleGases.CO2; //Steps.Media.CO2; Steps.Media.CO2; //
  package medium_main = Steps.Media.SCO2;
  package medium_heater = SolarTherm.Media.Sodium.Sodium_pT;
  // Modelica.Media.IdealGases.SingleGases.CO2;
  package medium_cooler = ThermoPower.Water.StandardWater;
  
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
    mdot_heater = mdot_heater_hot,
    T_heater_hot_in = T_heater_hot,
    T_heater_hot_out = T_heater_hot - 200,
    T_heater_cold_out = T_heater_hot - 90,
    T_cooler_cold_in = T_cooler_cold
  );       
  
  // set boundary conditions
  parameter HEBoundaryCondition bc_HTR = cfg.bc_HTR;  
  parameter HEBoundaryCondition bc_LTR = cfg.bc_LTR;
  parameter HEBoundaryCondition bc_heater = cfg.bc_heater;
  parameter HEBoundaryCondition bc_cooler = cfg.bc_cooler;  
  parameter ThermoState st_bypass = cfg.st_bypass;
  
  // load thermodynamic parameters.     
  parameter EntityThermoParam thermo_heater_tube = cfg.cfg_heater_tube.thermo;
  parameter EntityThermoParam thermo_cooler_tube = cfg.cfg_cooler_tube.thermo;  
  parameter EntityThermoParam thermo_mixer = cfg.cfg_mixer.thermo;
    
  //Components
  inner ThermoPower.System system(allowFlowReversal = false, initOpt=ThermoPower.Choices.Init.Options.noInit) annotation(
    Placement(transformation(extent = {{100, 80}, {120, 100}})));

  // global init opition (system.initOpt) leads to order reduction error
  // use this flag to control the initialization of all components instead. 
  parameter Boolean SSInit = false "Steady-state initialization";
 
  ThermoPower.Water.SourceMassFlow source_heater_hot(
    redeclare package Medium = medium_heater, 
    w0 = bc_heater.st_hot_in.mdot,
    p0 = bc_heater.st_hot_in.p,
    h = bc_heater.st_hot_in.h,
    T = bc_heater.st_hot_in.T)  annotation(
    Placement(visible = true, transformation(origin = {-80, 22}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
    
  ThermoPower.Water.SinkPressure sink_heater_hot(
    redeclare package Medium = medium_heater, 
    p0 = bc_heater.st_hot_out.p, 
    T = bc_heater.st_hot_out.T,
    use_T = true) annotation(
    Placement(visible = true, transformation(origin = {-81, 63}, extent = {{-5, -5}, {5, 5}}, rotation = 180)));
    
  TPComponents.HE heater(
    redeclare package FluidMedium = medium_heater, 
    redeclare package FlueGasMedium = medium_main, 
    fluidFlow(fixedMassFlowSimplified = true, hstartin = bc_heater.st_hot_in.h, hstartout=bc_heater.st_hot_out.h), // set the fluid flow as fixed mdot for simplarity
    gasFlow(QuasiStatic = true, Tstartin = bc_heater.st_cold_in.T, Tstartout = bc_heater.st_cold_out.T),
    bc = bc_heater, 
    geo_hot = cfg.cfg_heater_hot.geo,
    geo_cold = cfg.cfg_heater_cold.geo,
    geo_tube = cfg.cfg_heater_tube.geo,  
    thermo_hot = cfg.cfg_heater_hot.thermo,
    thermo_cold = cfg.cfg_heater_cold.thermo,
    thermo_tube = cfg.cfg_heater_tube.thermo,      
    SSInit=SSInit,
    FluidPhaseStart=ThermoPower.Choices.FluidPhase.FluidPhases.Liquid,    
    redeclare replaceable model HeatTransfer_F =  ThermoPower.Thermal.HeatTransferFV.IdealHeatTransfer, // ConstantHeatTransferCoefficientTwoGrids(gamma = thermo_HTR_hot.gamma_he),    
    redeclare replaceable model HeatTransfer_G = ThermoPower.Thermal.HeatTransferFV.IdealHeatTransfer, // ConstantHeatTransferCoefficient(gamma = thermo_HTR_cold.gamma_he),     
    redeclare model HeatExchangerTopology = ThermoPower.Thermal.HeatExchangerTopologies.CounterCurrentFlow,
    metalTube(WallRes=false, Tstartbar=heater.Tstartbar_M)) annotation(
    Placement(visible = true, transformation(origin = {-61, 41}, extent = {{-7, -7}, {7, 7}}, rotation = 180)));

  // use FlowJoin to mix flow
  Gas.FlowJoin mixer(redeclare package Medium = medium_main) annotation(
    Placement(visible = true, transformation(origin = {-12, 60}, extent = {{-6, -6}, {6, 6}}, rotation = 180)));

  //Real kim_cor_coe[4] = {kim.a, kim.b, kim.c, kim.d};
  parameter SI.Length pitch = 12.3e-3 "pitch length";
  parameter Real phi = 35 "pitch angle °";

  TPComponents.PCHE HTR(
    redeclare package FluidMedium = medium_main, 
    redeclare package FlueGasMedium = medium_main, 
    gasFlow(QuasiStatic = true, Tstartin = bc_HTR.st_hot_in.T, Tstartout = bc_HTR.st_hot_out.T),
    fluidFlow(Tstartin = bc_HTR.st_cold_in.T, Tstartout = bc_HTR.st_cold_out.T),  
    redeclare replaceable model HeatTransfer_F = Steps.TPComponents.KimPCHEHeatTransferFV(
    pitch = pitch,
    phi = phi), 
    redeclare replaceable model HeatTransfer_G = Steps.TPComponents.KimPCHEHeatTransferFV(
    pitch = pitch,
    phi = phi), 
     //redeclare replaceable model HeatTransfer_G = ThermoPower.Thermal.HeatTransferFV.IdealHeatTransfer,
    redeclare model HeatExchangerTopology = ThermoPower.Thermal.HeatExchangerTopologies.CounterCurrentFlow,  
    bc = bc_HTR, 
    geo_hot = cfg.cfg_HTR_hot.geo,
    geo_cold = cfg.cfg_HTR_cold.geo,
    geo_tube = cfg.cfg_HTR_tube.geo,  
    thermo_hot = cfg.cfg_HTR_hot.thermo,
    thermo_cold = cfg.cfg_HTR_cold.thermo,
    thermo_tube = cfg.cfg_HTR_tube.thermo, 
    gasQuasiStatic = true,
    fluidQuasiStatic = true,
    SSInit=SSInit,
    metalTube(WallRes=false, Tstartbar=HTR.Tstartbar_M)) annotation(
    Placement(visible = true, transformation(origin = {-21, 33}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
    
  TPComponents.PCHE LTR(
    redeclare package FluidMedium = medium_main, 
    redeclare package FlueGasMedium = medium_main, 
    gasFlow(QuasiStatic = true, Tstartin = bc_LTR.st_hot_in.T, Tstartout = bc_LTR.st_hot_out.T),
    fluidFlow(Tstartin = bc_LTR.st_cold_in.T, Tstartout = bc_LTR.st_cold_out.T),   
    // redeclare replaceable model HeatTransfer_F = ThermoPower.Thermal.HeatTransferFV.IdealHeatTransfer, 
    redeclare replaceable model HeatTransfer_F = Steps.TPComponents.KimPCHEHeatTransferFV(
    pitch = pitch,
    phi = phi),
    // ThermoPower.Thermal.HeatTransferFV.ConstantHeatTransferCoefficient(gamma = thermo_LTR_cold.gamma_he),
    // redeclare replaceable model HeatTransfer_G = ThermoPower.Thermal.HeatTransferFV.IdealHeatTransfer, 
    redeclare replaceable model HeatTransfer_G = Steps.TPComponents.KimPCHEHeatTransferFV(
    pitch = pitch,
    phi = phi),  
    // ThermoPower.Thermal.HeatTransferFV.ConstantHeatTransferCoefficientTwoGrids(gamma = thermo_LTR_hot.gamma_he),
    redeclare model HeatExchangerTopology = ThermoPower.Thermal.HeatExchangerTopologies.CounterCurrentFlow,  
    bc = bc_LTR, 
    geo_hot = cfg.cfg_LTR_hot.geo,
    geo_cold = cfg.cfg_LTR_cold.geo,
    geo_tube = cfg.cfg_LTR_tube.geo,  
    thermo_hot = cfg.cfg_LTR_hot.thermo,
    thermo_cold = cfg.cfg_LTR_cold.thermo,
    thermo_tube = cfg.cfg_LTR_tube.thermo,  
    gasQuasiStatic = true,
    fluidQuasiStatic = true,
    SSInit=SSInit,
    metalTube(WallRes=false, Tstartbar=LTR.Tstartbar_M)) annotation(
    Placement(visible = true, transformation(origin = {21, 33}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
  
  ThermoPower.Gas.Turbine turbine(
  redeclare package Medium = medium_main, 
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
    p(nominal = turbine.pstart_in), 
    T(nominal = turbine.Tstart_in)),
  gas_iso(
    p(nominal = turbine.pstart_out), 
    T(nominal = turbine.Tstart_out))) annotation(
    Placement(visible = true, transformation(origin = {-69, -51}, extent = {{-11, -11}, {11, 11}}, rotation = 0)));
  
  Modelica.Mechanics.Rotational.Sources.ConstantSpeed const_speed_turbine(
      w_fixed=60000.0, useSupport=false) annotation(
    Placement(visible = true, transformation(origin = {-46, -50}, extent = {{6, -6}, {-6, 6}}, rotation = 0)));
  
  ThermoPower.Water.SourceMassFlow source_cooler_cold(
  redeclare package Medium = medium_cooler,
    T = bc_cooler.st_cold_in.T, 
    p0 = bc_cooler.st_cold_in.p, 
    use_T = true,
    use_in_T = false, 
    w0 = bc_cooler.st_cold_in.mdot) annotation(
    Placement(visible = true, transformation(origin = {29, -37}, extent = {{-5, -5}, {5, 5}}, rotation = 0)));
    
  ThermoPower.Water.SinkPressure sink_cooler_cold(
  redeclare package Medium = medium_cooler,
    p0 = bc_cooler.st_cold_out.p, 
    T = bc_cooler.st_cold_out.T,
    use_T = true) annotation(
    Placement(visible = true, transformation(origin = {31, -75}, extent = {{-5, -5}, {5, 5}}, rotation = 180)));
      
  TPComponents.HE cooler(
    //Components.HEG2G cooler(
      redeclare package FluidMedium = medium_cooler, 
      redeclare package FlueGasMedium = medium_main, 
      fluidFlow(fixedMassFlowSimplified = true, hstartin = bc_cooler.st_cold_in.h, hstartout=bc_cooler.st_cold_out.h), // set the fluid flow as fixed mdot for simplarity
      gasFlow(QuasiStatic = true, Tstartin = bc_cooler.st_hot_in.T, Tstartout = bc_cooler.st_hot_out.T),
      bc = bc_cooler, 
      geo_hot = cfg.cfg_cooler_hot.geo,
      geo_cold = cfg.cfg_cooler_cold.geo,
      geo_tube = cfg.cfg_cooler_tube.geo,  
      thermo_hot = cfg.cfg_cooler_hot.thermo,
      thermo_cold = cfg.cfg_cooler_cold.thermo,
      thermo_tube = cfg.cfg_cooler_tube.thermo, 
      SSInit=SSInit,
      FluidPhaseStart=ThermoPower.Choices.FluidPhase.FluidPhases.Liquid,    
      redeclare replaceable model HeatTransfer_F = ThermoPower.Thermal.HeatTransferFV.IdealHeatTransfer, //ConstantHeatTransferCoefficientTwoGrids(gamma = thermo_hot.gamma_he),     
      redeclare replaceable model HeatTransfer_G = ThermoPower.Thermal.HeatTransferFV.IdealHeatTransfer, //ConstantHeatTransferCoefficient(gamma =  thermo_cold.gamma_he),     
      redeclare model HeatExchangerTopology = ThermoPower.Thermal.HeatExchangerTopologies.CounterCurrentFlow,
      metalTube(WallRes=false,Tstartbar=bc_cooler.st_hot_in.T - 50)) annotation(
    Placement(visible = true, transformation(origin = {47, -57}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
    
  ThermoPower.Gas.FlowSplit splitter(redeclare package Medium = medium_main) annotation(
    Placement(visible = true, transformation(origin = {17, 1}, extent = {{-5, -5}, {5, 5}}, rotation = 0)));
    
  ThermoPower.Gas.Compressor compressor(
    redeclare package Medium = medium_main,
    pstart_in=bc_cooler.st_hot_out.p,
    pstart_out=bc_LTR.st_cold_in.p,
    Tstart_in=bc_cooler.st_hot_out.T,
    Tstart_out=bc_LTR.st_cold_in.T,
    tablePhic=tablePhic_comp_mc,
    tableEta=tableEta_comp_mc,
    tablePR=tablePR_comp_mc,
    Table=ThermoPower.Choices.TurboMachinery.TableTypes.matrix,
    Ndesign=523.3,
    Tdes_in=compressor.Tstart_in) annotation(
    Placement(visible = true, transformation(origin = {103, -11}, extent = {{-11, -11}, {11, 11}}, rotation = 0)));

  Modelica.Mechanics.Rotational.Sources.ConstantSpeed const_speed_comp(
      w_fixed=523.3, useSupport=false) annotation(
    Placement(visible = true, transformation(origin = {81, -7}, extent = {{-5, -5}, {5, 5}}, rotation = 0)));

  ThermoPower.Gas.Compressor recompressor(
    redeclare package Medium = medium_main,
    pstart_in=bc_HTR.st_hot_out.p,
    pstart_out=st_bypass.p,
    Tstart_in=bc_HTR.st_hot_out.T,
    Tstart_out=st_bypass.T,
    tablePhic=tablePhic_comp_rc,
    tableEta=tableEta_comp_rc,
    tablePR=tablePR_comp_rc,
    Table=ThermoPower.Choices.TurboMachinery.TableTypes.matrix,
    Ndesign=523.3,
    Tdes_in=recompressor.Tstart_in) annotation(
    Placement(visible = true, transformation(origin = {57, -5}, extent = {{-11, -11}, {11, 11}}, rotation = 0)));
 
  Modelica.Mechanics.Rotational.Sources.ConstantSpeed const_speed_recomp(
      w_fixed=523.3, useSupport=false) annotation(
    Placement(visible = true, transformation(origin = {33, -7}, extent = {{-5, -5}, {5, 5}}, rotation = 0)));

  ThermoPower.Gas.SourceMassFlow source(
  redeclare package Medium = medium_main,
    w0 = bc_LTR.st_hot_out.mdot,
    p0 = bc_LTR.st_hot_out.p,
    T = bc_LTR.st_hot_out.T,
    use_in_T = false,
    use_in_w0 = false,
    //use_in_p0 = true,
    gas(T(nominal = source.T), p(nominal=source.p0)),
    flange(p(nominal=source.p0), m_flow(nominal=source.w0))) annotation(
    Placement(visible = true, transformation(origin = {-91, -17}, extent = {{-5, -5}, {5, 5}}, rotation = -90)));
    
  ThermoPower.Gas.SinkPressure sink(
  redeclare package Medium = medium_main,
    use_in_T = true,
    use_in_p0 = true,  
    p0 = bc_LTR.st_hot_out.p,
    T = bc_LTR.st_hot_out.T) annotation(
    Placement(visible = true, transformation(origin = {-92, -4}, extent = {{-4, -4}, {4, 4}}, rotation = 270)));
  /*
  ThermoPower.Gas.Utility.ClosedSystemInit sys_init(
    redeclare package Medium = medium_main,
    //w_b(nominal = bc_heater.st_cold_out.mdot),
    pstart = bc_heater.st_cold_out.p);
  */
  /*
  ThermoPower.Gas.ThroughMassFlow source_through_cooler(
    redeclare package Medium = medium_main,
    w0 = bc_cooler.st_hot_out.mdot) annotation(
    Placement(visible = true, transformation(origin = {88, 12}, extent = {{-6, -6}, {6, 6}}, rotation = 90)));

  
  ThermoPower.Gas.ThroughMassFlow source_through_heater(
    redeclare package Medium = medium_main,
    w0 = bc_heater.st_cold_out.mdot);
*/ 

  TPComponents.GasStateReader r01(redeclare package Medium = medium_main) annotation(
    Placement(visible = true, transformation(origin = {66, -56}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
  TPComponents.GasStateReader r02(redeclare package Medium = medium_main) annotation(
    Placement(visible = true, transformation(origin = {96, 44}, extent = {{-6, -6}, {6, 6}}, rotation = 180)));
  TPComponents.GasStateReader r03(redeclare package Medium = medium_main) annotation(
    Placement(visible = true, transformation(origin = {-22, 50}, extent = {{-6, -6}, {6, 6}}, rotation = -90)));
  TPComponents.GasStateReader r04(redeclare package Medium = medium_main) annotation(
    Placement(visible = true, transformation(origin = {-32, 20}, extent = {{-6, -6}, {6, 6}}, rotation = 180)));
  TPComponents.GasStateReader r05(redeclare package Medium = medium_main) annotation(
    Placement(visible = true, transformation(origin = {-92, 18}, extent = {{-6, -6}, {6, 6}}, rotation = -90)));
  //Components.GasStateReader r05a(redeclare package Medium = medium_main) annotation(
  //  Placement(visible = true, transformation(origin = {-92, -32}, extent = {{-6, -6}, {6, 6}}, rotation = -90)));
  TPComponents.GasStateReader r06(redeclare package Medium = medium_main) annotation(
    Placement(visible = true, transformation(origin = {-46, -16}, extent = {{-6, -6}, {6, 6}}, rotation = 90)));
  TPComponents.GasStateReader r07(redeclare package Medium = medium_main) annotation(
    Placement(visible = true, transformation(origin = {-4, 34}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
  TPComponents.GasStateReader r08(redeclare package Medium = medium_main) annotation(
    Placement(visible = true, transformation(origin = {32, 18}, extent = {{-6, -6}, {6, 6}}, rotation = -90)));

  TPComponents.GasStateReader r08_source(redeclare package Medium = medium_main) annotation(
    Placement(visible = true, transformation(origin = {32, 18}, extent = {{-6, -6}, {6, 6}}, rotation = -90)));

  
  TPComponents.GasStateReader r08a(redeclare package Medium = medium_main) annotation(
    Placement(visible = true, transformation(origin = {20, -24}, extent = {{-6, -6}, {6, 6}}, rotation = -90)));
  TPComponents.GasStateReader r08b(redeclare package Medium = medium_main) annotation(
    Placement(visible = true, transformation(origin = {34, 4}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
  TPComponents.GasStateReader r09(redeclare package Medium = medium_main) annotation(
    Placement(visible = true, transformation(origin = {50, 62}, extent = {{-6, -6}, {6, 6}}, rotation = 180)));
  TPComponents.GasStateReader r10(redeclare package Medium = medium_main) annotation(
    Placement(visible = true, transformation(origin = {4, 58}, extent = {{-6, -6}, {6, 6}}, rotation = 180)));
  TPComponents.WaterStateReader rh1(redeclare package Medium = medium_heater) annotation(
    Placement(visible = true, transformation(origin = {-66, 22}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
  TPComponents.WaterStateReader rh2(redeclare package Medium = medium_heater) annotation(
    Placement(visible = true, transformation(origin = {-64, 64}, extent = {{-6, -6}, {6, 6}}, rotation = 180)));
  TPComponents.WaterStateReader rc1(redeclare package Medium = medium_cooler) annotation(
    Placement(visible = true, transformation(origin = {44, -36}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
  TPComponents.WaterStateReader rc2(redeclare package Medium = medium_cooler) annotation(
    Placement(visible = true, transformation(origin = {44, -74}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));

  // calculated variables
  // for turbine
  Modelica.SIunits.Power W_turb = (r05.h - r06.h) * r05.w / 1e6 "W->MW, net power for turbine";
  Modelica.SIunits.Efficiency eta_turb = turbine.eta * 100;
  // main compressor
  Modelica.SIunits.Power W_MC = (r02.h - r01.h) * r01.w / 1e6 "W->MW, net power for main compressor";
  Modelica.SIunits.Efficiency eta_MC = compressor.eta * 100;
  // re compressor
  Modelica.SIunits.Power W_RC = (r09.h - r08b.h) * r08b.w / 1e6 "W->MW, net power for recompressor";
  Modelica.SIunits.Efficiency eta_RC = recompressor.eta * 100;
  
  // power block
  Modelica.SIunits.Power W_net = W_turb - W_MC - W_RC "net power generated";
  Modelica.SIunits.Power Q_heater = (rh1.h - rh2.h) * rh1.w / 1e6 "heat input for the power block";
  Modelica.SIunits.Power Q_cooler = (rc2.h - rc1.h) * rc1.w / 1e6 "heat input for the power block";
  Modelica.SIunits.Efficiency eta_pb = W_net / Q_heater * 100 "power block efficiency";
  Real SR = r08b.w / r08.w * 100 "split ratio";
  
  // heat transfer coefficient for HTR and LTR
  // HTR
  Modelica.SIunits.Power Q_HTR = (r04.h - r03.h) * r03.w / 1e6 "W->MW, heat input for HTR";
  Modelica.SIunits.TemperatureDifference dT1_HTR = (r06.T - r04.T);
  Modelica.SIunits.TemperatureDifference dT2_HTR = (r07.T - r03.T);
  Real T_ltmd_HTR = if dT1_HTR  > 0 and dT2_HTR > 0 then (dT2_HTR - dT1_HTR) / Modelica.Math.log(abs(dT2_HTR / dT1_HTR)) else -1;
  Real UA_HTR = if T_ltmd_HTR > 0 then Q_HTR / T_ltmd_HTR else 0.0;
  // LTR
  Modelica.SIunits.Power Q_LTR = (r10.h - r02.h) * r02.w / 1e6 "W->MW";
  Modelica.SIunits.TemperatureDifference dT1_LTR = (r07.T - r10.T);
  Modelica.SIunits.TemperatureDifference dT2_LTR = (r08.T - r02.T);
  Real T_ltmd_LTR = if dT1_LTR > 0 and dT2_LTR > 0 then (dT2_LTR - dT1_LTR) / Modelica.Math.log(abs(dT2_LTR / dT1_LTR)) else -1;
  Real UA_LTR = if T_ltmd_LTR > 0 then Q_LTR / T_ltmd_LTR else 0.0;  
  
  // Liquid Na exit temperature
  Modelica.SIunits.Temperature T_heater_hot_out = rh2.T;

protected
  // performance map for main compressor
  parameter Real tableEta_comp_mc[5, 4]=[0, 95, 100, 105; 1, 0.85310219, 0.837591241, 0.832420925; 2, 0.868613139, 0.857238443, 0.851034063; 3, 0.860340633, 0.85, 0.842761557; 4, 0.85310219, 0.839659367, 0.816909976];
  
  parameter Real tablePhic_comp_mc[5, 4]=[0, 95, 100, 105; 1, 0.000134346, 0.000150832, 0.000164161; 2, 0.000137854, 0.000153638, 0.00016802; 3, 0.000142414, 0.000158549, 0.000169774; 4, 0.000145921, 0.000161706, 0.000171528];

  parameter Real tablePR_comp_mc[5, 4]=[0, 95, 100, 105; 1, 1.967529638, 2.350588505, 2.785882673; 2, 1.915294338, 2.315764972, 2.681412073; 3, 1.810823737, 2.220000255, 2.524706172; 4, 1.654117837, 2.115529655, 2.359294389];

  // performance map for re compressor      
  parameter Real tableEta_comp_rc[5, 4]=[0, 95, 100, 105; 1, 0.85310219, 0.837591241, 0.832420925; 2, 0.868613139, 0.857238443, 0.851034063; 3, 0.860340633, 0.85, 0.842761557; 4, 0.85310219, 0.839659367, 0.816909976];
  parameter Real tablePhic_comp_rc[5, 4]=[0, 95, 100, 105; 1, 7.17663E-05, 8.05731E-05, 8.76935E-05; 2, 7.36401E-05, 8.20721E-05, 8.97547E-05; 3, 7.6076E-05, 8.46954E-05, 9.06916E-05; 4, 7.79498E-05, 8.63819E-05, 9.16285E-05];

  parameter Real tablePR_comp_rc[5, 4]=[0, 95, 100, 105; 1, 1.967529638, 2.350588505, 2.785882673; 2, 1.915294338, 2.315764972, 2.681412073; 3, 1.810823737, 2.220000255, 2.524706172; 4, 1.654117837, 2.115529655, 2.359294389]; 
equation
  // open loop with compressor + recompressor
  // main path - source -> sink
//  connect(source.flange, r05a.inlet) annotation(
//    Line(points = {{-90, -22}, {-92, -22}, {-92, -27}}, color = {0, 0, 255}, thickness = 0.5));
//  connect(r05a.outlet, turbine.inlet) annotation(
//    Line(points = {{-92, -37}, {-93, -37}, {-93, -43}, {-78, -43}, {-78, -42}}, color = {0, 0, 255}, thickness = 0.5));

  connect(source.flange, r08_source.inlet);
  connect(r08_source.outlet, splitter.inlet);
  connect(splitter.outlet2, r08a.inlet) annotation(
    Line(points = {{20, 0}, {20, -19}}, color = {159, 159, 223}));
  connect(r08a.outlet, cooler.gasIn) annotation(
    Line(points = {{20, -29}, {20, -57}, {40, -57}}, color = {159, 159, 223}));
  connect(cooler.gasOut, r01.inlet) annotation(
    Line(points = {{54, -57}, {54, -56}, {62, -56}}, color = {159, 159, 223}));
  connect(r01.outlet, compressor.inlet) annotation(
    Line(points = {{70, -56}, {72, -56}, {72, 4}, {95, 4}, {95, -2}, {94, -2}}, color = {159, 159, 223}));
  connect(compressor.outlet, r02.inlet) annotation(
    Line(points = {{112, -2}, {112, 44}, {101, 44}}, color = {0, 0, 255}, thickness = 0.5));
  connect(r02.outlet, LTR.waterIn) annotation(
    Line(points = {{92, 44}, {20, 44}, {20, 40}, {22, 40}}, color = {0, 0, 255}, thickness = 0.5));
  connect(mixer.outlet, r03.inlet) annotation(
    Line(points = {{-16, 60}, {-22, 60}, {-22, 54}, {-22, 54}}, color = {0, 0, 255}, thickness = 0.5));
  connect(r03.outlet, HTR.waterIn) annotation(
    Line(points = {{-22, 46}, {-20, 46}, {-20, 40}, {-20, 40}}, color = {0, 0, 255}, thickness = 0.5));
  connect(HTR.waterOut, r04.inlet) annotation(
    Line(points = {{-20, 26}, {-21.5, 26}, {-21.5, 20}, {-27, 20}}, color = {0, 0, 255}, thickness = 0.5));
  connect(r04.outlet, heater.gasIn) annotation(
    Line(points = {{-37, 20}, {-50, 20}, {-50, 41}, {-54, 41}}, color = {0, 0, 255}, thickness = 0.5));
  connect(heater.gasOut, r05.inlet) annotation(
    Line(points = {{-68, 41}, {-92, 41}, {-92, 23}}, color = {0, 0, 255}, thickness = 0.5));
  connect(r05.outlet, turbine.inlet);
  connect(turbine.outlet, r06.inlet) annotation(
    Line(points = {{-60, -42}, {-46, -42}, {-46, -21}}, color = {159, 159, 223}));
  connect(r06.outlet, HTR.gasIn) annotation(
    Line(points = {{-46, -11}, {-45, -11}, {-45, 34}, {-28, 34}}, color = {159, 159, 223}));
  connect(HTR.gasOut, r07.inlet) annotation(
    Line(points = {{-14, 34}, {-8, 34}}, color = {159, 159, 223}));
  connect(r07.outlet, LTR.gasIn) annotation(
    Line(points = {{0, 34}, {14, 34}}, color = {159, 159, 223}));
  connect(LTR.gasOut, r08.inlet) annotation(
    Line(points = {{28, 34}, {32, 34}, {32, 22}, {32, 22}}, color = {159, 159, 223}));
  connect(r08.outlet, sink.flange) annotation(
    Line(points = {{32, 14}, {14, 14}, {14, 2}, {14, 2}}, color = {159, 159, 223}));
  
  connect(r08_source.T, sink.in_T);
  connect(r08_source.p, sink.in_p0);
//  connect(r05.outlet, sink.flange) annotation(
//    Line(points = {{-92, 13}, {-92, 0}}, color = {0, 0, 255}, thickness = 0.5));
  
//  connect(r05.T, source.in_T);
//  connect(r05.w, source.in_w0);    
  //connect(sys_init.flange, turbine.inlet);
  
  // recompressor path
  connect(splitter.outlet1, r08b.inlet) annotation(
    Line(points = {{20, 4}, {30, 4}, {30, 4}, {30, 4}}, color = {159, 159, 223}));
  connect(r08b.outlet, recompressor.inlet) annotation(
    Line(points = {{38, 4}, {48, 4}}, color = {159, 159, 223}));
  connect(recompressor.outlet, r09.inlet) annotation(
    Line(points = {{66, 4}, {66, 62}, {54, 62}}, color = {0, 0, 255}, thickness = 0.5));
  connect(r09.outlet, mixer.inlet2) annotation(
    Line(points = {{46, 62}, {-8, 62}, {-8, 62}, {-8, 62}}, color = {0, 0, 255}, thickness = 0.5));
  connect(r10.outlet, mixer.inlet1) annotation(
    Line(points = {{0, 58}, {-8, 58}}, color = {0, 0, 255}, thickness = 0.5));
  connect(LTR.waterOut, r10.inlet) annotation(
    Line(points = {{22, 26}, {21, 26}, {21, 22}, {8, 22}, {8, 58}}, color = {0, 0, 255}, thickness = 0.5));

  // turbine, compressor, recompressor
  connect(const_speed_turbine.flange, turbine.shaft_a) annotation(
    Line(points = {{-62, -51}, {-51, -51}, {-51, -50}, {-52, -50}}));
  connect(const_speed_comp.flange, compressor.shaft_a) annotation(
    Line(points = {{86, -6}, {96, -6}, {96, -10}, {96, -10}}));
  connect(const_speed_recomp.flange, recompressor.shaft_a) annotation(
    Line(points = {{38, -6}, {50, -6}, {50, -4}, {50, -4}}));

  // hot fluid path for heater
  connect(source_heater_hot.flange, rh1.inlet) annotation(
    Line(points = {{-74, 22}, {-70, 22}}, color = {255, 0, 0}, thickness = 0.5));
  connect(rh1.outlet, heater.waterIn) annotation(
    Line(points = {{-62, 22}, {-61, 22}, {-61, 34}, {-60, 34}}, color = {255, 0, 0}, thickness = 0.5));
  connect(heater.waterOut, rh2.inlet) annotation(
    Line(points = {{-60, 48}, {-60, 64}}, color = {255, 0, 0}, thickness = 0.5));
  connect(rh2.outlet, sink_heater_hot.flange) annotation(
    Line(points = {{-68, 64}, {-76, 64}}, color = {255, 0, 0}, thickness = 0.5));
    
  // cold fluid path for cooler  
  connect(source_cooler_cold.flange, rc1.inlet) annotation(
    Line(points = {{34, -36}, {40, -36}, {40, -36}, {40, -36}}, color = {0, 255, 0}, thickness = 0.5));
  connect(rc1.outlet, cooler.waterIn) annotation(
    Line(points = {{48, -36}, {48, -36}, {48, -50}, {48, -50}}, color = {0, 255, 0}, thickness = 0.5));
  connect(cooler.waterOut, rc2.inlet) annotation(
    Line(points = {{36, -74}, {40, -74}, {40, -74}, {40, -74}}, color = {0, 255, 0}, thickness = 0.5));
  connect(rc2.outlet, sink_cooler_cold.flange) annotation(
    Line(points = {{48, -74}, {48, -74}, {48, -64}, {48, -64}}, color = {0, 255, 0}, thickness = 0.5));

  annotation(
    Diagram(coordinateSystem(extent = {{-100, -100}, {120, 100}})),
    experiment(StartTime = 0, StopTime = 100, Tolerance = 1e-3, Interval = 10),
    __OpenModelica_commandLineOptions = "--matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection -d=initialization,NLSanalyticJacobian,bltdump",
    __OpenModelica_simulationFlags(lv = "LOG_DEBUG,LOG_NLS,LOG_NLS_V,LOG_STATS,LOG_INIT,LOG_STDOUT, -w", outputFormat = "mat", s = "dassl", nls = "homotopy"));
end TP_RCBCycle_PCHE;
