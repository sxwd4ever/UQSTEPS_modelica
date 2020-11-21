within Steps.Cycle;

model TP_RCBCycleGUI "Recompression Brayton Cycle with ThermoPower"
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
  // results of sscar's simulation for 10 Mw power block
  parameter Model.PBConfiguration cfg_default;
  parameter Model.PBConfiguration cfg_tune(r_i_HTR = 5e-3, r_t_HTR = cfg_default.r_i_HTR + 2e-3, r_o_HTR = cfg_default.r_t_HTR + 2e-3, N_ch_HTR = 200, L_HTR = 1 "Don't modify this, since L in HE model is fixed as 1m. Modify Nt instead", r_i_LTR = 2e-3, r_t_LTR = cfg_default.r_i_LTR + 1e-3, r_o_LTR = cfg_default.r_t_LTR + 1e-3, N_ch_LTR = 200, L_LTR = 1);
  // select the configuration of parameters
  parameter Model.PBConfiguration cfg = cfg_default;
  // set the values of parameters accordingly
  parameter HEBoundaryCondition bc_HTR = cfg.bc_HTR;
  parameter HEBoundaryCondition bc_LTR = cfg.bc_LTR;
  parameter HEBoundaryCondition bc_heater = cfg.bc_heater;
  parameter HEBoundaryCondition bc_cooler = cfg.bc_cooler;
  parameter ThermoState st_bypass = cfg.st_bypass;
  parameter EntityGeoParam geo_HTR_hot = cfg.cfg_HTR_hot.geo;
  parameter EntityGeoParam geo_HTR_cold = cfg.cfg_HTR_cold.geo;
  parameter EntityGeoParam geo_HTR_tube = cfg.cfg_HTR_tube.geo;
  parameter EntityGeoParam geo_LTR_hot = cfg.cfg_LTR_hot.geo;
  parameter EntityGeoParam geo_LTR_cold = cfg.cfg_LTR_cold.geo;
  parameter EntityGeoParam geo_LTR_tube = cfg.cfg_LTR_tube.geo;
  // use HTR's geo parameters as default
  parameter EntityGeoParam geo_heater_hot = cfg.cfg_HTR_hot.geo;
  parameter EntityGeoParam geo_heater_cold = cfg.cfg_HTR_cold.geo;
  parameter EntityGeoParam geo_heater_tube = cfg.cfg_HTR_tube.geo;
  parameter EntityGeoParam geo_cooler_hot = cfg.cfg_cooler_hot.geo;
  parameter EntityGeoParam geo_cooler_cold = cfg.cfg_cooler_cold.geo;
  parameter EntityGeoParam geo_cooler_tube = cfg.cfg_cooler_tube.geo;
  parameter EntityThermoParam thermo_HTR_hot = cfg.cfg_HTR_hot.thermo;
  parameter EntityThermoParam thermo_HTR_cold = cfg.cfg_HTR_cold.thermo;
  parameter EntityThermoParam thermo_HTR_tube = cfg.cfg_HTR_tube.thermo;
  parameter EntityThermoParam thermo_LTR_hot = cfg.cfg_LTR_hot.thermo;
  parameter EntityThermoParam thermo_LTR_cold = cfg.cfg_LTR_cold.thermo;
  parameter EntityThermoParam thermo_LTR_tube = cfg.cfg_LTR_tube.thermo;
  parameter EntityThermoParam thermo_mixer = cfg.cfg_mixer.thermo;
  parameter EntityThermoParam thermo_heater_tube = cfg.cfg_heater_tube.thermo;
  parameter EntityThermoParam thermo_cooler_tube = cfg.cfg_cooler_tube.thermo;
  //Components
  inner ThermoPower.System system(allowFlowReversal = false, initOpt = ThermoPower.Choices.Init.Options.noInit) annotation(
    Placement(visible = true, transformation(extent = {{110, 80}, {130, 100}}, rotation = 0)));
  parameter Boolean SSInit = false "Steady-state initialization";
  ThermoPower.Water.SourceMassFlow source_heater_hot(redeclare package Medium = medium_heater, w0 = bc_heater.st_hot_in.mdot, p0 = bc_heater.st_hot_in.p, h = bc_heater.st_hot_in.h, T = bc_heater.st_hot_in.T) annotation(
    Placement(visible = true, transformation(origin = {-80, 22}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
  ThermoPower.Water.SinkPressure sink_heater_hot(redeclare package Medium = medium_heater, p0 = bc_heater.st_hot_out.p, T = bc_heater.st_hot_out.T, use_T = true) annotation(
    Placement(visible = true, transformation(origin = {-81, 63}, extent = {{-5, -5}, {5, 5}}, rotation = 180)));
  ThermoPower.PowerPlants.HRSG.Components.HE heater(redeclare package FluidMedium = medium_heater, redeclare package FlueGasMedium = medium_main, fluidFlow(fixedMassFlowSimplified = true, hstartin = bc_heater.st_hot_in.h, hstartout = bc_heater.st_hot_out.h), gasFlow(redeclare package Medium = medium_main, QuasiStatic = true), N_G = geo_heater_cold.N_seg, N_F = geo_heater_hot.N_seg, Nw_G = geo_heater_tube.N_seg, gasNomFlowRate = bc_heater.st_cold_in.mdot, gasNomPressure = bc_heater.st_cold_in.p, fluidNomFlowRate = bc_heater.st_hot_in.mdot, fluidNomPressure = bc_heater.st_hot_in.p, exchSurface_G = geo_heater_cold.A_ex, exchSurface_F = geo_heater_hot.A_ex, extSurfaceTub = geo_heater_tube.A_ex, gasVol = geo_heater_cold.V, fluidVol = geo_heater_hot.V, metalVol = geo_heater_tube.V, SSInit = false, rhomcm = thermo_heater_tube.rho_mcm, lambda = thermo_heater_tube.lambda, Tstartbar_G = bc_heater.st_cold_in.T, Tstartbar_M = bc_heater.st_hot_in.T - 50, pstart_F = bc_heater.st_hot_in.p, FluidPhaseStart = ThermoPower.Choices.FluidPhase.FluidPhases.Liquid, redeclare replaceable model HeatTransfer_F = ThermoPower.Thermal.HeatTransferFV.IdealHeatTransfer, redeclare replaceable model HeatTransfer_G = ThermoPower.Thermal.HeatTransferFV.IdealHeatTransfer, redeclare model HeatExchangerTopology = ThermoPower.Thermal.HeatExchangerTopologies.CounterCurrentFlow, metalTube(WallRes = false, Tstartbar = heater.Tstartbar_M)) annotation(
    Placement(visible = true, transformation(origin = {-61, 41}, extent = {{-7, -7}, {7, 7}}, rotation = 180)));
  // set the fluid flow as fixed mdot for simplarity
  // ConstantHeatTransferCoefficientTwoGrids(gamma = thermo_HTR_hot.gamma_he),
  // ConstantHeatTransferCoefficient(gamma = thermo_HTR_cold.gamma_he),
  // set the fluid flow as fixed mdot for simplarity
  // ConstantHeatTransferCoefficientTwoGrids(gamma = thermo_HTR_hot.gamma_he),
  // ConstantHeatTransferCoefficient(gamma = thermo_HTR_cold.gamma_he),
  // use FlowJoin to mix flow
  // use FlowJoin to mix flow
  ThermoPower.Gas.FlowJoin mixer(redeclare package Medium = medium_main) annotation(
    Placement(visible = true, transformation(origin = {-12, 60}, extent = {{-6, -6}, {6, 6}}, rotation = 180)));
  /*
                                          ThermoPower.Gas.SensT T_waterOut(
                                            redeclare package Medium = medium_main) 
                                          annotation(
                                            Placement(transformation(origin = {4, -50}, extent = {{-10, -10}, {10, 10}}, rotation = 270)));
                                          
                                          ThermoPower.Gas.SensT T_gasOut(
                                            redeclare package Medium = medium_main) 
                                          annotation(
                                            Placement(transformation(extent = {{30, -6}, {50, 14}}, rotation = 0)));
                                          */
  Components.HEG2G HTR(redeclare package FluidMedium = medium_main, redeclare package FlueGasMedium = medium_main, redeclare replaceable model HeatTransfer_F = ThermoPower.Thermal.HeatTransferFV.ConstantHeatTransferCoefficient(gamma = thermo_HTR_cold.gamma_he), redeclare replaceable model HeatTransfer_G = ThermoPower.Thermal.HeatTransferFV.ConstantHeatTransferCoefficientTwoGrids(gamma = thermo_HTR_hot.gamma_he), redeclare model HeatExchangerTopology = ThermoPower.Thermal.HeatExchangerTopologies.CounterCurrentFlow, N_F = geo_HTR_cold.N_seg, N_G = geo_HTR_hot.N_seg, Tstartbar_G = bc_HTR.st_hot_in.T, Tstartbar_F = bc_HTR.st_cold_in.T, exchSurface_F = geo_HTR_cold.A_ex, exchSurface_G = geo_HTR_hot.A_ex, extSurfaceTub = geo_HTR_tube.A_ex, fluidNomFlowRate = bc_HTR.st_cold_in.mdot, fluidNomPressure = bc_HTR.st_cold_in.p, fluidVol = geo_HTR_cold.V, gasNomFlowRate = bc_HTR.st_hot_in.mdot, gasNomPressure = bc_HTR.st_hot_in.p, gasVol = geo_HTR_hot.V, lambda = thermo_HTR_tube.lambda, metalVol = geo_HTR_tube.V, pstart_F = bc_HTR.st_cold_in.p, pstart_G = bc_HTR.st_hot_in.p, rhomcm = thermo_HTR_tube.rho_mcm, gasQuasiStatic = true, fluidQuasiStatic = true, metalTube(WallRes = false, Tstartbar = HTR.Tstartbar_M)) annotation(
    Placement(visible = true, transformation(origin = {-21, 33}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
  // redeclare replaceable model HeatTransfer_F = ThermoPower.Thermal.HeatTransferFV.IdealHeatTransfer,
  //redeclare replaceable model HeatTransfer_G = ThermoPower.Thermal.HeatTransferFV.IdealHeatTransfer,
  // redeclare replaceable model HeatTransfer_F = ThermoPower.Thermal.HeatTransferFV.IdealHeatTransfer,
  //redeclare replaceable model HeatTransfer_G = ThermoPower.Thermal.HeatTransferFV.IdealHeatTransfer,
  Components.HEG2G LTR(redeclare package FluidMedium = medium_main, redeclare package FlueGasMedium = medium_main, redeclare replaceable model HeatTransfer_F = ThermoPower.Thermal.HeatTransferFV.ConstantHeatTransferCoefficient(gamma = thermo_LTR_cold.gamma_he), redeclare replaceable model HeatTransfer_G = ThermoPower.Thermal.HeatTransferFV.ConstantHeatTransferCoefficientTwoGrids(gamma = thermo_LTR_hot.gamma_he), redeclare model HeatExchangerTopology = ThermoPower.Thermal.HeatExchangerTopologies.CounterCurrentFlow, N_F = geo_LTR_cold.N_seg, N_G = geo_LTR_hot.N_seg, Tstartbar_G = bc_LTR.st_hot_in.T, Tstartbar_F = bc_LTR.st_cold_in.T, exchSurface_F = geo_LTR_cold.A_ex, exchSurface_G = geo_LTR_hot.A_ex, extSurfaceTub = geo_LTR_tube.A_ex, fluidNomFlowRate = bc_LTR.st_cold_in.mdot, fluidNomPressure = bc_LTR.st_cold_in.p, fluidVol = geo_LTR_cold.V, gasNomFlowRate = bc_LTR.st_hot_in.mdot, gasNomPressure = bc_LTR.st_hot_in.p, gasVol = geo_LTR_hot.V, lambda = thermo_LTR_tube.lambda, metalVol = geo_LTR_tube.V, pstart_F = bc_LTR.st_cold_in.p, pstart_G = bc_LTR.st_hot_in.p, rhomcm = thermo_LTR_tube.rho_mcm, gasQuasiStatic = true, fluidQuasiStatic = true, metalTube(WallRes = false, Tstartbar = LTR.Tstartbar_M)) annotation(
    Placement(visible = true, transformation(origin = {21, 33}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
  // redeclare replaceable model HeatTransfer_F = ThermoPower.Thermal.HeatTransferFV.IdealHeatTransfer,
  // redeclare replaceable model HeatTransfer_G = ThermoPower.Thermal.HeatTransferFV.IdealHeatTransfer,
  //SSInit = SSInit,
  ThermoPower.Gas.Turbine turbine(redeclare package Medium = medium_main, fileName = Modelica.Utilities.Files.loadResource("modelica://Steps/Resources/Data/turbine_map.txt"), tablePhic = fill(0.0, 14, 12), tableEta = fill(0.0, 14, 12), pstart_in = bc_heater.st_cold_out.p, pstart_out = bc_HTR.st_hot_in.p, Tstart_in = bc_heater.st_cold_out.T, Tstart_out = bc_HTR.st_hot_in.T, Ndesign = 60000.0, Tdes_in = bc_heater.st_cold_out.T, Table = ThermoPower.Choices.TurboMachinery.TableTypes.file, explicitIsentropicEnthalpy = false, gas_in(p(nominal = turbine.pstart_in), T(nominal = turbine.Tstart_in)), gas_iso(p(nominal = turbine.pstart_out), T(nominal = turbine.Tstart_out))) annotation(
    Placement(visible = true, transformation(origin = {-69, -51}, extent = {{-11, -11}, {11, 11}}, rotation = 0)));
  //tablePhic,
  //tableEta,
  Modelica.Mechanics.Rotational.Sources.Speed speed1 annotation(
    Placement(visible = true, transformation(origin = {-46, -50}, extent = {{6, -6}, {-6, 6}}, rotation = 0)));
  Modelica.Blocks.Sources.Constant const1(k = 60000) annotation(
    Placement(visible = true, transformation(origin = {-26, -50}, extent = {{6, -6}, {-6, 6}}, rotation = 0)));
  ThermoPower.Water.SourceMassFlow source_cooler_cold(redeclare package Medium = medium_cooler, T = bc_cooler.st_cold_in.T, p0 = bc_cooler.st_cold_in.p, use_T = true, use_in_T = false, w0 = bc_cooler.st_cold_in.mdot) annotation(
    Placement(visible = true, transformation(origin = {29, -37}, extent = {{-5, -5}, {5, 5}}, rotation = 0)));
  ThermoPower.Water.SinkPressure sink_cooler_cold(redeclare package Medium = medium_cooler, p0 = bc_cooler.st_cold_out.p, T = bc_cooler.st_cold_out.T, use_T = true) annotation(
    Placement(visible = true, transformation(origin = {31, -75}, extent = {{-5, -5}, {5, 5}}, rotation = 180)));
  ThermoPower.PowerPlants.HRSG.Components.HE cooler(redeclare package FluidMedium = medium_cooler, redeclare package FlueGasMedium = medium_main, fluidFlow(fixedMassFlowSimplified = true, hstartin = bc_cooler.st_cold_in.h, hstartout = bc_cooler.st_cold_out.h), gasFlow(QuasiStatic = true), N_G = geo_cooler_hot.N_seg, N_F = geo_cooler_cold.N_seg, Nw_G = geo_cooler_tube.N_seg, gasNomFlowRate = bc_cooler.st_hot_in.mdot, gasNomPressure = bc_cooler.st_hot_in.p, fluidNomFlowRate = bc_cooler.st_cold_in.mdot, fluidNomPressure = bc_cooler.st_cold_in.p, exchSurface_G = geo_cooler_hot.A_ex, exchSurface_F = geo_cooler_cold.A_ex, extSurfaceTub = geo_cooler_tube.A_ex, gasVol = geo_cooler_hot.V, fluidVol = geo_cooler_cold.V, metalVol = geo_cooler_tube.V, SSInit = false, rhomcm = thermo_cooler_tube.rho_mcm "use thermo props of heater for simplicity", lambda = thermo_cooler_tube.lambda "use thermo props of heater for simplicity", Tstartbar_G = bc_cooler.st_hot_in.T, Tstartbar_M = bc_cooler.st_hot_in.T - 50, pstart_F = bc_cooler.st_cold_in.p, FluidPhaseStart = ThermoPower.Choices.FluidPhase.FluidPhases.Liquid, redeclare replaceable model HeatTransfer_F = ThermoPower.Thermal.HeatTransferFV.IdealHeatTransfer, redeclare replaceable model HeatTransfer_G = ThermoPower.Thermal.HeatTransferFV.IdealHeatTransfer, redeclare model HeatExchangerTopology = ThermoPower.Thermal.HeatExchangerTopologies.CounterCurrentFlow, metalTube(WallRes = false, Tstartbar = bc_cooler.st_hot_in.T - 50)) annotation(
    Placement(visible = true, transformation(origin = {47, -57}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
  //Components.HEG2G cooler(
  // set the fluid flow as fixed mdot for simplarity
  //ConstantHeatTransferCoefficientTwoGrids(gamma = thermo_hot.gamma_he),
  //ConstantHeatTransferCoefficient(gamma =  thermo_cold.gamma_he),
  //Components.HEG2G cooler(
  // set the fluid flow as fixed mdot for simplarity
  //ConstantHeatTransferCoefficientTwoGrids(gamma = thermo_hot.gamma_he),
  //ConstantHeatTransferCoefficient(gamma =  thermo_cold.gamma_he),
  ThermoPower.Gas.FlowSplit splitter(redeclare package Medium = medium_main) annotation(
    Placement(visible = true, transformation(origin = {17, 1}, extent = {{-5, -5}, {5, 5}}, rotation = 0)));
  ThermoPower.Gas.Compressor compressor(redeclare package Medium = medium_main, pstart_in = bc_cooler.st_hot_out.p, pstart_out = bc_LTR.st_cold_in.p, Tstart_in = bc_cooler.st_hot_out.T, Tstart_out = bc_LTR.st_cold_in.T, tablePhic = tablePhic_comp_mc, tableEta = tableEta_comp_mc, tablePR = tablePR_comp_mc, Table = ThermoPower.Choices.TurboMachinery.TableTypes.matrix, Ndesign = 523.3, Tdes_in = compressor.Tstart_in) annotation(
    Placement(visible = true, transformation(origin = {103, -11}, extent = {{-11, -11}, {11, 11}}, rotation = 0)));
  Modelica.Mechanics.Rotational.Sources.ConstantSpeed ConstantSpeed1(w_fixed = 523.3, useSupport = false) annotation(
    Placement(visible = true, transformation(origin = {81, -7}, extent = {{-5, -5}, {5, 5}}, rotation = 0)));
  ThermoPower.Gas.Compressor compressor_bypass(redeclare package Medium = medium_main, pstart_in = bc_HTR.st_hot_out.p, pstart_out = st_bypass.p, Tstart_in = bc_HTR.st_hot_out.T, Tstart_out = st_bypass.T, tablePhic = tablePhic_comp_rc, tableEta = tableEta_comp_rc, tablePR = tablePR_comp_rc, Table = ThermoPower.Choices.TurboMachinery.TableTypes.matrix, Ndesign = 523.3, Tdes_in = compressor_bypass.Tstart_in) annotation(
    Placement(visible = true, transformation(origin = {57, -5}, extent = {{-11, -11}, {11, 11}}, rotation = 0)));
  Modelica.Mechanics.Rotational.Sources.ConstantSpeed ConstantSpeed2(w_fixed = 523.3, useSupport = false) annotation(
    Placement(visible = true, transformation(origin = {33, -7}, extent = {{-5, -5}, {5, 5}}, rotation = 0)));
  ThermoPower.Gas.SourceMassFlow source(redeclare package Medium = medium_main, w0 = bc_heater.st_cold_out.mdot, p0 = bc_heater.st_cold_out.p, T = bc_heater.st_cold_out.T, gas(T(nominal = source.T), p(nominal = source.p0))) annotation(
    Placement(visible = true, transformation(origin = {-91, -17}, extent = {{-5, -5}, {5, 5}}, rotation = -90)));
  ThermoPower.Gas.SinkPressure sink(redeclare package Medium = medium_main, p0 = bc_heater.st_cold_out.p, T = bc_heater.st_cold_out.T) annotation(
    Placement(visible = true, transformation(origin = {-92, -4}, extent = {{-4, -4}, {4, 4}}, rotation = 270)));
  /*
                                          ThermoPower.Gas.SensT sensT_sink(redeclare package Medium = medium_main) annotation(
                                            Placement(visible = true, transformation(origin = {-91, 23}, extent = {{-5, -5}, {5, 5}}, rotation = -90)));
                                          ThermoPower.Gas.SensW sensW_sink(redeclare package Medium = medium_main) annotation(
                                            Placement(visible = true, transformation(origin = {-92, 10}, extent = {{-6, -6}, {6, 6}}, rotation = -90)));
                                        */
  /*
                                          ThermoPower.Gas.ThroughMassFlow source_through_cooler(
                                            redeclare package Medium = medium_main,
                                            w0 = bc_cooler.st_hot_out.mdot) annotation(
                                            Placement(visible = true, transformation(origin = {88, 12}, extent = {{-6, -6}, {6, 6}}, rotation = 90)));

                                          
                                          ThermoPower.Gas.ThroughMassFlow source_through_heater(
                                            redeclare package Medium = medium_main,
                                            w0 = bc_heater.st_cold_out.mdot);     
                                             

                                          ThermoPower.Gas.Utility.ClosedSystemInit sys_init(
                                            redeclare package Medium = medium_main,
                                            //w_b(nominal = bc_heater.st_cold_out.mdot),
                                            pstart = bc_heater.st_cold_out.p);
                                        */
  /*
                                          ThermoPower.Gas.SourceMassFlow source_bypass(
                                          redeclare package Medium = medium_main,
                                            w0 = st_bypass.mdot,
                                            p0 = st_bypass.p,
                                            T = st_bypass.T) 
                                            annotation(
                                            Placement(transformation(origin = {0, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 270)));
                                            
                                          ThermoPower.Gas.SinkPressure sink_bypass(
                                          redeclare package Medium = medium_main,
                                            p0 = bc_LTR.st_hot_out.p,
                                            T = bc_LTR.st_hot_out.T) 
                                            annotation(
                                            Placement(transformation(origin = {0, -80}, extent = {{-10, -10}, {10, 10}}, rotation = 270)));
                                          

                                          ThermoPower.Gas.ThroughMassFlow source_through_bypass(
                                            redeclare package Medium = medium_main,
                                            w0 = st_bypass.mdot)
                                            annotation(
                                            Placement(visible = true, transformation(origin = {61, 29}, extent = {{-5, -5}, {5, 5}}, rotation = 90)));    
                                                      */
  ThermoPower.PowerPlants.HRSG.Components.StateReader_gas r01(redeclare package Medium = medium_main) annotation(
    Placement(visible = true, transformation(origin = {66, -56}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
  ThermoPower.PowerPlants.HRSG.Components.StateReader_gas r02(redeclare package Medium = medium_main) annotation(
    Placement(visible = true, transformation(origin = {96, 44}, extent = {{-6, -6}, {6, 6}}, rotation = 180)));
  ThermoPower.PowerPlants.HRSG.Components.StateReader_gas r03(redeclare package Medium = medium_main) annotation(
    Placement(visible = true, transformation(origin = {-22, 50}, extent = {{-6, -6}, {6, 6}}, rotation = -90)));
  ThermoPower.PowerPlants.HRSG.Components.StateReader_gas r04(redeclare package Medium = medium_main) annotation(
    Placement(visible = true, transformation(origin = {-32, 20}, extent = {{-6, -6}, {6, 6}}, rotation = 180)));
  ThermoPower.PowerPlants.HRSG.Components.StateReader_gas r05(redeclare package Medium = medium_main) annotation(
    Placement(visible = true, transformation(origin = {-92, 18}, extent = {{-6, -6}, {6, 6}}, rotation = -90)));
  ThermoPower.PowerPlants.HRSG.Components.StateReader_gas r05a(redeclare package Medium = medium_main) annotation(
    Placement(visible = true, transformation(origin = {-92, -32}, extent = {{-6, -6}, {6, 6}}, rotation = -90)));
  ThermoPower.PowerPlants.HRSG.Components.StateReader_gas r06(redeclare package Medium = medium_main) annotation(
    Placement(visible = true, transformation(origin = {-46, -16}, extent = {{-6, -6}, {6, 6}}, rotation = 90)));
  ThermoPower.PowerPlants.HRSG.Components.StateReader_gas r07(redeclare package Medium = medium_main) annotation(
    Placement(visible = true, transformation(origin = {-4, 34}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
  ThermoPower.PowerPlants.HRSG.Components.StateReader_gas r08(redeclare package Medium = medium_main) annotation(
    Placement(visible = true, transformation(origin = {32, 18}, extent = {{-6, -6}, {6, 6}}, rotation = -90)));
  ThermoPower.PowerPlants.HRSG.Components.StateReader_gas r08a(redeclare package Medium = medium_main) annotation(
    Placement(visible = true, transformation(origin = {20, -24}, extent = {{-6, -6}, {6, 6}}, rotation = -90)));
  ThermoPower.PowerPlants.HRSG.Components.StateReader_gas r08b(redeclare package Medium = medium_main) annotation(
    Placement(visible = true, transformation(origin = {34, 4}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
  ThermoPower.PowerPlants.HRSG.Components.StateReader_gas r09(redeclare package Medium = medium_main) annotation(
    Placement(visible = true, transformation(origin = {50, 62}, extent = {{-6, -6}, {6, 6}}, rotation = 180)));
  ThermoPower.PowerPlants.HRSG.Components.StateReader_gas r10(redeclare package Medium = medium_main) annotation(
    Placement(visible = true, transformation(origin = {4, 58}, extent = {{-6, -6}, {6, 6}}, rotation = 180)));
  ThermoPower.PowerPlants.HRSG.Components.StateReader_water rh1(redeclare package Medium = medium_heater) annotation(
    Placement(visible = true, transformation(origin = {-66, 22}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
  ThermoPower.PowerPlants.HRSG.Components.StateReader_water rh2(redeclare package Medium = medium_heater) annotation(
    Placement(visible = true, transformation(origin = {-64, 64}, extent = {{-6, -6}, {6, 6}}, rotation = 180)));
  ThermoPower.PowerPlants.HRSG.Components.StateReader_water rc1(redeclare package Medium = medium_cooler) annotation(
    Placement(visible = true, transformation(origin = {44, -36}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
  ThermoPower.PowerPlants.HRSG.Components.StateReader_water rc2(redeclare package Medium = medium_cooler) annotation(
    Placement(visible = true, transformation(origin = {44, -74}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
protected
  // performance map for main compressor
  parameter Real tableEta_comp_mc[5, 4] = [0, 95, 100, 105; 1, 0.85310219, 0.837591241, 0.832420925; 2, 0.868613139, 0.857238443, 0.851034063; 3, 0.860340633, 0.85, 0.842761557; 4, 0.85310219, 0.839659367, 0.816909976];
  parameter Real tablePhic_comp_mc[5, 4] = [0, 95, 100, 105; 1, 0.000134346, 0.000150832, 0.000164161; 2, 0.000137854, 0.000153638, 0.00016802; 3, 0.000142414, 0.000158549, 0.000169774; 4, 0.000145921, 0.000161706, 0.000171528];
  parameter Real tablePR_comp_mc[5, 4] = [0, 95, 100, 105; 1, 1.967529638, 2.350588505, 2.785882673; 2, 1.915294338, 2.315764972, 2.681412073; 3, 1.810823737, 2.220000255, 2.524706172; 4, 1.654117837, 2.115529655, 2.359294389];
  // performance map for re compressor
  parameter Real tableEta_comp_rc[5, 4] = [0, 95, 100, 105; 1, 0.85310219, 0.837591241, 0.832420925; 2, 0.868613139, 0.857238443, 0.851034063; 3, 0.860340633, 0.85, 0.842761557; 4, 0.85310219, 0.839659367, 0.816909976];
  parameter Real tablePhic_comp_rc[5, 4] = [0, 95, 100, 105; 1, 7.17663E-05, 8.05731E-05, 8.76935E-05; 2, 7.36401E-05, 8.20721E-05, 8.97547E-05; 3, 7.6076E-05, 8.46954E-05, 9.06916E-05; 4, 7.79498E-05, 8.63819E-05, 9.16285E-05];
  parameter Real tablePR_comp_rc[5, 4] = [0, 95, 100, 105; 1, 1.967529638, 2.350588505, 2.785882673; 2, 1.915294338, 2.315764972, 2.681412073; 3, 1.810823737, 2.220000255, 2.524706172; 4, 1.654117837, 2.115529655, 2.359294389];
equation
/*  
  // Close Loop with a through mass flow source 
  connect(source.flange, splitter.inlet);
  connect(splitter.outlet1, cooler.gasIn);
  connect(cooler.gasOut, compressor.inlet);
  connect(compressor.outlet, LTR.waterIn);
  connect(compressor.shaft_a, ConstantSpeed1.flange);
  
  connect(LTR.waterOut, mixer.inlet1);

  connect(splitter.outlet2, compressor_bypass.inlet);
  connect(compressor_bypass.outlet, mixer.inlet2);  
  connect(compressor_bypass.shaft_a, ConstantSpeed2.flange);
  
  connect(mixer.outlet, HTR.waterIn);
  connect(HTR.waterOut, heater.gasIn);  
  connect(heater.gasOut, turbine.inlet) annotation(
    Line(points = {{1.83697e-015, -70}, {1.83697e-015, -56}, {-8.88178e-016, -56}}, thickness = 0.5, color = {0, 0, 255}));    

  connect(turbine.shaft_b, speed1.flange) annotation(
    Line(points = {{30, 0}, {74, 0}, {74, 0}, {74, 0}}));
  connect(speed1.w_ref, const1.y) annotation(
    Line(points = {{96, 0}, {120, 0}, {120, 0}, {118, 0}}, color = {0, 0, 127}));

  connect(turbine.outlet, sens_turbine.inlet) annotation(
   Line(points = {{-50, 0}, {-20, 0}}, color = {159, 159, 223}, thickness = 0.5, smooth = Smooth.None));  
  connect(sens_turbine.outlet, HTR.gasIn);  
  connect(HTR.gasOut, LTR.gasIn);  
  connect(LTR.gasOut, sink.flange) annotation(
    Line(points = {{46, 0}, {46, 0}, {60, 0}}, color = {159, 159, 223}, thickness = 0.5));
   
  // hot stream for heater
  connect(source_heater_hot.flange, heater.waterIn);
  connect(heater.waterOut, sink_heater_hot.flange);

  // cold stream for cooler
  connect(source_cooler_cold.flange, cooler.waterIn);
  connect(cooler.waterOut, sink_cooler_cold.flange);
*/
/*
  // open loop without compressor
  connect(source.flange, LTR.waterIn) annotation(
    Line(points = {{83, 36}, {83, 61}, {21, 61}, {21, 40}}, color = {0, 0, 255}, thickness = 0.5));
     
  connect(LTR.waterOut, mixer.inlet1) annotation(
    Line(points = {{21, 26}, {21.25, 26}, {21.25, 20}, {-0.5, 20}, {-0.5, 48}, {-8, 48}}, color = {0, 0, 255}, thickness = 0.5));
  connect(source_through_bypass.outlet, mixer.inlet2) annotation(
    Line(points = {{61, 34}, {61, 52}, {-8, 52}}, color = {0, 0, 255}, thickness = 0.5));

  connect(mixer.outlet, HTR.waterIn) annotation(
    Line(points = {{-16, 50}, {-21, 50}, {-21, 40}}, color = {0, 0, 255}, thickness = 0.5));
    
  connect(HTR.waterOut, heater.gasIn) annotation(
    Line(points = {{-21, 26}, {-27.5, 26}, {-27.5, 20}, {-50, 20}, {-50, 34}, {-62, 34}}, color = {0, 0, 255}, thickness = 0.5));  
  connect(heater.gasOut, turbine.inlet) annotation(
    Line(points = {{-76, 33}, {-76, 33.75}, {-92, 33.75}, {-92, -41.5}, {-78, -41.5}, {-78, -42}}, color = {0, 0, 255}, thickness = 0.5));
    
  connect(turbine.shaft_b, speed1.flange) annotation(
    Line(points = {{-62, -51}, {-51, -51}, {-51, -50}, {-52, -50}}));
    
  connect(speed1.w_ref, const1.y) annotation(
    Line(points = {{-39, -50}, {-33, -50}}, color = {0, 0, 127}));
    
  connect(turbine.outlet, sens_turbine.inlet) annotation(
    Line(points = {{-60, -42}, {-42, -42}, {-42, -14}, {-42, -14}}, color = {159, 159, 223}));
  connect(sens_turbine.outlet, HTR.gasIn) annotation(
    Line(points = {{-42, -6}, {-42, 34}, {-28, 34}}, color = {159, 159, 223}));  
  connect(HTR.gasOut, LTR.gasIn) annotation(
    Line(points = {{-14, 33}, {14, 33}}, color = {159, 159, 223})); 
  connect(LTR.gasOut, splitter.inlet) annotation(
    Line(points = {{36, 10}, {36, 34}, {28, 34}}, color = {159, 159, 223}));
  
  connect(splitter.outlet1, cooler.gasIn) annotation(
    Line(points = {{37, 4}, {61, 4}, {61, 24}}, color = {159, 159, 223}));
  connect(cooler.gasOut, sink.flange) annotation(
    Line(points = {{64, -47}, {82, -47}, {82, 22}}, color = {159, 159, 223})); 
  //connect(sys_init.flange, mixer.inlet1);    
  
  connect(splitter.outlet2, source_through_bypass.inlet) annotation(
    Line(points = {{33, 4}, {33, -47.5}, {50, -47.5}, {50, -47}}, color = {159, 159, 223}));
   
  // hot stream for heater
  connect(source_heater_hot.flange, heater.waterIn) annotation(
    Line(points = {{-70, 14}, {-68.5, 14}, {-68.5, 26}, {-69, 26}}, color = {0, 0, 255}));
  connect(heater.waterOut, sink_heater_hot.flange) annotation(
    Line(points = {{-69, 40}, {-69, 57}, {-72, 57}}, color = {0, 0, 255}));

  // cold stream for cooler
  connect(source_cooler_cold.flange, cooler.waterIn) annotation(
    Line(points = {{56, -26}, {60, -26}, {60, -40}, {58, -40}}, color = {0, 0, 255}));
  connect(cooler.waterOut, sink_cooler_cold.flange) annotation(
    Line(points = {{57, -54}, {57, -66}, {58, -66}, {58, -67}}, color = {0, 0, 255}));
*/
// close loop without compressor
/*
//connect(ConstantSpeed1.flange, compressor.shaft_a);
  //connect(source_through_heater.outlet, turbine.inlet);
  connect(heater.gasOut, turbine.inlet);
  connect(sys_init.flange, turbine.inlet);
  */
  connect(turbine.shaft_a, speed1.flange) annotation(
    Line(points = {{-62, -51}, {-51, -51}, {-51, -50}, {-52, -50}}));
  connect(speed1.w_ref, const1.y) annotation(
    Line(points = {{-39, -50}, {-33, -50}}, color = {0, 0, 127}));
/*
  connect(LTR.gasOut, splitter.inlet) annotation(
    Line(points = {{36, 10}, {36, 14}, {35, 14}}, color = {159, 159, 223}));
  */
//connect(sys_init.flange, mixer.inlet1);
// hot stream for heater
  connect(compressor_bypass.shaft_a, ConstantSpeed2.flange) annotation(
    Line);
  connect(compressor.shaft_a, ConstantSpeed1.flange) annotation(
    Line);
/*
  connect(heater.gasOut, sensT_sink.inlet) annotation(
    Line(points = {{-72, 41}, {-92, 41}, {-92, 26}, {-93, 26}}, color = {159, 159, 223}));
  connect(sensT_sink.outlet, sensW_sink.inlet) annotation(
    Line(points = {{-93, 20}, {-94, 20}, {-94, 14}}, color = {159, 159, 223}));
  */
/*
  connect(sensW_sink.w, source.in_w0) annotation(
    Line(points = {{-88, 6}, {-86.5, 6}, {-86.5, -22}}, color = {0, 0, 127}));
  connect(sensT_sink.T, source.in_T) annotation(
    Line(points = {{-88, 19.5}, {-88, 19.125}, {-82, 19.125}, {-82, -24.625}, {-86, -24.625}, {-86, -24}}, color = {0, 0, 127}));   
  */
// connect(heater.gasOut, source_through_heater.inlet);
// cold stream for cooler
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
  connect(r05.outlet, sink.flange) annotation(
    Line(points = {{-92, 13}, {-92, 0}}, color = {0, 0, 255}, thickness = 0.5));
  connect(source.flange, r05a.inlet) annotation(
    Line(points = {{-90, -22}, {-92, -22}, {-92, -27}}, color = {0, 0, 255}, thickness = 0.5));
  connect(r05a.outlet, turbine.inlet) annotation(
    Line(points = {{-92, -37}, {-93, -37}, {-93, -43}, {-78, -43}, {-78, -42}}, color = {0, 0, 255}, thickness = 0.5));
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
  connect(r08.outlet, splitter.inlet) annotation(
    Line(points = {{32, 14}, {14, 14}, {14, 2}, {14, 2}}, color = {159, 159, 223}));
  connect(splitter.outlet2, r08a.inlet) annotation(
    Line(points = {{20, 0}, {20, -19}}, color = {159, 159, 223}));
  connect(r08a.outlet, cooler.gasIn) annotation(
    Line(points = {{20, -29}, {20, -57}, {40, -57}}, color = {159, 159, 223}));
  connect(splitter.outlet1, r08b.inlet) annotation(
    Line(points = {{20, 4}, {30, 4}, {30, 4}, {30, 4}}, color = {159, 159, 223}));
  connect(r08b.outlet, compressor_bypass.inlet) annotation(
    Line(points = {{38, 4}, {48, 4}}, color = {159, 159, 223}));
  connect(source_heater_hot.flange, rh1.inlet) annotation(
    Line(points = {{-74, 22}, {-70, 22}}, color = {255, 0, 0}, thickness = 0.5));
  connect(rh1.outlet, heater.waterIn) annotation(
    Line(points = {{-62, 22}, {-61, 22}, {-61, 34}, {-60, 34}}, color = {255, 0, 0}, thickness = 0.5));
  connect(heater.waterOut, rh2.inlet) annotation(
    Line(points = {{-60, 48}, {-60, 64}}, color = {255, 0, 0}, thickness = 0.5));
  connect(rh2.outlet, sink_heater_hot.flange) annotation(
    Line(points = {{-68, 64}, {-76, 64}}, color = {255, 0, 0}, thickness = 0.5));
  connect(source_cooler_cold.flange, rc1.inlet) annotation(
    Line(points = {{34, -36}, {40, -36}, {40, -36}, {40, -36}}, color = {0, 255, 0}, thickness = 0.5));
  connect(rc1.outlet, cooler.waterIn) annotation(
    Line(points = {{48, -36}, {48, -36}, {48, -50}, {48, -50}}, color = {0, 255, 0}, thickness = 0.5));
  connect(sink_cooler_cold.flange, rc2.inlet) annotation(
    Line(points = {{36, -74}, {40, -74}, {40, -74}, {40, -74}}, color = {0, 255, 0}, thickness = 0.5));
  connect(rc2.outlet, cooler.waterOut) annotation(
    Line(points = {{48, -74}, {48, -74}, {48, -64}, {48, -64}}, color = {0, 255, 0}, thickness = 0.5));
  connect(compressor_bypass.outlet, r09.inlet) annotation(
    Line(points = {{66, 4}, {66, 62}, {54, 62}}, color = {0, 0, 255}, thickness = 0.5));
  connect(r09.outlet, mixer.inlet2) annotation(
    Line(points = {{46, 62}, {-8, 62}, {-8, 62}, {-8, 62}}, color = {0, 0, 255}, thickness = 0.5));
  connect(r10.outlet, mixer.inlet1) annotation(
    Line(points = {{0, 58}, {-8, 58}}, color = {0, 0, 255}, thickness = 0.5));
  connect(LTR.waterOut, r10.inlet) annotation(
    Line(points = {{22, 26}, {21, 26}, {21, 22}, {8, 22}, {8, 58}}, color = {0, 0, 255}, thickness = 0.5));
  connect(ConstantSpeed1.flange, compressor.shaft_a) annotation(
    Line(points = {{86, -6}, {96, -6}, {96, -10}, {96, -10}}));
  connect(ConstantSpeed2.flange, compressor_bypass.shaft_a) annotation(
    Line(points = {{38, -6}, {50, -6}, {50, -4}, {50, -4}}));
  annotation(
    Diagram(coordinateSystem(extent = {{-100, -100}, {130, 100}})),
    experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-3, Interval = 1),
    __OpenModelica_commandLineOptions = "--matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection -d=initialization,NLSanalyticJacobian,bltdump",
    __OpenModelica_simulationFlags(lv = "LOG_DEBUG,LOG_NLS,LOG_NLS_V,LOG_STATS,LOG_INIT,LOG_STDOUT, -w", outputFormat = "mat", s = "dassl", nls = "homotopy"),
    Icon(coordinateSystem(extent = {{-100, -100}, {130, 100}})));
end TP_RCBCycleGUI;
