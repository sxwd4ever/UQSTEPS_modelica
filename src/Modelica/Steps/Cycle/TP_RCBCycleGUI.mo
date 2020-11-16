within Steps.Cycle;

model TP_RCBCycleGUI "Recompression Brayton Cycle with ThermoPower for layout"
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
    Placement(transformation(extent = {{80, 80}, {100, 100}})));
  parameter Boolean SSInit = false "Steady-state initialization";
  ThermoPower.Water.SourceMassFlow source_heater_hot(redeclare package Medium = medium_heater, w0 = bc_heater.st_hot_in.mdot, p0 = bc_heater.st_hot_in.p, T = bc_heater.st_hot_in.T, use_T = true) annotation(
    Placement(visible = true, transformation(origin = {-76, 14}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
  ThermoPower.Water.SinkPressure sink_heater_hot(redeclare package Medium = medium_heater, p0 = bc_heater.st_hot_out.mdot, T = bc_heater.st_hot_out.T, use_T = true) annotation(
    Placement(visible = true, transformation(origin = {-79, 49}, extent = {{-5, -5}, {5, 5}}, rotation = 180)));
  ThermoPower.PowerPlants.HRSG.Components.HE heater(redeclare package FluidMedium = medium_heater, redeclare package FlueGasMedium = medium_main, fluidFlow(fixedMassFlowSimplified = true, hstartin = bc_heater.st_hot_in.h, hstartout = bc_heater.st_hot_out.h), gasFlow(QuasiStatic = true), N_G = geo_heater_cold.N_seg, N_F = geo_heater_hot.N_seg, Nw_G = geo_heater_tube.N_seg, gasNomFlowRate = bc_heater.st_cold_in.mdot, gasNomPressure = bc_heater.st_cold_in.p, fluidNomFlowRate = bc_heater.st_hot_in.mdot, fluidNomPressure = bc_heater.st_hot_in.p, exchSurface_G = geo_heater_cold.A_ex, exchSurface_F = geo_heater_hot.A_ex, extSurfaceTub = geo_heater_tube.A_ex, gasVol = geo_heater_cold.V, fluidVol = geo_heater_hot.V, metalVol = geo_heater_tube.V, SSInit = false, rhomcm = thermo_heater_tube.rho_mcm, lambda = thermo_heater_tube.lambda, Tstartbar_G = bc_heater.st_cold_in.T, Tstartbar_M = bc_heater.st_hot_in.T - 50, pstart_F = bc_heater.st_hot_in.p, FluidPhaseStart = ThermoPower.Choices.FluidPhase.FluidPhases.Liquid, redeclare replaceable model HeatTransfer_F = ThermoPower.Thermal.HeatTransferFV.IdealHeatTransfer, redeclare replaceable model HeatTransfer_G = ThermoPower.Thermal.HeatTransferFV.IdealHeatTransfer, redeclare model HeatExchangerTopology = ThermoPower.Thermal.HeatExchangerTopologies.CounterCurrentFlow, metalTube(WallRes = false, Tstartbar = bc_heater.st_hot_in.T - 50)) annotation(
    Placement(visible = true, transformation(origin = {-69, 33}, extent = {{-7, -7}, {7, 7}}, rotation = 180)));
  // set the fluid flow as fixed mdot for simplarity
  // ConstantHeatTransferCoefficientTwoGrids(gamma = thermo_HTR_hot.gamma_he),
  // ConstantHeatTransferCoefficient(gamma = thermo_HTR_cold.gamma_he),
  // use FlowJoin to mix flow
  ThermoPower.Gas.FlowJoin mixer(redeclare package Medium = medium_main) annotation(
    Placement(visible = true, transformation(origin = {-12, 50}, extent = {{-6, -6}, {6, 6}}, rotation = 180)));
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
  Steps.Components.HEG2G HTR(redeclare package FluidMedium = medium_main, redeclare package FlueGasMedium = medium_main, redeclare replaceable model HeatTransfer_F = ThermoPower.Thermal.HeatTransferFV.ConstantHeatTransferCoefficient(gamma = thermo_HTR_cold.gamma_he), redeclare replaceable model HeatTransfer_G = ThermoPower.Thermal.HeatTransferFV.ConstantHeatTransferCoefficientTwoGrids(gamma = thermo_HTR_hot.gamma_he), redeclare model HeatExchangerTopology = ThermoPower.Thermal.HeatExchangerTopologies.CounterCurrentFlow, N_F = geo_HTR_cold.N_seg, N_G = geo_HTR_hot.N_seg, Tstartbar_G = bc_HTR.st_hot_in.T, Tstartbar_F = bc_HTR.st_cold_in.T, exchSurface_F = geo_HTR_cold.A_ex, exchSurface_G = geo_HTR_hot.A_ex, extSurfaceTub = geo_HTR_tube.A_ex, fluidNomFlowRate = bc_HTR.st_cold_in.mdot, fluidNomPressure = bc_HTR.st_cold_in.p, fluidVol = geo_HTR_cold.V, gasNomFlowRate = bc_HTR.st_hot_in.mdot, gasNomPressure = bc_HTR.st_hot_in.p, gasVol = geo_HTR_hot.V, lambda = thermo_HTR_tube.lambda, metalVol = geo_HTR_tube.V, pstart_F = bc_HTR.st_cold_in.p, pstart_G = bc_HTR.st_hot_in.T, rhomcm = thermo_HTR_tube.rho_mcm, gasQuasiStatic = true, fluidQuasiStatic = true, metalTube(WallRes = false)) annotation(
    Placement(visible = true, transformation(origin = {-21, 33}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
  // redeclare replaceable model HeatTransfer_F = ThermoPower.Thermal.HeatTransferFV.IdealHeatTransfer,
  //redeclare replaceable model HeatTransfer_G = ThermoPower.Thermal.HeatTransferFV.IdealHeatTransfer,
  Steps.Components.HEG2G LTR(redeclare package FluidMedium = medium_main, redeclare package FlueGasMedium = medium_main, redeclare replaceable model HeatTransfer_F = ThermoPower.Thermal.HeatTransferFV.ConstantHeatTransferCoefficient(gamma = thermo_LTR_cold.gamma_he), redeclare replaceable model HeatTransfer_G = ThermoPower.Thermal.HeatTransferFV.ConstantHeatTransferCoefficientTwoGrids(gamma = thermo_LTR_hot.gamma_he), redeclare model HeatExchangerTopology = ThermoPower.Thermal.HeatExchangerTopologies.CounterCurrentFlow, N_F = geo_LTR_cold.N_seg, N_G = geo_LTR_hot.N_seg, Tstartbar_G = bc_LTR.st_hot_in.T, Tstartbar_F = bc_LTR.st_cold_in.T, exchSurface_F = geo_LTR_cold.A_ex, exchSurface_G = geo_LTR_hot.A_ex, extSurfaceTub = geo_LTR_tube.A_ex, fluidNomFlowRate = bc_LTR.st_cold_in.mdot, fluidNomPressure = bc_LTR.st_cold_in.p, fluidVol = geo_LTR_cold.V, gasNomFlowRate = bc_LTR.st_hot_in.mdot, gasNomPressure = bc_LTR.st_hot_in.p, gasVol = geo_LTR_hot.V, lambda = thermo_LTR_tube.lambda, metalVol = geo_LTR_tube.V, pstart_F = bc_LTR.st_cold_in.p, pstart_G = bc_LTR.st_hot_in.T, rhomcm = thermo_LTR_tube.rho_mcm, gasQuasiStatic = true, fluidQuasiStatic = true, metalTube(WallRes = false)) annotation(
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
  ThermoPower.Gas.SensT sens_turbine(redeclare package Medium = medium_main) annotation(
    Placement(visible = true, transformation(origin = {-44, -10}, extent = {{-6, -6}, {6, 6}}, rotation = 90)));
  ThermoPower.Water.SourceMassFlow source_cooler_cold(redeclare package Medium = medium_cooler, T = bc_cooler.st_cold_in.T, p0 = bc_cooler.st_cold_in.p, use_T = true, use_in_T = false, w0 = bc_cooler.st_cold_in.mdot) annotation(
    Placement(visible = true, transformation(origin = {51, -27}, extent = {{-5, -5}, {5, 5}}, rotation = 0)));
  ThermoPower.Water.SinkPressure sink_cooler_cold(redeclare package Medium = medium_cooler, p0 = bc_cooler.st_cold_out.p, T = bc_cooler.st_cold_out.T, use_T = true) annotation(
    Placement(visible = true, transformation(origin = {44, -64}, extent = {{-6, -6}, {6, 6}}, rotation = 180)));
  ThermoPower.PowerPlants.HRSG.Components.HE cooler(redeclare package FluidMedium = medium_cooler, redeclare package FlueGasMedium = medium_main, fluidFlow(fixedMassFlowSimplified = true, hstartin = bc_cooler.st_cold_in.h, hstartout = bc_cooler.st_cold_out.h), gasFlow(QuasiStatic = true), N_G = geo_cooler_hot.N_seg, N_F = geo_cooler_cold.N_seg, Nw_G = geo_cooler_tube.N_seg, gasNomFlowRate = bc_cooler.st_hot_in.mdot, gasNomPressure = bc_cooler.st_hot_in.p, fluidNomFlowRate = bc_cooler.st_cold_in.mdot, fluidNomPressure = bc_cooler.st_cold_in.p, exchSurface_G = geo_cooler_hot.A_ex, exchSurface_F = geo_cooler_cold.A_ex, extSurfaceTub = geo_cooler_tube.A_ex, gasVol = geo_cooler_hot.V, fluidVol = geo_cooler_cold.V, metalVol = geo_cooler_tube.V, SSInit = false, rhomcm = thermo_cooler_tube.rho_mcm "use thermo props of heater for simplicity", lambda = thermo_cooler_tube.lambda "use thermo props of heater for simplicity", Tstartbar_G = bc_cooler.st_hot_in.T, Tstartbar_M = bc_cooler.st_hot_in.T - 50, pstart_F = bc_cooler.st_cold_in.p, FluidPhaseStart = ThermoPower.Choices.FluidPhase.FluidPhases.Liquid, redeclare replaceable model HeatTransfer_F = ThermoPower.Thermal.HeatTransferFV.IdealHeatTransfer, redeclare replaceable model HeatTransfer_G = ThermoPower.Thermal.HeatTransferFV.IdealHeatTransfer, redeclare model HeatExchangerTopology = ThermoPower.Thermal.HeatExchangerTopologies.CounterCurrentFlow, metalTube(WallRes = false, Tstartbar = bc_cooler.st_hot_in.T - 50)) annotation(
    Placement(visible = true, transformation(origin = {57, -47}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
  //Components.HEG2G cooler(
  // set the fluid flow as fixed mdot for simplarity
  //ConstantHeatTransferCoefficientTwoGrids(gamma = thermo_hot.gamma_he),
  //ConstantHeatTransferCoefficient(gamma =  thermo_cold.gamma_he),
  ThermoPower.Gas.FlowSplit splitter(redeclare package Medium = medium_main) annotation(
    Placement(visible = true, transformation(origin = {35, 7}, extent = {{-5, -5}, {5, 5}}, rotation = -90)));
  /*
        ThermoPower.Gas.Compressor compressor(
          redeclare package Medium = medium_main,
          pstart_in=bc_cooler.st_hot_out.p,
          pstart_out=bc_LTR.st_cold_in.p,
          Tstart_in=bc_cooler.st_hot_out.T,
          Tstart_out=bc_LTR.st_cold_in.T,
          tablePhic=tablePhic_comp,
          tableEta=tableEta_comp,
          tablePR=tablePR_comp,
          Table=ThermoPower.Choices.TurboMachinery.TableTypes.matrix,
          Ndesign=523.3,
          Tdes_in=244.4) annotation (Placement(transformation(extent={{-20,-20},{
                  20,20}}, rotation=0)));
                  
        Modelica.Mechanics.Rotational.Sources.ConstantSpeed ConstantSpeed1(
            w_fixed=523.3, useSupport=false) annotation (Placement(transformation(
                extent={{-50,-10},{-30,10}}, rotation=0)));

        ThermoPower.Gas.Compressor compressor_bypass(
          redeclare package Medium = medium_main,
          pstart_in=bc_HTR.st_hot_out.p,
          pstart_out=st_bypass.p,
          Tstart_in=bc_HTR.st_hot_out.T,
          Tstart_out=st_bypass.T,
          tablePhic=tablePhic_comp,
          tableEta=tableEta_comp,
          tablePR=tablePR_comp,
          Table=ThermoPower.Choices.TurboMachinery.TableTypes.matrix,
          Ndesign=523.3,
          Tdes_in=244.4) annotation (Placement(transformation(extent={{-20,-20},{
                  20,20}}, rotation=0)));
                  
        Modelica.Mechanics.Rotational.Sources.ConstantSpeed ConstantSpeed2(
            w_fixed=523.3, useSupport=false) annotation (Placement(transformation(
                extent={{-50,-10},{-30,10}}, rotation=0)));
      */
  ThermoPower.Gas.SourceMassFlow source(redeclare package Medium = medium_main, w0 = bc_LTR.st_cold_in.mdot, p0 = bc_LTR.st_cold_in.p, T = bc_LTR.st_cold_in.T) annotation(
    Placement(visible = true, transformation(origin = {83, 31}, extent = {{-5, -5}, {5, 5}}, rotation = 90)));
  ThermoPower.Gas.SinkPressure sink(redeclare package Medium = medium_main, p0 = bc_cooler.st_hot_out.p, T = bc_cooler.st_hot_out.T) annotation(
    Placement(visible = true, transformation(origin = {82, 18}, extent = {{-4, -4}, {4, 4}}, rotation = 270)));
  /*
    ThermoPower.Gas.ThroughMassFlow source_through(redeclare package Medium = medium_main, w0 = bc_cooler.st_hot_out.mdot) annotation(
      Placement(visible = true, transformation(extent = {{50, 10}, {70, 30}}, rotation = 0)));
    ThermoPower.Gas.Utility.ClosedSystemInit sys_init(redeclare package Medium = medium_main, pstart = bc_LTR.st_cold_out.p);
    */
  //w_b(nominal = bc_heater.st_cold_out.mdot),
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
        */
  ThermoPower.Gas.ThroughMassFlow source_through_bypass(redeclare package Medium = medium_main, w0 = st_bypass.mdot) annotation(
    Placement(visible = true, transformation(origin = {61, 29}, extent = {{-5, -5}, {5, 5}}, rotation = 90)));
protected
  parameter Real tableEta_comp[6, 4] = [0, 95, 100, 105; 1, 82.5e-2, 81e-2, 80.5e-2; 2, 84e-2, 82.9e-2, 82e-2; 3, 83.2e-2, 82.2e-2, 81.5e-2; 4, 82.5e-2, 81.2e-2, 79e-2; 5, 79.5e-2, 78e-2, 76.5e-2];
  parameter Real tablePhic_comp[6, 4] = [0, 95, 100, 105; 1, 38.3e-3 / 400, 43e-3 / 400, 46.8e-3 / 400; 2, 39.3e-3 / 400, 43.8e-3 / 400, 47.9e-3 / 400; 3, 40.6e-3 / 400, 45.2e-3 / 400, 48.4e-3 / 400; 4, 41.6e-3 / 400, 46.1e-3 / 400, 48.9e-3 / 400; 5, 42.3e-3 / 400, 46.6e-3 / 400, 49.3e-3 / 400];
  parameter Real tablePR_comp[6, 4] = [0, 95, 100, 105; 1, 22.6 / 12, 27 / 12, 32 / 12; 2, 22 / 12, 26.6 / 12, 30.8 / 12; 3, 20.8 / 12, 25.5 / 12, 29 / 12; 4, 19 / 12, 24.3 / 12, 27.1 / 12; 5, 17 / 12, 21.5 / 12, 24.2 / 12];
  parameter Real tablePhic_turbine[5, 4] = [1, 37, 80, 100; 1.5, 7.10E-05, 7.10E-05, 7.10E-05; 2, 8.40E-05, 8.40E-05, 8.40E-05; 2.5, 8.70E-05, 8.70E-05, 8.70E-05; 3, 1.04E-04, 1.04E-04, 1.04E-04];
  parameter Real tableEta_turbine[5, 4] = [1, 37, 80, 100; 1.5, 0.57, 0.89, 0.81; 2, 0.46, 0.82, 0.88; 2.5, 0.41, 0.76, 0.85; 3, 0.38, 0.72, 0.82];
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
// Close loop without compressor
  connect(source.flange, LTR.waterIn) annotation(
    Line(points = {{83, 36}, {83, 61}, {21, 61}, {21, 40}}, color = {0, 0, 255}, thickness = 0.5));
  connect(LTR.waterOut, mixer.inlet1) annotation(
    Line(points = {{21, 26}, {21.25, 26}, {21.25, 20}, {-0.5, 20}, {-0.5, 48}, {-8, 48}}, color = {0, 0, 255}, thickness = 0.5));
  connect(source_through_bypass.outlet, mixer.inlet2) annotation(
    Line(points = {{61, 34}, {61, 52}, {-8, 52}}, color = {0, 0, 255}, thickness = 0.5));
  connect(heater.gasOut, turbine.inlet) annotation(
    Line(points = {{-76, 33}, {-76, 33.75}, {-92, 33.75}, {-92, -41.5}, {-78, -41.5}, {-78, -42}}, color = {0, 0, 255}, thickness = 0.5));
  connect(turbine.shaft_b, speed1.flange) annotation(
    Line(points = {{-62, -51}, {-51, -51}, {-51, -50}, {-52, -50}}));
  connect(speed1.w_ref, const1.y) annotation(
    Line(points = {{-39, -50}, {-33, -50}}, color = {0, 0, 127}));
  connect(sens_turbine.outlet, HTR.gasIn) annotation(
    Line(points = {{-42, -6}, {-42, 34}, {-28, 34}}, color = {159, 159, 223}));
  connect(HTR.gasOut, LTR.gasIn) annotation(
    Line(points = {{-14, 33}, {14, 33}}, color = {159, 159, 223}));
  connect(splitter.outlet1, cooler.gasIn) annotation(
    Line);
  connect(cooler.gasOut, sink.flange) annotation(
    Line(points = {{64, -47}, {82, -47}, {82, 22}}, color = {159, 159, 223}));
  connect(sys_init.flange, mixer.inlet1) annotation(
    Line);
  connect(splitter.outlet2, source_through_bypass.inlet) annotation(
    Line);
// hot stream for heater
  connect(source_heater_hot.flange, heater.waterIn) annotation(
    Line(points = {{-70, 14}, {-68.5, 14}, {-68.5, 26}, {-69, 26}}, color = {0, 0, 255}));
  connect(heater.waterOut, sink_heater_hot.flange) annotation(
    Line);
// cold stream for cooler
  connect(cooler.waterOut, sink_cooler_cold.flange) annotation(
    Line(points = {{50, -64}, {50, -60.5}, {58, -60.5}, {58, -54}}, color = {0, 0, 255}));
/*
  // open loop without compressor
  connect(source.flange, LTR.waterIn);    
  connect(LTR.waterOut, mixer.inlet1);
  connect(source_through_bypass.outlet, mixer.inlet2);

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
  connect(LTR.gasOut, splitter.inlet);
  
  connect(splitter.outlet1, cooler.gasIn);
  connect(cooler.gasOut, sink.flange);  
  connect(splitter.outlet2, source_through_bypass.inlet);
   
  // hot stream for heater
  connect(source_heater_hot.flange, heater.waterIn);
  connect(heater.waterOut, sink_heater_hot.flange);

  // cold stream for cooler
  connect(source_cooler_cold.flange, cooler.waterIn);
  connect(cooler.waterOut, sink_cooler_cold.flange);
*/
  connect(HTR.gasOut, LTR.gasIn) annotation(
    Line(points = {{-14, 33}, {14, 33}}, color = {159, 159, 223}));
  connect(source_heater_hot.flange, heater.waterIn) annotation(
    Line(points = {{-70, 14}, {-68.5, 14}, {-68.5, 26}, {-69, 26}}, color = {0, 0, 255}));
  connect(splitter.outlet2, cooler.gasIn) annotation(
    Line(points = {{33, 4}, {33, -47.5}, {50, -47.5}, {50, -47}}, color = {159, 159, 223}));
  connect(splitter.outlet1, source_through_bypass.inlet) annotation(
    Line(points = {{37, 4}, {61, 4}, {61, 24}}, color = {159, 159, 223}));
  connect(source_through_bypass.outlet, mixer.inlet2) annotation(
    Line(points = {{61, 34}, {61, 52}, {-8, 52}}, color = {0, 0, 255}, thickness = 0.5));
  connect(LTR.waterOut, mixer.inlet1) annotation(
    Line(points = {{21, 26}, {21.25, 26}, {21.25, 20}, {-0.5, 20}, {-0.5, 48}, {-8, 48}}, color = {0, 0, 255}, thickness = 0.5));
  connect(cooler.gasOut, sink.flange) annotation(
    Line(points = {{64, -47}, {82, -47}, {82, 22}}, color = {159, 159, 223}));
  connect(source.flange, LTR.waterIn) annotation(
    Line(points = {{83, 36}, {83, 61}, {21, 61}, {21, 40}}, color = {0, 0, 255}, thickness = 0.5));
  connect(HTR.waterOut, heater.gasIn) annotation(
    Line(points = {{-21, 26}, {-27.5, 26}, {-27.5, 20}, {-50, 20}, {-50, 34}, {-62, 34}}, color = {0, 0, 255}, thickness = 0.5));
  connect(mixer.outlet, HTR.waterIn) annotation(
    Line(points = {{-16, 50}, {-21, 50}, {-21, 40}}, color = {0, 0, 255}, thickness = 0.5));
  connect(turbine.outlet, sens_turbine.inlet) annotation(
    Line(points = {{-60, -42}, {-42, -42}, {-42, -14}, {-42, -14}}, color = {159, 159, 223}));
  connect(sens_turbine.outlet, HTR.gasIn) annotation(
    Line(points = {{-42, -6}, {-42, 34}, {-28, 34}}, color = {159, 159, 223}));
  connect(splitter.inlet, LTR.gasOut) annotation(
    Line(points = {{36, 10}, {36, 34}, {28, 34}}, color = {159, 159, 223}));
  connect(sink_cooler_cold.flange, cooler.waterOut) annotation(
    Line(points = {{50, -64}, {50, -60.5}, {58, -60.5}, {58, -54}}, color = {0, 0, 255}));
  connect(sink_heater_hot.flange, heater.waterOut) annotation(
    Line(points = {{-74, 50}, {-68, 50}, {-68, 40}, {-68, 40}}, color = {0, 0, 255}));
  connect(source_cooler_cold.flange, cooler.waterIn) annotation(
    Line(points = {{56, -26}, {60, -26}, {60, -40}, {58, -40}}, color = {0, 0, 255}));
  annotation(
    Diagram(coordinateSystem(initialScale = 0.5)),
    experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-3, Interval = 1),
    __OpenModelica_commandLineOptions = "--matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection -d=initialization,NLSanalyticJacobian,bltdump",
    __OpenModelica_simulationFlags(lv = "LOG_DEBUG,LOG_NLS,LOG_NLS_V,LOG_STATS,LOG_INIT,LOG_STDOUT, -w", outputFormat = "mat", s = "dassl", nls = "homotopy"),
    Icon(coordinateSystem(initialScale = 0.5)));
end TP_RCBCycleGUI;
