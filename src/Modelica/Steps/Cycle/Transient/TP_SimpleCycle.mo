within Steps.Cycle.Transient;

model TP_SimpleCycle
  "simple cycle built referring Bone's Study [bone2021] (with an extra cooler) comp-hx-turbine-cooler"

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
  // package medium_main = Modelica.Media.IdealGases.SingleGases.CO2;
  package medium_main =  ExternalMedia.Examples.CO2CoolProp;
  // package medium_hot  = Steps.Media.SCO2;
  // package medium_main = Steps.Media.SCO2;
  
  // geometry parameters
  constant  Real pi           = Modelica.Constants.pi;

  parameter Integer N_ch      = integer(2400) "channel number";
  parameter Integer N_seg     = 10 "number of segments in one tube";
  parameter SI.Length D_ch    = 1.72e-3  "channel diameter, semi circular tube";
  parameter SI.Length r_ch    = D_ch / 2 "channel radiaus";
  parameter SI.Length L_fp    = 250 * 1e-3 "equivalent valid channel flow path length";
  parameter SI.Length L_pitch = 12e-3 "pitch length";
  parameter Real a_phi        = 36 "pitch angle degree";
  parameter SI.Length H_ch    = 4.17e-3 "Height of the solid domain, containing one cold tube and one hot tube";
  parameter SI.Length W_ch    = 2.3e-3 "Width of the solid domain";
  parameter SI.Length T_wall  = 0.51e-3 "Wall thinckness";
  parameter SI.Length L_wall  = 420e-3 "Length of wall, not necessarily equals to length of flow path";
  parameter SI.Area A         = pi * r_ch ^2 / 2 "Area of cross section of semi circular tube";

  // boundary conditon  
  // zigzag higher T
  // parameter SI.Velocity u_hot_in = 7.564 "hot inlet velocity m/s";
  // parameter SI.Velocity u_cold_in = 1.876 "cold inlet velocity m/s";
  parameter SI.Pressure p_hot_in                       = 4e6;      //from_bar(10) "hot inlet pressure - Irrelevant for Incompressible Thermia Oil" ;
  parameter SI.Pressure p_source                       = 8.65e6 ;
  parameter SI.Pressure p_comp_out                     = 12e6;     // 12.567478e6 "cold inlet pressure";  
  parameter SI.Pressure p_sink                         = 9.2e6;    // Not used by Config. Record only;
  parameter Modelica.SIunits.Temperature T_source      = 322;
  parameter Modelica.SIunits.Temperature T_comp_out    = 348;
  parameter Modelica.SIunits.Temperature T_heater_hin  = 595;
  parameter Modelica.SIunits.Temperature T_heater_hout = 478;
  parameter Modelica.SIunits.Temperature T_heater_cin  = 350;
  parameter Modelica.SIunits.Temperature T_heater_cout = 571;
  parameter Modelica.SIunits.Temperature T_turb_in     = 570;
  parameter Modelica.SIunits.Temperature T_sink        = 545;
  
  // pressure drop correction coefficient 
  parameter Real Cf_C1 = 1, Cf_C2 = 1, Cf_C3 = 1;
  //parameter Real Cf_a_cold = 1, Cf_b_cold = 1, Cf_c_cold = 1;  
  
  // meshram's cp and rho for alloy Inconel 617
  // parameter Modelica.SIunits.Density rho_wall = 8360 "density of wall, kg/m3";
  // parameter Modelica.SIunits.SpecificHeatCapacity cp_wall = 417 "cp of wall, J/kg-K";  
   
  // Stainless 316, 316L, 317, 317L
  parameter Modelica.SIunits.Density rho_wall             = 8030 "density of wall, kg/m3";
  parameter Modelica.SIunits.SpecificHeatCapacity cp_wall = 485 "cp of wall, J/kg-K";
  // thermal conductivity (T in K) https://www.theworldmaterial.com/aisi-316-ss316-stainless-steel-properties-composition/
  // parameter Real table_k_metalwall[:,:] = [20, 12.1; 100, 16.3; 500, 21.5];
  parameter Real table_k_metalwall[:,:] = [293.15, 12.1; 373.15, 16.3; 773.15, 21.5];

  parameter SI.MassFlowRate mdot_heater = 10;    //1.218944;
  parameter SI.MassFlowRate mdot_main   = 10.5;  // 0.0854299999999999;
  
  parameter Real N_s_compressor = 2100 "rotational speed of compressor";

  //Components
  // for transient simulation, set initOpt = steadyState and HE's SSInit = true
  inner ThermoPower.System system(allowFlowReversal = false, initOpt = ThermoPower.Choices.Init.Options.steadyState) annotation(
    Placement(transformation(extent = {{80, 80}, {100, 100}})));  
  
  // local reference of the config objects
  parameter Model.SimpleCycleConfig cfg(
  redeclare package medium_main = medium_main,
  redeclare package medium_heater = medium_hot,
    mdot_heater = mdot_heater
  );
  
  parameter Model.HeatExchangerConfig cfg_heater = cfg.cfg_heater;
  parameter Model.TurbomachineryConfig cfg_comp  = cfg.cfg_comp;
  parameter Model.TurbomachineryConfig cfg_turb  = cfg.cfg_turb;

  // key point thermodynamic state, may be altered for different cases
  parameter Model.ThermoState st_source = cfg.st_source;  
  parameter Model.ThermoState st_sink = cfg.st_comp_in;   
  
  ThermoPower.Gas.SourceMassFlow source_cold(
    redeclare package Medium = medium_main,
    use_in_T = false,
    T        = st_source.T,
    p0       = st_source.p,
    w0       = st_source.mdot,
    // T = 322,
    // p0 = 8.65e6,
    // w0 = 10.5,
    // gas(p(nominal = bc_HE.st_cold_in.p), 
    gas(
      p(nominal = st_source.p), 
      // p(nominal = 8.65e6), 
      // T(nominal=bc_HE.st_cold_in.T))) 
      T(nominal = st_source.T))
      // T(start = 322, nominal=322))
    ) 
  annotation(
    Placement(transformation(origin = {0, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 270)));
  
  ThermoPower.Gas.SinkPressure sink_cold(
    redeclare package Medium = medium_main,
    // p0 = 12e6, 
    // T = 348,
    p0 = st_sink.p,
    T  = st_sink.T,
    gas(
      p(nominal = st_sink.p), 
      T(nominal = st_sink.T))
      // T(start = 348, nominal=348))    
  );
/*
  ThermoPower.Gas.SinkPressure sink_hot(
    redeclare package Medium = medium_hot,
    T  = cfg.st_heater_hout.T,
    p0 = cfg.st_heater_hout.p)
  annotation(
    Placement(transformation(extent = {{60, -10}, {80, 10}}, rotation = 0)));
  
  ThermoPower.Gas.SourceMassFlow source_hot(
    redeclare package Medium = medium_hot,
    T         = cfg.st_heater_hin.T,
    p0        = cfg.st_heater_hin.p,
    w0        = cfg.st_heater_hin.mdot,
    use_in_T  = false,
    use_in_w0 = true,
    gas(
      p(nominal = cfg.st_heater_hin.p), 
      T(nominal=cfg.st_heater_hin.T))) 
  annotation(
    Placement(transformation(extent = {{-70, -10}, {-50, 10}}, rotation = 0))); 
*/

  ThermoPower.Gas.SensT T_waterIn(redeclare package Medium = medium_main) annotation(
    Placement(transformation(origin = {4, -50}, extent = {{-10, -10}, {10, 10}}, rotation = 270)));
  ThermoPower.Gas.SensT T_waterOut(redeclare package Medium = medium_main) annotation(
    Placement(transformation(origin = {4, -50}, extent = {{-10, -10}, {10, 10}}, rotation = 270)));

/*
  ThermoPower.Gas.SensT T_gasIn(redeclare package Medium = medium_hot);
  ThermoPower.Gas.SensT T_gasOut(redeclare package Medium = medium_hot);
*/
/*
  Steps.TPComponents.PCHE HE(
    redeclare package FluidMedium   = medium_main,
    redeclare package FlueGasMedium = medium_hot,     
    // use Marchionni PCHE HeatTransfer
    // slow but can have a result - set a_phi = 0 to use Gnielinski's correlation 
    redeclare replaceable model HeatTransfer_F = Steps.TPComponents.MarchionniPCHEHeatTransferFV(),
    redeclare replaceable model HeatTransfer_G = Steps.TPComponents.MarchionniPCHEHeatTransferFV(),
    gasFlow(heatTransfer(pitch = cfg.pitch, phi = cfg.phi, Cf_C1 = Cf_C1, Cf_C2 = Cf_C2, Cf_C3 = Cf_C3)),
    fluidFlow(heatTransfer(pitch = cfg.pitch, phi = cfg.phi, Cf_C1 = Cf_C1, Cf_C2 = Cf_C2, Cf_C3 = Cf_C3)),    
    
    // fast and works fine for now. Error occurs when mass flow rate is zero, i.e. one flow is shut down. 
    // redeclare replaceable model HeatTransfer_F = ThermoPower.Thermal.HeatTransferFV.IdealHeatTransfer,
    // redeclare replaceable model HeatTransfer_G = ThermoPower.Thermal.HeatTransferFV.IdealHeatTransfer,   
    
    cfg = cfg.cfg_heater,
    N_G = cfg.gp_heater_flow.N_seg, // N_G and N_F has to be assigned explicitly since they set dimension of Gas and Fluid cells as gasflow.gas[N_G] and fluidflow.gas[N_F]
    N_F = cfg.gp_heater_flow.N_seg,
    SSInit            = true,
    gasQuasiStatic    = false,
    fluidQuasiStatic  = false,
    metalWall(L = L_wall, w_ch = W_ch, h_ch = H_ch, dx = T_wall)
    // metalQuasiStatic = true
    // override the values of Am and L of metaltubeFV
    // to make them agree with semi-circular tube of PCHE
    // ('final' modifier of Am in metalTubeFv was removed as well)
    //metalTube(WallRes=false, L = 1, rhomcm=200, Am = HE.metalVol / 1) 
  )
  annotation(
    Placement(transformation(extent = {{-20, -20}, {20, 20}}, rotation = 0)));  
*/
  ThermoPower.Gas.Compressor compressor(
    redeclare package Medium = medium_main,
    // pstart_in  = 8.65e6,
    // pstart_out = 12e6,
    // Tstart_in  = 322,
    // Tstart_out = 348,
    pstart_in  = cfg_comp.st_in.p,
    Tstart_in  = cfg_comp.st_in.T,
    pstart_out = cfg_comp.st_out.p,
    Tstart_out = cfg_comp.st_out.T,
    tablePhic  = tablePhic_comp_mc,
    tableEta   = tableEta_comp_mc,
    tablePR    = tablePR_comp_mc,
    Table      = ThermoPower.Choices.TurboMachinery.TableTypes.matrix,
    Ndesign    = cfg_comp.N,
    explicitIsentropicEnthalpy = false, 
    Tdes_in    = compressor.Tstart_in)
    annotation(
    Placement(visible = true, transformation(origin = {103, -11}, extent = {{-11, -11}, {11, 11}}, rotation = 0)));

  Modelica.Mechanics.Rotational.Sources.Torque torque_comp(useSupport=false);
  Modelica.Mechanics.Rotational.Components.Inertia Inertia1(J = 0.7);

/*  
  Modelica.Mechanics.Rotational.Sources.ConstantSpeed const_speed_comp(
      w_fixed=N_s_compressor, useSupport=false) annotation(
    Placement(visible = true, transformation(origin = {81, -7}, extent = {{-5, -5}, {5, 5}}, rotation = 0)));  
*/  
/*  
  ThermoPower.Gas.Turbine turbine(
  redeclare package Medium = medium_main, 
  fileName = Modelica.Utilities.Files.loadResource("modelica://Steps/Resources/Data/turbine_map.txt"),   
  tablePhic = fill(0.0, 14, 12), //tablePhic, un
  tableEta = fill(0.0, 14, 12), //tableEta, 
  pstart_in = cfg_turb.st_in.p, 
  pstart_out = cfg_turb.st_out.p, 
  Tstart_in = cfg_turb.st_in.T, 
  Tstart_out = cfg_turb.st_out.T, 
  Ndesign = cfg_turb.N, 
  Tdes_in = cfg_turb.st_in.T,  
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
      w_fixed=cfg_turb.N, useSupport=false) annotation(
    Placement(visible = true, transformation(origin = {-46, -50}, extent = {{6, -6}, {-6, 6}}, rotation = 0)));
   
*/ 
/*  
  Modelica.Blocks.Sources.TimeTable tt_mdot_hot_in(
    startTime = 0, 
    table = [
      0, 10; 20.68, 9.951; 30.59, 11.422; 46.46, 12.843;
      64.87, 13.284; 90.08, 13.382; 95.47, 9.167; 96.6, 8.627;
      104.53, 7.745; 116.71, 6.814; 131.73, 6.814; 145.89, 6.961;
      160.34, 7.01; 169.69, 8.873; 180.45, 9.853; 200.28, 10.049;
      230.31, 9.951; 250.71, 13.922; 262.04, 14.657; 272.8, 15.049;
      288.67, 15.392; 300, 15.49
    ]) annotation(
    Placement(visible = true, transformation(origin = {76, 78}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
*/

  // hot inlet temperature in K
  Modelica.Blocks.Sources.TimeTable tt_T_motor(    
    table = [
      0, 84.615; 2.55, 85.714; 20.09, 84.615; 23.21, 124.725;
      26.89, 98.901; 29.15, 100; 90, 94.505; 93.96, 35.165;
      98.49, 72.527; 101.32, 69.231; 159.91, 71.978; 163.3, 109.341;
      167.26, 85.714; 171.23, 86.813; 229.81, 84.615; 233.77, 135.714;
      237.45, 107.692; 241.42, 108.791; 243.96, 100; 246.79, 106.044;
      250.47, 102.747; 299.72, 100.549
    ]) annotation(
    Placement(visible = true, transformation(origin = {76, 44}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
   
  
  /*
  // variable for validation
  Modelica.SIunits.Power Q_out    = (HE.gasIn.h_outflow - HE.gasOut.h_outflow) * HE.gasIn.m_flow;
  Modelica.SIunits.Power Q_in     = (HE.waterOut.h_outflow - HE.waterIn.h_outflow) * HE.waterIn.m_flow;
  Boolean                isQMatch = abs(Q_out - Q_in) < 1e-3;
  */
  
protected
  // performance map for main compressor
  parameter Real tableEta_comp_mc[:, :] = [
  0,85.0,90.0,95.0,100.0,105.0,110.0;
  1,0.771,0.777,0.785,0.796,0.804,0.808;
  2,0.8,0.804,0.809,0.812,0.815,0.818;
  3,0.813,0.816,0.819,0.822,0.824,0.827;
  4,0.818,0.823,0.826,0.83,0.83,0.83;
  5,0.818,0.829,0.83,0.83,0.83,0.83;
  6,0.818,0.827,0.83,0.83,0.83,0.83;
  7,0.815,0.82,0.824,0.829,0.83,0.83;
  8,0.8,0.808,0.814,0.817,0.824,0.829;
  9,0.766,0.78,0.787,0.804,0.809,0.815;
  10,0.725,0.745,0.762,0.771,0.779,0.796];

  parameter Real tablePhic_comp_mc[:, :] = [
  0,85.0,90.0,95.0,100.0,105.0,110.0;
  1,1.54e-05,1.66e-05,1.78e-05,1.9e-05,2.02e-05,2.14e-05;
  2,1.66e-05,1.78e-05,1.9e-05,2.02e-05,2.14e-05,2.26e-05;
  3,1.78e-05,1.9e-05,2.02e-05,2.14e-05,2.26e-05,2.38e-05;
  4,1.9e-05,2.02e-05,2.14e-05,2.26e-05,2.38e-05,2.5e-05;
  5,2.02e-05,2.14e-05,2.26e-05,2.38e-05,2.5e-05,2.61e-05;
  6,2.14e-05,2.26e-05,2.38e-05,2.5e-05,2.61e-05,2.73e-05;
  7,2.26e-05,2.38e-05,2.5e-05,2.61e-05,2.73e-05,2.85e-05;
  8,2.38e-05,2.5e-05,2.61e-05,2.73e-05,2.85e-05,2.97e-05;
  9,2.5e-05,2.61e-05,2.73e-05,2.85e-05,2.97e-05,3.09e-05;
  10,2.61e-05,2.73e-05,2.85e-05,2.97e-05,3.09e-05,3.21e-05];

  parameter Real tablePR_comp_mc[:, :] = [
  0,85.0,90.0,95.0,100.0,105.0,110.0;
  1,1.63,1.72,1.83,1.95,2.08,2.21;
  2,1.65,1.74,1.85,1.96,2.09,2.22;
  3,1.65,1.75,1.85,1.96,2.09,2.22;
  4,1.65,1.75,1.85,1.96,2.09,2.21;
  5,1.64,1.74,1.84,1.95,2.07,2.2;
  6,1.62,1.72,1.83,1.93,2.05,2.17;
  7,1.6,1.7,1.79,1.9,2.02,2.14;
  8,1.56,1.66,1.75,1.85,1.97,2.1;
  9,1.5,1.6,1.68,1.79,1.92,2.04;
  10,1.45,1.53,1.62,1.72,1.83,1.97];  

initial equation
  Inertia1.w = N_s_compressor;
  
equation
/*
  // HE alone
  connect(source_cold.flange, T_waterIn.inlet);
  connect(T_waterIn.outlet, HE.waterIn) annotation(
    Line(points = {{-1.83697e-015, 50}, {-1.83697e-015, 20}, {0, 20}}, color = {0, 0, 255}, thickness = 0.5, smooth = Smooth.None));
  connect(HE.waterOut, T_waterOut.inlet);
  connect(T_waterOut.outlet, sink_cold.flange) annotation(
    Line(points = {{1.83697e-015, -70}, {1.83697e-015, -56}, {-8.88178e-016, -56}}, thickness = 0.5, color = {0, 0, 255}));
  
  connect(tt_mdot_hot_in.y, source_hot.in_w0);
  connect(source_hot.flange, T_gasIn.inlet);
  connect(T_gasIn.outlet, HE.gasIn) annotation(
        Line(points = {{-50, 0}, {-20, 0}}, color = {159, 159, 223}, thickness = 0.5, smooth = Smooth.None)); 
  connect(HE.gasOut, T_gasOut.inlet);
  connect(T_gasOut.outlet, sink_hot.flange) annotation(
    Line(points = {{46, 0}, {46, 0}, {60, 0}}, color = {159, 159, 223}, thickness = 0.5));    
*/   


  // compressor alone
  connect(source_cold.flange, T_waterIn.inlet);
  connect(T_waterIn.outlet, compressor.inlet);    
  connect(compressor.outlet, T_waterOut.inlet);
  connect(T_waterOut.outlet, sink_cold.flange);

  connect(tt_T_motor.y, torque_comp.tau);
  connect(torque_comp.flange, Inertia1.flange_b);
  connect(Inertia1.flange_a, compressor.shaft_b) annotation(
    Line(points = {{86, -6}, {96, -6}, {96, -10}, {96, -10}}));

  //connect(const_speed_comp.flange, compressor.shaft_a);
  
/*
  // heater + compressor
  connect(source_cold.flange, T_waterIn.inlet);
  connect(T_waterIn.outlet, compressor.inlet);
  
  connect(compressor.outlet, HE.waterIn) annotation(
    Line(points = {{-1.83697e-015, 50}, {-1.83697e-015, 20}, {0, 20}}, color = {0, 0, 255}, thickness = 0.5, smooth = Smooth.None));
  connect(HE.waterOut, T_waterOut.inlet) annotation(
    Line(points = {{8.88178e-016, -44}, {8.88178e-016, -20}, {0, -20}}, thickness = 0.5, color = {0, 0, 255}));      
  connect(T_waterOut.outlet, sink_cold.flange) annotation(
    Line(points = {{1.83697e-015, -70}, {1.83697e-015, -56}, {-8.88178e-016, -56}}, thickness = 0.5, color = {0, 0, 255}));
  
  connect(const_speed_comp.flange, compressor.shaft_a) annotation(
    Line(points = {{86, -6}, {96, -6}, {96, -10}, {96, -10}}));
  
  connect(source_hot.flange, T_gasIn.inlet);
  connect(T_gasIn.outlet, HE.gasIn) annotation(
        Line(points = {{-50, 0}, {-20, 0}}, color = {159, 159, 223}, thickness = 0.5, smooth = Smooth.None)); 
  connect(HE.gasOut, T_gasOut.inlet) annotation(
    Line(points = {{34, 0}, {34, 0}, {20, 0}}, color = {159, 159, 223}, thickness = 0.5));
  connect(T_gasOut.outlet, sink_hot.flange) annotation(
    Line(points = {{46, 0}, {46, 0}, {60, 0}}, color = {159, 159, 223}, thickness = 0.5));     
*/
/*
  // heater + compressor + turbine
  connect(source_cold.flange, T_waterIn.inlet);
  connect(T_waterIn.outlet, compressor.inlet);
  
  connect(compressor.outlet, HE.waterIn) annotation(
    Line(points = {{-1.83697e-015, 50}, {-1.83697e-015, 20}, {0, 20}}, color = {0, 0, 255}, thickness = 0.5, smooth = Smooth.None));
  connect(HE.waterOut, turbine.inlet);
  
  connect(turbine.outlet, T_waterOut.inlet) annotation(
    Line(points = {{8.88178e-016, -44}, {8.88178e-016, -20}, {0, -20}}, thickness = 0.5, color = {0, 0, 255}));      
  connect(T_waterOut.outlet, sink_cold.flange) annotation(
    Line(points = {{1.83697e-015, -70}, {1.83697e-015, -56}, {-8.88178e-016, -56}}, thickness = 0.5, color = {0, 0, 255}));
  

  connect(tt_T_motor_hot_in.y, torque_comp.tau);
  connect(torque_comp.flange, compressor.shaft_a) annotation(
    Line(points = {{86, -6}, {96, -6}, {96, -10}, {96, -10}}));
 
  connect(const_speed_turbine.flange, turbine.shaft_a);
  
  connect(tt_mdot_hot_in.y, source_hot.in_w0);
  connect(source_hot.flange, T_gasIn.inlet);
  connect(T_gasIn.outlet, HE.gasIn) annotation(
        Line(points = {{-50, 0}, {-20, 0}}, color = {159, 159, 223}, thickness = 0.5, smooth = Smooth.None)); 
  connect(HE.gasOut, T_gasOut.inlet) annotation(
    Line(points = {{34, 0}, {34, 0}, {20, 0}}, color = {159, 159, 223}, thickness = 0.5));
  connect(T_gasOut.outlet, sink_hot.flange) annotation(
    Line(points = {{46, 0}, {46, 0}, {60, 0}}, color = {159, 159, 223}, thickness = 0.5));     
*/   
annotation(
    Diagram(graphics),
    // for steady-state simulation - value check
    experiment(StartTime = 0, StopTime = 10, Tolerance = 1e-3, Interval = 2),
    // for complete transient simulation
    // experiment(StartTime = 0, StopTime = 600, Tolerance = 1e-3, Interval = 10),
    __OpenModelica_commandLineOptions = "--matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection -d=initialization,NLSanalyticJacobian,aliasConflicts",
    // remove the option flag --matchingAlgorithm=PFPlusExt, which may lead to 'Internal error - IndexReduction.dynamicStateSelectionWork failed!' during Translation
    // __OpenModelica_commandLineOptions = "--matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection -d=initialization,NLSanalyticJacobian,aliasConflicts",    
    __OpenModelica_simulationFlags(lv = "LOG_DEBUG,LOG_NLS,LOG_NLS_V,LOG_STATS,LOG_INIT,LOG_STDOUT, -w", outputFormat = "mat", s = "dassl", nls = "homotopy")
    );
end TP_SimpleCycle;
