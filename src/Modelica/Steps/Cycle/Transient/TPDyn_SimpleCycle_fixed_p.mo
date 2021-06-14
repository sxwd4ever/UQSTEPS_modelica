within Steps.Cycle.Transient;

model TPDyn_SimpleCycle_fixed_p
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
  // package medium_main =  ExternalMedia.Examples.CO2CoolProp;
  // package medium_hot  = Steps.Media.SCO2;
  package medium_main = Steps.Media.SCO2(
    inputChoice = ExternalMedia.Common.InputChoice.pT
  );
  
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
  parameter Modelica.SIunits.Temperature T_comp_out    = 428;
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
    mdot_heater = mdot_heater,
    p_comp_out = p_comp_out
  );
  
  parameter Model.HeatExchangerConfig cfg_heater = cfg.cfg_heater;
  parameter Model.TurbomachineryConfig cfg_comp  = cfg.cfg_comp;
  parameter Model.TurbomachineryConfig cfg_turb  = cfg.cfg_turb;

  // key point thermodynamic state, may be altered for different cases
  parameter Model.ThermoState st_source = cfg_comp.st_in;  
  parameter Model.ThermoState st_sink = cfg_turb.st_out;   
  
  /*   
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
  */ 

  ThermoPower.Gas.SourceMassFlow source_hot(
    redeclare package Medium = medium_hot,
    T         = cfg_heater.cfg_hot.st_in.T,
    p0        = cfg_heater.cfg_hot.st_in.p,
    w0        = cfg_heater.cfg_hot.st_in.mdot,
    use_in_T  = false,
    use_in_w0 = true,
    gas(
      p(nominal = cfg_heater.cfg_hot.st_in.p), 
      T(nominal = cfg_heater.cfg_hot.st_in.T))) 
  annotation(
    Placement(transformation(extent = {{-70, -10}, {-50, 10}}, rotation = 0)));       
  
  ThermoPower.Gas.SinkPressure sink_hot(
    redeclare package Medium = medium_hot,
    T  = cfg_heater.cfg_hot.st_out.T,
    p0 = cfg_heater.cfg_hot.st_out.p,
    gas(
      p(nominal = cfg_heater.cfg_hot.st_out.p), 
      T(nominal = cfg_heater.cfg_hot.st_out.T))
    )
  annotation(
    Placement(transformation(extent = {{60, -10}, {80, 10}}, rotation = 0)));
 
  ThermoPower.Gas.SourcePressure source_cold(
    redeclare package Medium = medium_main,
    use_in_T = false,
    T        = st_source.T,
    p0       = st_source.p,
    // w0       = st_source.mdot,
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
  
  ThermoPower.Gas.SensT T_waterIn(redeclare package Medium = medium_main) annotation(
    Placement(transformation(origin = {4, -50}, extent = {{-10, -10}, {10, 10}}, rotation = 270)));
  ThermoPower.Gas.SensT T_waterOut(redeclare package Medium = medium_main) annotation(
    Placement(transformation(origin = {4, -50}, extent = {{-10, -10}, {10, 10}}, rotation = 270)));

  ThermoPower.Gas.SensT T_gasIn(redeclare package Medium = medium_hot);
  ThermoPower.Gas.SensT T_gasOut(redeclare package Medium = medium_hot);

  Steps.TPComponents.PCHE HE(
    redeclare package FluidMedium   = medium_main,
    redeclare package FlueGasMedium = medium_hot,     
    // use Marchionni PCHE HeatTransfer
    // slow but can have a result - set a_phi = 0 to use Gnielinski's correlation 
    redeclare replaceable model HeatTransfer_F = Steps.TPComponents.MarchionniPCHEHeatTransferFV(),
    redeclare replaceable model HeatTransfer_G = Steps.TPComponents.MarchionniPCHEHeatTransferFV(),
    gasFlow(
      heatTransfer(
        pitch = cfg.pitch,
        phi   = cfg.phi,
        Cf_C1 = Cf_C1,
        Cf_C2 = Cf_C2,
        Cf_C3 = Cf_C3)),
    fluidFlow(
      heatTransfer(
        pitch = cfg.pitch,
        phi   = cfg.phi,
        Cf_C1 = Cf_C1,
        Cf_C2 = Cf_C2,
        Cf_C3 = Cf_C3)),
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

  ThermoPower.Gas.Compressor compressor(
    redeclare package Medium = medium_main,
    pstart_in                  = cfg_comp.st_in.p,
    Tstart_in                  = cfg_comp.st_in.T,
    pstart_out                 = cfg_comp.st_out.p,
    Tstart_out                 = cfg_comp.st_out.T,
    tablePhic                  = tablePhic_comp_mc,
    tableEta                   = tableEta_comp_mc,
    tablePR                    = tablePR_comp_mc,
    Table                      = ThermoPower.Choices.TurboMachinery.TableTypes.matrix,
    Ndesign                    = cfg_comp.N,
    explicitIsentropicEnthalpy = false,
    Tdes_in                    = cfg_comp.st_in.T)
    annotation(
    Placement(visible = true, transformation(origin = {103, -11}, extent = {{-11, -11}, {11, 11}}, rotation = 0)));
 
  Steps.TPComponents.GasStateReader r_comp_in(redeclare package Medium = Medium); 
  
  Steps.TPComponents.FixedPController FPC_comp(
    redeclare package Medium = Medium,
    tablePhic                  = tablePhic,
    tablePR                    = tablePR,
    Table                      = ThermoPower.Choices.TurboMachinery.TableTypes.matrix,
    Ndesign                    = cfg_comp.N,
    Tdes_in                    = cfg_comp.st_in.T
  );
  
  Modelica.Mechanics.Rotational.Sources.Speed speed_comp(
    exact = false,
    w_ref(nominal(cfg_comp.N))
  );

/*
  Modelica.Mechanics.Rotational.Sources.ConstantSpeed const_speed_comp(
      w_fixed=N_s_compressor, useSupport=false) annotation(
    Placement(visible = true, transformation(origin = {81, -7}, extent = {{-5, -5}, {5, 5}}, rotation = 0)));  
*/

  ThermoPower.Gas.Turbine turbine(
  redeclare package Medium = medium_main, 
  fileName                   = Modelica.Utilities.Files.loadResource("modelica://Steps/Resources/Data/turbine_map.txt"),
  tablePhic                  = fill(0.0, 14, 12),                                                                          //tablePhic, un
  tableEta                   = fill(0.0, 14, 12),                                                                          //tableEta, 
  pstart_in                  = cfg_turb.st_in.p,
  pstart_out                 = cfg_turb.st_out.p,
  Tstart_in                  = cfg_turb.st_in.T,
  Tstart_out                 = cfg_turb.st_out.T,
  Ndesign                    = cfg_turb.N,
  Tdes_in                    = cfg_turb.st_in.T,
  Table                      = ThermoPower.Choices.TurboMachinery.TableTypes.file,
  explicitIsentropicEnthalpy = false,
  gas_in(
    p(nominal = cfg_turb.st_in.p), 
    T(nominal = cfg_turb.st_in.T)),
  gas_iso(
    p(nominal = cfg_turb.st_out.p), 
    T(nominal = cfg_turb.st_out.T),
    h(start = cfg_turb.st_out.h))) annotation(
    Placement(visible = true, transformation(origin = {-69, -51}, extent = {{-11, -11}, {11, 11}}, rotation = 0))); 
 
  Steps.TPComponents.GasStateReader r_turb_in(redeclare package Medium = Medium) annotation(
    Placement(visible = true, transformation(origin = {-92, -38}, extent = {{-4, -4}, {4, 4}}, rotation = 0)));

  Steps.TPComponents.TurbineFixedPController FPC_turb(
    redeclare package Medium   = Medium,
    fileName  = Modelica.Utilities.Files.loadResource("modelica://Steps/Resources/Data/turbine_map_2.txt"),
    tablePhic = fill(0.0, 14, 12), 
    Table     = ThermoPower.Choices.TurboMachinery.TableTypes.file,
    Ndesign   = cfg_turb.N,
    Tdes_in   = cfg_turb.st_in.T,
    use_in_p0 = true,    
    in_T1(start = cfg_turb.st_in.T, nominal = cfg_turb.st_in.T),
    in_p1(start = cfg_turb.st_in.p, nominal = cfg_turb.st_in.p),
    in_w1(start = cfg_turb.st_in.mdot, nominal = cfg_turb.st_in.mdot)
  ) "Fixed pressure controller for main compressor";
  
  Modelica.Mechanics.Rotational.Sources.ConstantSpeed const_speed_turb(
      w_fixed=cfg_turb.N, useSupport=false) annotation(
    Placement(visible = true, transformation(origin = {-49, -47}, extent = {{5, -5}, {-5, 5}}, rotation = 0)));
   
  
  Modelica.Mechanics.Rotational.Sources.Speed speed_turb(
    exact = false,
    w_ref(start = cfg_turb.N, nominal= cfg_turb.N),
    w(start = cfg_turb.N, nominal = cfg_turb.N)
  );  
 
  // value to accelerate the simulation
  constant Integer SEC_PER_HOUR = integer(60 * 60 / time_scaling_factor); 
  constant Real time_scaling_factor = 12; // 5 min = 1 hour
  // constant Real time_scaling_factor = 1; // 1 hour = 1 hour  

  // ramp input to simulate the ramp change in ASPEN+ simulation
  // ramp change for case I in Guan's report
  Steps.TPComponents.RampSeq_W ramp_mdot(    
    time_start = 3 * SEC_PER_HOUR,
    interval   = 1 * SEC_PER_HOUR,
    duration   = 0.15 * SEC_PER_HOUR,
    offset     = st_source.mdot
  );
  
  Modelica.Blocks.Sources.Constant p_constant_source(k(start = st_source.p, nominal = st_source.p));
  
  Modelica.Blocks.Sources.Constant p_constant_sink(k(start = st_sink.p, nominal = st_sink.p));

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
    10,0.725,0.745,0.762,0.771,0.779,0.796;
    11,0.671,0.682,0.699,0.724,0.751,0.766;
    12,0.624,0.641,0.657,0.672,0.69,0.71;
    13,0.564,0.592,0.61,0.637,0.65,0.665;
    14,0.443,0.493,0.535,0.588,0.6,0.627;
    15,0.361,0.39,0.424,0.472,0.526,0.566
   ];

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
    10,2.61e-05,2.73e-05,2.85e-05,2.97e-05,3.09e-05,3.21e-05;
    11,2.73e-05,2.85e-05,2.97e-05,3.09e-05,3.21e-05,3.33e-05;
    12,2.85e-05,2.97e-05,3.09e-05,3.21e-05,3.33e-05,3.45e-05;
    13,2.97e-05,3.09e-05,3.21e-05,3.33e-05,3.45e-05,3.56e-05;
    14,3.09e-05,3.21e-05,3.33e-05,3.45e-05,3.56e-05,3.68e-05;
    15,3.21e-05,3.33e-05,3.45e-05,3.56e-05,3.68e-05,3.8e-05
  ];

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
    10,1.45,1.53,1.62,1.72,1.83,1.97;
    11,1.38,1.45,1.54,1.64,1.76,1.89;
    12,1.32,1.39,1.47,1.55,1.65,1.76;
    13,1.25,1.32,1.39,1.48,1.56,1.66;
    14,1.17,1.23,1.3,1.39,1.46,1.58;
    15,1.11,1.15,1.2,1.27,1.36,1.45
  ];  

initial equation
  Inertia1.w = N_s_compressor;
  
equation

  // heater + compressor + turbine
  connect(source_cold.flange, T_waterIn.inlet);
  connect(T_waterIn.outlet, compressor.inlet);  
  connect(compressor.outlet, HE.waterIn);
  connect(HE.waterOut, turbine.inlet);
  connect(turbine.outlet, T_waterOut.inlet);    
  connect(T_waterOut.outlet, sink_cold.flange);    

  connect(tt_mdot_hot_in.y, source_hot.in_w0);
  connect(source_hot.flange, T_gasIn.inlet);
  connect(T_gasIn.outlet, HE.gasIn) annotation(
        Line(points = {{-50, 0}, {-20, 0}}, color = {159, 159, 223}, thickness = 0.5, smooth = Smooth.None)); 
  connect(HE.gasOut, T_gasOut.inlet) annotation(
    Line(points = {{34, 0}, {34, 0}, {20, 0}}, color = {159, 159, 223}, thickness = 0.5));
  connect(T_gasOut.outlet, sink_hot.flange) annotation(
    Line(points = {{46, 0}, {46, 0}, {60, 0}}, color = {159, 159, 223}, thickness = 0.5));     

  connect(tt_T_motor.y, torque_comp.tau);
  connect(torque_comp.flange, Inertia1.flange_b);
  connect(Inertia1.flange_a, compressor.shaft_b) annotation(
    Line(points = {{86, -6}, {96, -6}, {96, -10}, {96, -10}}));    

  connect(const_speed_turbine.flange, turbine.shaft_a);

annotation(
    Diagram(graphics),
    // for steady-state simulation - value check
    experiment(StartTime = 0, StopTime = 300, Tolerance = 1e-3, Interval = 5),
    // for complete transient simulation
    // experiment(StartTime = 0, StopTime = 600, Tolerance = 1e-3, Interval = 10),
    __OpenModelica_commandLineOptions = "--matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection -d=initialization,NLSanalyticJacobian,aliasConflicts",
    // remove the option flag --matchingAlgorithm=PFPlusExt, which may lead to 'Internal error - IndexReduction.dynamicStateSelectionWork failed!' during Translation
    // __OpenModelica_commandLineOptions = "--matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection -d=initialization,NLSanalyticJacobian,aliasConflicts",    
    __OpenModelica_simulationFlags(lv = "LOG_DEBUG,LOG_NLS,LOG_NLS_V,LOG_STATS,LOG_INIT,LOG_STDOUT, -w", outputFormat = "mat", s = "dassl", nls = "homotopy")
    );
end TPDyn_SimpleCycle_fixed_p;
