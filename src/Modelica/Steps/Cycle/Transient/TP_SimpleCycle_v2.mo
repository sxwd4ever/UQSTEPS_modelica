within Steps.Cycle.Transient;

model TP_SimpleCycle_v2
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
  // package medium_cold = Modelica.Media.IdealGases.SingleGases.CO2;
  package medium_cold = Steps.Media.SCO2;
  
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

  // parameter SI.Density rho_hot_in = medium_hot.density_pT(p_hot_in, T_hot_in);
  // parameter SI.Density rho_cold_in = medium_cold.density_pT(p_cold_in, T_cold_in);
  parameter SI.MassFlowRate mdot_heater = 10;    //1.218944;
  parameter SI.MassFlowRate mdot_main   = 10.5;  // 0.0854299999999999;
  
  parameter Real N_s_compressor = 2100 "rotational speed of compressor";

  // use configuration of LTR for this test since the mdot are different for hot and cold side
  parameter Model.PBConfig_PCHE cfg_test(
    p_pump_in         = p_source,
    p_ATM             = p_source,
    p_pump_out        = p_comp_out,
    p_heater          = p_hot_in,
    mdot_main         = mdot_main,
    mdot_pump         = mdot_main,
    mdot_heater       = mdot_heater,
    T_HTR_cold_out    = T_heater_cin,
    T_heater_cold_out = T_heater_cout,
    T_heater_hot_in   = T_heater_hin,    // T_LTR_hot_in = T_HTR_hot_out,
    T_heater_hot_out  = T_heater_hout,
    T_cooler_cold_out = T_source,
    T_HTR_hot_in      = T_sink,
    r_heater          = r_ch,
    L_h               = L_fp,
    N_ch_h            = N_ch,
    N_seg_heater      = N_seg,
    pitch             = L_pitch,
    phi               = a_phi,
    rho_wall          = rho_wall,
    cp_wall           = cp_wall
  ); 

  // local reference of the config objects
  parameter Model.SimpleCycleConfig cfg;
  parameter Model.HeatExchangerConfig cfg_heater = cfg.cfg_heater;
  parameter Model.TurbomachineryConfig cfg_comp = cfg.cfg_comp;
  parameter Model.TurbomachineryConfig cfg_turb = cfg.cfg_turb;
  
 
  //Components
  // for transient simulation, set initOpt = steadyState
  inner ThermoPower.System system(allowFlowReversal = false, initOpt = ThermoPower.Choices.Init.Options.noInit) annotation(
    Placement(transformation(extent = {{80, 80}, {100, 100}})));  
  
  ThermoPower.Gas.SourceMassFlow source_cold(
    redeclare package Medium = medium_cold,
    use_in_T = false,
    T        = cfg.st_source.T,
    p0       = cfg.st_source.p,
    w0       = cfg.st_source.mdot,
    // T = 322,
    // p0 = 8.65e6,
    // w0 = 10.5,
    // gas(p(nominal = bc_HE.st_cold_in.p), 
    gas(
      p(nominal = cfg.st_source.p), 
      // p(nominal = 8.65e6), 
      // T(nominal=bc_HE.st_cold_in.T))) 
      T(nominal = cfg.st_source.T))
      // T(start = 322, nominal=322))
    ) 
  annotation(
    Placement(transformation(origin = {0, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 270)));
  
  ThermoPower.Gas.SinkPressure sink_cold(
    redeclare package Medium = medium_cold,
    // p0 = 12e6, 
    // T = 348,
    p0 = cfg.st_sink.p,
    T  = cfg.st_sink.T,
    gas(
      p(nominal = cfg.st_sink.p), 
      T(nominal= cfg.st_sink.T))
      // T(start = 348, nominal=348))    
  );

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
    use_in_w0 = false,
    gas(
      p(nominal = cfg.st_heater_hin.p), 
      T(nominal=cfg.st_heater_hin.T))) 
  annotation(
    Placement(transformation(extent = {{-70, -10}, {-50, 10}}, rotation = 0))); 


  ThermoPower.Gas.SensT T_waterIn(redeclare package Medium = medium_cold) annotation(
    Placement(transformation(origin = {4, -50}, extent = {{-10, -10}, {10, 10}}, rotation = 270)));
  ThermoPower.Gas.SensT T_waterOut(redeclare package Medium = medium_cold) annotation(
    Placement(transformation(origin = {4, -50}, extent = {{-10, -10}, {10, 10}}, rotation = 270)));


  ThermoPower.Gas.SensT T_gasIn(redeclare package Medium = medium_hot);
  ThermoPower.Gas.SensT T_gasOut(redeclare package Medium = medium_hot);


  Steps.TPComponents.PCHE HE(
    redeclare package FluidMedium   = medium_cold,
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
    SSInit            = false,
    gasQuasiStatic    = true,
    fluidQuasiStatic  = true,
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
    redeclare package Medium = medium_cold,
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
    Tdes_in    = compressor.Tstart_in)
    annotation(
    Placement(visible = true, transformation(origin = {103, -11}, extent = {{-11, -11}, {11, 11}}, rotation = 0)));

  Modelica.Mechanics.Rotational.Sources.ConstantSpeed const_speed_comp(
      w_fixed=N_s_compressor, useSupport=false) annotation(
    Placement(visible = true, transformation(origin = {81, -7}, extent = {{-5, -5}, {5, 5}}, rotation = 0)));  
  /*
  // variable for validation
  Modelica.SIunits.Power Q_out    = (HE.gasIn.h_outflow - HE.gasOut.h_outflow) * HE.gasIn.m_flow;
  Modelica.SIunits.Power Q_in     = (HE.waterOut.h_outflow - HE.waterIn.h_outflow) * HE.waterIn.m_flow;
  Boolean                isQMatch = abs(Q_out - Q_in) < 1e-3;
  
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
  
  Modelica.Blocks.Math.MultiSwitch multiSwitch1(nu=2, expr = {mdot_hot_in, mdot_hot_in}) annotation(
    Placement(visible = true, transformation(origin = {-96, 80}, extent = {{-10, -10}, {30, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.BooleanStep booleanStep(startTime = 220, startValue = true) annotation(
    Placement(visible = true, transformation(origin = {-148, 98}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.BooleanStep booleanStep1(startTime = 480) annotation(
    Placement(visible = true, transformation(origin = {-150, 58}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  
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

  // performance map for re compressor      
  parameter Real tableEta_comp_rc[5, 4]  = [0, 95, 100, 105; 1, 0.85310219, 0.837591241, 0.832420925; 2, 0.868613139, 0.857238443, 0.851034063; 3, 0.860340633, 0.85, 0.842761557; 4, 0.85310219, 0.839659367, 0.816909976];
  parameter Real tablePhic_comp_rc[5, 4] = [0, 95, 100, 105; 1, 7.17663E-05, 8.05731E-05, 8.76935E-05; 2, 7.36401E-05, 8.20721E-05, 8.97547E-05; 3, 7.6076E-05, 8.46954E-05, 9.06916E-05; 4, 7.79498E-05, 8.63819E-05, 9.16285E-05];
  parameter Real tablePR_comp_rc[5, 4]   = [0, 95, 100, 105; 1, 1.967529638, 2.350588505, 2.785882673; 2, 1.915294338, 2.315764972, 2.681412073; 3, 1.810823737, 2.220000255, 2.524706172; 4, 1.654117837, 2.115529655, 2.359294389];
  
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

/*
  // compressor alone
//  connect(source_cold.flange, T_waterIn.inlet);
//  connect(T_waterIn.outlet, compressor.inlet);  
  
  connect(source_cold.flange, compressor.inlet);  
  
  
//   connect(compressor.outlet, T_waterOut.inlet);
//   connect(T_waterOut.outlet, sink_cold.flange);
  
  connect(compressor.outlet, sink_cold.flange);
  
  connect(const_speed_comp.flange, compressor.shaft_a) annotation(
    Line(points = {{86, -6}, {96, -6}, {96, -10}, {96, -10}}));

*/


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



  /*
  // temperature input1
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
  connect(booleanStep.y, multiSwitch1.u[1]) annotation(
    Line(points = {{-137, 98}, {-122.5, 98}, {-122.5, 80}, {-106, 80}}, color = {255, 0, 255}));
  connect(booleanStep1.y, multiSwitch1.u[2]) annotation(
    Line(points = {{-139, 58}, {-106, 58}, {-106, 80}}, color = {255, 0, 255}));
  connect(multiSwitch1.y, source_hot.in_w0) annotation(
    Line(points = {{-64, 80}, {-34, 80}, {-34, 82}, {-34, 82}}, color = {255, 127, 0}));
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
    // for steady-state simulation - value check
    experiment(StartTime = 0, StopTime = 2, Tolerance = 1e-3, Interval = 1),
    // for complete transient simulation
    // experiment(StartTime = 0, StopTime = 600, Tolerance = 1e-3, Interval = 10),
    __OpenModelica_commandLineOptions = "--matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection -d=initialization,NLSanalyticJacobian,aliasConflicts",
    // remove the option flag --matchingAlgorithm=PFPlusExt, which may lead to 'Internal error - IndexReduction.dynamicStateSelectionWork failed!' during Translation
    // __OpenModelica_commandLineOptions = "--matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection -d=initialization,NLSanalyticJacobian,aliasConflicts,bltdump",    
    __OpenModelica_simulationFlags(lv = "LOG_DEBUG,LOG_NLS,LOG_NLS_V,LOG_STATS,LOG_INIT,LOG_STDOUT, -w", outputFormat = "mat", s = "dassl", nls = "homotopy")
    );
end TP_SimpleCycle_v2;
