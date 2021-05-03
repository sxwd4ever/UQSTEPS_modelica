within Steps.Cycle.Transient;

model TP_SimpleCycle
  "simple cycle built referring Bone's Study [bone2021] (with an extra cooler) comp-hx-turbine-cooler"

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
  // package medium_hot = Modelica.Media.IdealGases.SingleGases.CO2;
  package medium_hot = Steps.Media.ThermiaOilD;
  // package medium_cold = Modelica.Media.IdealGases.SingleGases.CO2;
  package medium_cold = Steps.Media.SCO2;  
  
  // geometry parameters
  constant Real pi = Modelica.Constants.pi;
  parameter Integer N_ch = integer(2400) "channel number";
  parameter Integer N_seg = 10 "number of segments in one tube";
  // parameter SI.Length D_ch = 2e-3 "channel diameter, semi circular tube";
  parameter SI.Length D_ch = 1.72e-3  "channel diameter, semi circular tube";
  parameter SI.Length r_ch = D_ch / 2 "channel radiaus";
  // parameter SI.Length L_fp = 270e-3 "channel flow path length";  
  // parameter SI.Length L_fp = (200 + 186) * 1e-3 "equivalent valid channel flow path length";  
  parameter SI.Length L_fp = 250 * 1e-3 "equivalent valid channel flow path length"; 
  parameter SI.Length L_pitch = 12e-3 "pitch length";
  parameter Real a_phi = 36 "pitch angle degree";
  parameter SI.Length H_ch = 4.17e-3 "Height of the solid domain, containing one cold tube and one hot tube";
  parameter SI.Length W_ch = 2.3e-3 "Width of the solid domain";
  parameter SI.Length T_wall = 0.51e-3 "Wall thinckness";
  parameter SI.Length L_wall = 420e-3 "Length of wall, not necessarily equals to length of flow path";
  parameter SI.Area A = pi * r_ch ^2 / 2 "Area of cross section of semi circular tube";

  // boundary conditon
  
  // zigzag higher T
  // parameter SI.Velocity u_hot_in = 7.564 "hot inlet velocity m/s";
  // parameter SI.Velocity u_cold_in = 1.876 "cold inlet velocity m/s";
  parameter SI.Pressure p_hot_in =  from_bar(10) "hot inlet pressure - Irrelevant for Incompressible Thermia Oil" ;
  parameter SI.Pressure p_cold_in = 12.567478e6 "cold inlet pressure";
  parameter SI.Temperature T_hot_in = from_degC(103.222748) "hot inlet temperature, K";
  parameter SI.Temperature T_hot_out = from_degC(96.145935) "cold outlet temperature, K";
  parameter SI.Temperature T_cold_in = from_degC(28.910231) "cold inlet temperature, K";
  parameter SI.Temperature T_cold_out = from_degC(99.666342) "cold outlet temperature, K";
  
  // pressure drop correction coefficient 
  // parameter Real kc_dp = 1.0;  
  /*
  parameter Real kc_cf_hot = 1;  
  parameter Real kc_cf_cold = 1;
  parameter Real Cf_a_hot = 1, Cf_b_hot = 1, Cf_c_hot = 1;
  parameter Real Cf_a_cold = 1, Cf_b_cold = 1, Cf_c_cold = 1;
  */
  
  parameter Real Cf_C1 = 1, Cf_C2 = 1, Cf_C3 = 1;
  //parameter Real Cf_a_cold = 1, Cf_b_cold = 1, Cf_c_cold = 1;  
  
  // meshram's cp and rho for alloy Inconel 617
  // parameter Modelica.SIunits.Density rho_wall = 8360 "density of wall, kg/m3";
  // parameter Modelica.SIunits.SpecificHeatCapacity cp_wall = 417 "cp of wall, J/kg-K";  
   
  // Stainless 316, 316L, 317, 317L
  parameter Modelica.SIunits.Density rho_wall = 8030 "density of wall, kg/m3";
  parameter Modelica.SIunits.SpecificHeatCapacity cp_wall = 485 "cp of wall, J/kg-K";  
  // thermal conductivity (T in K) https://www.theworldmaterial.com/aisi-316-ss316-stainless-steel-properties-composition/
  // parameter Real table_k_metalwall[:,:] = [20, 12.1; 100, 16.3; 500, 21.5];
  parameter Real table_k_metalwall[:,:] = [293.15, 12.1; 373.15, 16.3; 773.15, 21.5];

  // parameter SI.Density rho_hot_in = medium_hot.density_pT(p_hot_in, T_hot_in);
  // parameter SI.Density rho_cold_in = medium_cold.density_pT(p_cold_in, T_cold_in);
  parameter SI.MassFlowRate mdot_hot_in = 1.218944;
  parameter SI.MassFlowRate mdot_cold_in = 0.0854299999999999;

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
  //parameter HEBoundaryCondition bc_HE = cfg_test.bc_heater;

  //Components
  // for transient simulation, set initOpt = steadyState
  inner ThermoPower.System system(allowFlowReversal = false, initOpt = ThermoPower.Choices.Init.Options.noInit) annotation(
    Placement(transformation(extent = {{80, 80}, {100, 100}})));  
  
  parameter Model.SimpleCycleConfig cfg_test(
    redeclare package medium_main = medium_cold,
    redeclare package medium_heater = medium_hot,
    mdot_main = mdot_cold_in,
    mdot_heater = mdot_hot_in,     
    N_ch_h = N_ch,
    r_h = r_ch,    
    L_h = L_fp,    
    N_seg = N_seg,    
    pitch = L_pitch,
    phi = a_phi,
    T_heater_hin = T_hot_in,
    T_heater_hout = T_hot_out,
    T_heater_cin = T_cold_in,
    T_heater_cout = T_cold_out,
    p_heater_hin =  p_hot_in,
    p_comp_out = p_cold_in   
  );    

 
  Steps.TPComponents.PCHE heater(
    redeclare package FluidMedium = medium_cold, 
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
  
    bc = bc_HE,     
    // geo_hot = cfg.cfg_heater_hot.geo,
    // geo_cold = cfg.cfg_heater_cold.geo,
    // geo_tube = cfg.cfg_heater_tube.geo,  
    // thermo_hot = cfg.cfg_heater_hot.thermo,
    // thermo_cold = cfg.cfg_heater_cold.thermo,
    // thermo_tube = cfg.cfg_heater_tube.thermo, 
    geo_hot = cfg.cfg_LTR_hot.geo,
    geo_cold = cfg.cfg_LTR_cold.geo,
    geo_tube = cfg.cfg_LTR_tube.geo,  
    thermo_hot = cfg.cfg_LTR_hot.thermo,
    thermo_cold = cfg.cfg_LTR_cold.thermo,
    thermo_tube = cfg.cfg_LTR_tube.thermo,     
    table_k_metalwall =   table_k_metalwall,
    L = L_fp,
    SSInit = false,
    gasQuasiStatic = true,
    fluidQuasiStatic = true,
    metalWall(L = L_wall, w_ch = W_ch, h_ch = H_ch, dx = T_wall)
    // metalQuasiStatic = true
    // override the values of Am and L of metaltubeFV
    // to make them agree with semi-circular tube of PCHE
    // ('final' modifier of Am in metalTubeFv was removed as well)
    //metalTube(WallRes=false, L = 1, rhomcm=200, Am = HE.metalVol / 1) 
  )
  annotation(
    Placement(transformation(extent = {{-20, -20}, {20, 20}}, rotation = 0)));

/*
  Steps.TPComponents.PCHE heater(
    redeclare package FluidMedium = medium_cold, 
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
  
    bc = bc_HE, 
    geo_hot = cfg.cfg_LTR_hot.geo,
    geo_cold = cfg.cfg_LTR_cold.geo,
    geo_tube = cfg.cfg_LTR_tube.geo,  
    thermo_hot = cfg.cfg_LTR_hot.thermo,
    thermo_cold = cfg.cfg_LTR_cold.thermo,
    thermo_tube = cfg.cfg_LTR_tube.thermo, 
    table_k_metalwall =   table_k_metalwall,
    L = L_fp,
    SSInit = false,
    gasQuasiStatic = true,
    fluidQuasiStatic = true,
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
/*   
  ThermoPower.Gas.Compressor compressor(
    redeclare package Medium = medium_main,
    pstart_in=cfg.st_comp_in.p,
    pstart_out=cfg.st_comp_out.p,
    Tstart_in=cfg.st_comp_in.T,
    Tstart_out=cfg.st_comp_out.T,
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
  */

  ThermoPower.Gas.SourceMassFlow source(
    redeclare package Medium = medium_cold, 
    T = bc_HE.st_cold_in.T, 
    p0 = bc_HE.st_cold_in.p, 
    use_in_T = false, 
    w0 = bc_HE.st_cold_in.mdot,
    gas(p(nominal = bc_HE.st_cold_in.p), 
    T(nominal=bc_HE.st_cold_in.T))) 
  annotation(
    Placement(transformation(origin = {0, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 270)));
  
  ThermoPower.Gas.SinkPressure sink(
    redeclare package Medium = medium_cold, 
    p0 = bc_HE.st_cold_out.p, 
    T = bc_HE.st_cold_out.T) 
  annotation(
    Placement(transformation(origin = {0, -80}, extent = {{-10, -10}, {10, 10}}, rotation = 270)));

  ThermoPower.Gas.SinkPressure sink_heater_hot(
    redeclare package Medium = medium_hot,
    T = bc_HE.st_hot_out.T, 
    p0 = bc_HE.st_hot_out.p) 
  annotation(
    Placement(transformation(extent = {{60, -10}, {80, 10}}, rotation = 0)));
  
  ThermoPower.Gas.SourceMassFlow source_heater_hot(
    redeclare package Medium = medium_hot, 
    T = bc_HE.st_hot_in.T, 
    p0 = bc_HE.st_hot_in.p, 
    w0 = bc_HE.st_hot_in.mdot,
    use_in_T = false,
    use_in_w0 = false,
    gas(p(nominal = bc_HE.st_hot_in.p), 
    T(nominal=bc_HE.st_hot_in.T))) 
  annotation(
    Placement(transformation(extent = {{-70, -10}, {-50, 10}}, rotation = 0))); 

/*
  TPComponents.GasStateReader r01(redeclare package Medium = medium_main) annotation(
    Placement(visible = true, transformation(origin = {66, -56}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
*/  
  TPComponents.GasStateReader r02(redeclare package Medium = medium_cold) annotation(
    Placement(visible = true, transformation(origin = {96, 44}, extent = {{-6, -6}, {6, 6}}, rotation = 180)));
  TPComponents.GasStateReader r03(redeclare package Medium = medium_cold) annotation(
    Placement(visible = true, transformation(origin = {-22, 50}, extent = {{-6, -6}, {6, 6}}, rotation = -90)));
/*
  TPComponents.GasStateReader r04(redeclare package Medium = medium_main) annotation(
    Placement(visible = true, transformation(origin = {-32, 20}, extent = {{-6, -6}, {6, 6}}, rotation = 180)));
*/

  TPComponents.WaterStateReader rh1(redeclare package Medium = medium_hot) annotation(
    Placement(visible = true, transformation(origin = {-66, 22}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
  TPComponents.WaterStateReader rh2(redeclare package Medium = medium_hot) annotation(
    Placement(visible = true, transformation(origin = {-64, 64}, extent = {{-6, -6}, {6, 6}}, rotation = 180)));

  
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
/*
  // open loop with 
  // main path: source -> compressor -> heater -> turbine -> sink
  connect(source.flange, compressor.inlet) annotation(
    Line(points = {{70, -56}, {72, -56}, {72, 4}, {95, 4}, {95, -2}, {94, -2}}, color = {159, 159, 223}));
  connect(compressor.outlet, heater.gasIn) annotation(
    Line(points = {{-37, 20}, {-50, 20}, {-50, 41}, {-54, 41}}, color = {0, 0, 255}, thickness = 0.5));
  connect(heater.gasOut, turbine.inlet);
  connect(turbine.outlet, sink.flange) annotation(
    Line(points = {{32, 14}, {14, 14}, {14, 2}, {14, 2}}, color = {159, 159, 223}));
  
  // turbine, compressor
  connect(const_speed_turbine.flange, turbine.shaft_a) annotation(
    Line(points = {{-62, -51}, {-51, -51}, {-51, -50}, {-52, -50}}));
  connect(const_speed_comp.flange, compressor.shaft_a) annotation(
    Line(points = {{86, -6}, {96, -6}, {96, -10}, {96, -10}}));
  
  // hot side of heater
  connect(source_heater_hot.flange, heater.waterIn) annotation(
    Line(points = {{-62, 22}, {-61, 22}, {-61, 34}, {-60, 34}}, color = {255, 0, 0}, thickness = 0.5));
  connect(heater.waterOut, sink_heater_hot.flange) annotation(
    Line(points = {{-68, 64}, {-76, 64}}, color = {255, 0, 0}, thickness = 0.5));
*/

  // heater alone
  connect(source.flange, r02.inlet);
  connect(r02.outlet, heater.gasIn);
  connect(heater.gasOut, r03.inlet);
  connect(r03.outlet, sink.flange);
     
  connect(source_heater_hot.flange, rh1.inlet);
  connect(rh1.inlet, heater.waterIn);
  connect(heater.waterOut, rh2.inlet);
  connect(rh2.inlet, sink_heater_hot.flange);
   
  annotation(
    Diagram(coordinateSystem(extent = {{-100, -100}, {120, 100}})),
    experiment(StartTime = 0, StopTime = 100, Tolerance = 1e-3, Interval = 10),
    __OpenModelica_commandLineOptions = "--matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection -d=initialization,NLSanalyticJacobian,bltdump",
    __OpenModelica_simulationFlags(lv = "LOG_DEBUG,LOG_NLS,LOG_NLS_V,LOG_STATS,LOG_INIT,LOG_STDOUT, -w", outputFormat = "mat", s = "dassl", nls = "homotopy"));
end TP_SimpleCycle;
