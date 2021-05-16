within Steps.TPComponents;

model PCHE "PCHE model based on Thermo Power"
  extends Interfaces.HeatExchangerG2G;
    
  import SI = Modelica.SIunits;
  import Gas = ThermoPower.Gas;  
  import Choices = ThermoPower.Choices;
  import Thermal = ThermoPower.Thermal;
  import Options = ThermoPower.Choices.Init.Options;
  import Steps.Components.KimCorrelations;
  import Steps.Components.MaterialConductivity;  
  import Steps.Model.{EntityConfig, EntityGeoParam, EntityThermoParam, HEBoundaryCondition};
  
  replaceable model HeatTransfer_F                           = Thermal.HeatTransferFV.IdealHeatTransfer constrainedby ThermoPower.Thermal.BaseClasses.DistributedHeatTransferFV annotation(choicesAllMatching = true);
  replaceable model HeatTransfer_G                           = Thermal.HeatTransferFV.IdealHeatTransfer constrainedby ThermoPower.Thermal.BaseClasses.DistributedHeatTransferFV annotation(choicesAllMatching = true);
  replaceable model HeatExchangerTopology                    = Thermal.HeatExchangerTopologies.CounterCurrentFlow constrainedby ThermoPower.Thermal.BaseClasses.HeatExchangerTopologyData annotation(choicesAllMatching = true);
  parameter   Choices.Flow1D.FFtypes FFtype_G                = ThermoPower.Choices.Flow1D.FFtypes.NoFriction "Friction Factor Type, gas side";
  parameter   Real Kfnom_G                                   = 0 "Nominal hydraulic resistance coefficient, gas side";
  parameter   SI.Pressure dpnom_G                            = 0 "Nominal pressure drop, gas side (friction term only!)";
  parameter   SI.Density rhonom_G                            = 0 "Nominal inlet density, gas side";
  parameter   Real Cfnom_G                                   = 0 "Nominal Fanning friction factor, gas side";
  parameter   Choices.Flow1D.FFtypes FFtype_F                = ThermoPower.Choices.Flow1D.FFtypes.NoFriction "Friction Factor Type, fluid side";
  parameter   Real Kfnom_F                                   = 0 "Nominal hydraulic resistance coefficient";
  parameter   SI.Pressure dpnom_F                            = 0 "Nominal pressure drop fluid side (friction term only!)";
  parameter   SI.Density rhonom_F                            = 0 "Nominal inlet density fluid side";
  parameter   SI.PerUnit Cfnom_F                             = 0 "Nominal Fanning friction factor";
  parameter   Choices.Flow1D.HCtypes HCtype_F                = ThermoPower.Choices.Flow1D.HCtypes.Downstream "Location of the hydraulic capacitance, fluid side";
  parameter   Boolean gasQuasiStatic                         = false "Quasi-static model of the flue gas (mass, energy and momentum static balances";
  parameter   Boolean fluidQuasiStatic                       = false "Quasi-static model of the fluid (mass, energy and momentum static balances";
  constant    Real pi                                        = Modelica.Constants.pi;
  // parameter   SI.Distance L                                  = 1 "Tube length";
  parameter   Choices.FluidPhase.FluidPhases FluidPhaseStart = Choices.FluidPhase.FluidPhases.Liquid "Fluid phase (only for initialization!)" annotation(Dialog(tab = "Initialization"));
  parameter   Boolean SSInit                                 = false "Steady State initialization";
  
  // Thermal conductivity of LTR's metal wall 
  // material inconel_750
  // parameter Real table_k_metalwall[:, :] = [149, 16.9; 316, 20.5; 538, 26.5; 649, 28.7; 760, 31.4; 871, 35.3];  
  
  //IMPORTANT:
  // To date, (because of a bug in Gas.Flow1DFV
  // use TPComponents.GasFlow1DFV for transient simulation when two PCHEs connected (LTR and HTR),
  // use GasFlow1DFV otherwise, which is much quicker due to the property calculation
  
  // Gas.Flow1DFV fluidFlow(
  TPComponents.GasFlow1DFV fluidFlow(
  redeclare package Medium     = FluidMedium,
  redeclare model HeatTransfer = HeatTransfer_F,
  Nt          = Nt,                                                       //1, 
  N           = N_F,
  Nw          = Nw_F,
  wnom        = fluidNomFlowRate,
  initOpt     = if SSInit then Options.steadyState else Options.noInit,
  QuasiStatic = fluidQuasiStatic,
  pstart      = pstart_F,
  L           = fluidLength,                                              // Should be L = exchSurface_F ^ 2 / (fluidVol * pi * 4), instead of fixed L = 1 
  A           = fluidVol / fluidLength,                                   // fluidVol is account for single tube
  omega       = exchSurface_F / fluidLength,                              // exchSurface_F is account for single tube
  Dhyd        = fluidVol*4 / exchSurface_F,
  FFtype      = FFtype_F,
  Kfnom       = Kfnom_F,
  dpnom       = dpnom_F,
  rhonom      = rhonom_F,
  Cfnom       = Cfnom_F,
  Tstartbar   = Tstartbar_F,
  Tstartin    = Tstartbar_F,                                              //bc.st_cold_in.T, 
  Tstartout   = Tstartbar_F) annotation(
    Placement(transformation(extent = {{-10, -66}, {10, -46}}, rotation = 0)));
    
  //changed Medium=FlueGasMedium to Medium=FluidMedium
  
  //IMPORTANT:
  // To date, (because of a bug in Gas.Flow1DFV
  // use TPComponents.GasFlow1DFV for transient simulation when two PCHEs connected (LTR and HTR),
  // use GasFlow1DFV otherwise, which is much quicker due to the property calculation  
  
  TPComponents.GasFlow1DFV gasFlow(
  // Gas.Flow1DFV gasFlow(
  redeclare package Medium = FlueGasMedium, 
  redeclare model HeatTransfer = HeatTransfer_G,
  Nt          = Nt,                                                       //1,  
  N           = N_G,
  Nw          = Nw_G,
  wnom        = gasNomFlowRate,                                           // wnom(total) of gasFlow is different from wnom(for single tube, = gasNomFlowRate /Nt) of HeatTransfer_G  
  initOpt     = if SSInit then Options.steadyState else Options.noInit,
  QuasiStatic = gasQuasiStatic,
  pstart      = pstart_G,
  L           = gasLength,                                                // Should be L = exchSurface_G ^ 2 / (gasVol * pi * 4), instead of fixed L = 1
  A           = gasVol / gasLength,                                       // gasVol is account for single tube, 
  omega       = exchSurface_G / gasLength,                                // exchSurface_G is account for single tube, 
  Dhyd        = gasVol*4 / exchSurface_G,
  FFtype      = FFtype_G,
  Kfnom       = Kfnom_G,
  dpnom       = dpnom_G,
  rhonom      = rhonom_G,
  Cfnom       = Cfnom_G,
  Tstartbar   = Tstartbar_G,
  Tstartin    = Tstartbar_G,                                              //bc.st_hot_in.T, 
  Tstartout   = Tstartbar_G)
   annotation(
    Placement(transformation(extent = {{-12, 66}, {12, 46}}, rotation = 0)));
  
  PCHEMetalWallFV metalWall(
    L                     = cfg_wall.geo_wall.l,
    r_c                   = cfg_wall.geo_area.d / 2,
    w_ch                  = cfg_wall.geo_area.w,
    h_ch                  = cfg_wall.geo_area.h,
    dx                    = cfg_wall.geo_area.p1,
    //L_wall                = cfg_wall.geo_wall.p1,
    Nw                    = Nw_F,
    Nt                    = Nt * 2,
    Tstartbar             = (Tstartbar_G + Tstartbar_F) / 2,
    Tstart1               = Tstartbar_G,                       //bc.st_hot_out.T, 
    TstartN               = Tstartbar_F,                       //bc.st_hot_in.T,   
    WallRes               = false,
    table_th_conductivity = cfg_wall.table_k,
    rhomcm                = rhomcm
  ) annotation(
    Placement(transformation(extent = {{-10, -24}, {10, -4}})));
  
  Thermal.HeatExchangerTopologyFV heatExchangerTopology(Nw = Nw_F, redeclare model HeatExchangerTopology = HeatExchangerTopology) annotation(
    Placement(transformation(extent = {{-10, 6}, {10, 26}})));
equation
  connect(gasFlow.infl, gasIn) annotation(
    Line(points = {{-12, 56}, {-100, 56}, {-100, 0}}, color = {159, 159, 223}, thickness = 0.5));
  connect(gasFlow.outfl, gasOut) annotation(
    Line(points = {{12, 56}, {100, 56}, {100, 0}}, color = {159, 159, 223}, thickness = 0.5));
  connect(fluidFlow.outfl, waterOut) annotation(
    Line(points = {{10, -56}, {40, -56}, {40, -100}, {0, -100}}, thickness = 0.5, color = {0, 0, 255}));
  connect(fluidFlow.infl, waterIn) annotation(
    Line(points = {{-10, -56}, {-40, -56}, {-40, 100}, {0, 100}}, thickness = 0.5, color = {0, 0, 255}));
  connect(metalWall.ext, fluidFlow.wall) annotation(
    Line(points = {{0, -17.1}, {0, -51}}, color = {255, 127, 0}, smooth = Smooth.None));
  connect(heatExchangerTopology.side2, metalWall.int) annotation(
    Line(points = {{0, 12.9}, {0, -11}}, color = {255, 127, 0}, smooth = Smooth.None));
  connect(heatExchangerTopology.side1, gasFlow.wall) annotation(
    Line(points = {{0, 19}, {0, 51}}, color = {255, 127, 0}, smooth = Smooth.None));
  annotation(
    Diagram(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}}), graphics));
end PCHE;
