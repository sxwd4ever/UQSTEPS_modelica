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
  
  replaceable model HeatTransfer_F = Thermal.HeatTransferFV.IdealHeatTransfer constrainedby ThermoPower.Thermal.BaseClasses.DistributedHeatTransferFV annotation(
     choicesAllMatching = true);
  replaceable model HeatTransfer_G = Thermal.HeatTransferFV.IdealHeatTransfer constrainedby ThermoPower.Thermal.BaseClasses.DistributedHeatTransferFV annotation(
     choicesAllMatching = true);
  replaceable model HeatExchangerTopology = Thermal.HeatExchangerTopologies.CounterCurrentFlow constrainedby ThermoPower.Thermal.BaseClasses.HeatExchangerTopologyData annotation(
     choicesAllMatching = true);
  parameter Choices.Flow1D.FFtypes FFtype_G = ThermoPower.Choices.Flow1D.FFtypes.NoFriction "Friction Factor Type, gas side"; 
  
  parameter Real Kfnom_G = 0 "Nominal hydraulic resistance coefficient, gas side";
  parameter SI.Pressure dpnom_G = 0 "Nominal pressure drop, gas side (friction term only!)";
  parameter SI.Density rhonom_G = 0 "Nominal inlet density, gas side";
  parameter Real Cfnom_G = 0 "Nominal Fanning friction factor, gas side";
  parameter Choices.Flow1D.FFtypes FFtype_F = ThermoPower.Choices.Flow1D.FFtypes.NoFriction "Friction Factor Type, fluid side";
  parameter Real Kfnom_F = 0 "Nominal hydraulic resistance coefficient";
  parameter SI.Pressure dpnom_F = 0 "Nominal pressure drop fluid side (friction term only!)";
  parameter SI.Density rhonom_F = 0 "Nominal inlet density fluid side";
  parameter SI.PerUnit Cfnom_F = 0 "Nominal Fanning friction factor";
  parameter Choices.Flow1D.HCtypes HCtype_F = ThermoPower.Choices.Flow1D.HCtypes.Downstream "Location of the hydraulic capacitance, fluid side";
  parameter Boolean gasQuasiStatic = false "Quasi-static model of the flue gas (mass, energy and momentum static balances";
  parameter Boolean fluidQuasiStatic = false "Quasi-static model of the fluid (mass, energy and momentum static balances";
  // parameter Boolean metalQuasiStatic = false "Quasi-static model of the metalWall (energy static balances)";
  constant Real pi = Modelica.Constants.pi;
  parameter SI.Distance L = 1 "Tube length";
  parameter Choices.FluidPhase.FluidPhases FluidPhaseStart = Choices.FluidPhase.FluidPhases.Liquid "Fluid phase (only for initialization!)" annotation(
    Dialog(tab = "Initialization"));
  parameter Boolean SSInit = false "Steady State initialization";  
  parameter Real [:, :] table_k_metalwall;  
  
  //MaterialConductivity mc(name_material = name_material);  
  /*
  Gas.Flow1DFV fluidFlow(
   
  redeclare package Medium = FluidMedium, 
  redeclare model HeatTransfer = HeatTransfer_F, 
  A = (fluidVol * 4 / exchSurface_F) ^ 2 / 4 * pi, 
  Cfnom = Cfnom_F, 
  Dhyd = fluidVol * 4 / exchSurface_F, 
  FFtype = FFtype_F, 
  HydraulicCapacitance = HCtype_F, 
  Kfnom = Kfnom_F, 
  L = exchSurface_F ^ 2 / (fluidVol * pi * 4), F
  N = N_F,Nt = Nt, 
  Nw = Nw_F, 
  dpnom = dpnom_F, 
  initOpt = if SSInit then Options.steadyState else Options.noInit, 
  omega = fluidVol * 4 / exchSurface_F * pi, 
  pstart = pstart_F, 
  rhonom = rhonom_F, 
  wnom = fluidNomFlowRate) annotation(
    Placement(transformation(extent = {{-10, -66}, {10, -46}}, rotation = 0)));
  */
  Gas.Flow1DFV fluidFlow(
  Nt = Nt, //1, 
  N = N_F, 
  Nw = Nw_F,
  wnom = fluidNomFlowRate, 
  //initOpt = if SSInit then Options.steadyState else Options.noInit, 
  redeclare package Medium = FluidMedium, 
  initOpt = if SSInit then Options.steadyState else Options.noInit,
  QuasiStatic = fluidQuasiStatic, 
  pstart = pstart_F, 
  L = L, // Should be L = exchSurface_F ^ 2 / (fluidVol * pi * 4), instead of fixed L = 1 
  A = fluidVol / L, // fluidVol is account for single tube
  omega = exchSurface_F / L,// exchSurface_F is account for single tube
  Dhyd = fluidVol*4 / exchSurface_F,
  FFtype = FFtype_F, 
  Kfnom = Kfnom_F, 
  dpnom = dpnom_F, 
  rhonom = rhonom_F, 
  Cfnom = Cfnom_F, 
  Tstartbar = Tstartbar_F, 
  Tstartin = Tstartbar_F, //bc.st_cold_in.T, 
  Tstartout = Tstartbar_F, //bc.st_cold_out.T, 
  redeclare model HeatTransfer = HeatTransfer_F) annotation(
    Placement(transformation(extent = {{-10, -66}, {10, -46}}, rotation = 0)));
    
  //changed Medium=FlueGasMedium to Medium=FluidMedium
  Gas.Flow1DFV gasFlow(
  Nt = Nt, //1, 
  N = N_G, 
  Nw = Nw_G,
  wnom = gasNomFlowRate,  // wnom(total) of gasFlow is different from wnom(for single tube, = gasNomFlowRate /Nt) of HeatTransfer_G
  //initOpt = if SSInit then Options.steadyState else Options.noInit, 
  redeclare package Medium = FlueGasMedium, 
  initOpt = if SSInit then Options.steadyState else Options.noInit,
  QuasiStatic = gasQuasiStatic, 
  pstart = pstart_G, 
  L = L, // Should be L = exchSurface_G ^ 2 / (gasVol * pi * 4), instead of fixed L = 1
  A = gasVol / L, // gasVol is account for single tube, 
  omega = exchSurface_G / L, // exchSurface_G is account for single tube, 
  Dhyd = gasVol*4 / exchSurface_G,
  FFtype = FFtype_G, 
  Kfnom = Kfnom_G, 
  dpnom = dpnom_G, 
  rhonom = rhonom_G, 
  Cfnom = Cfnom_G, 
  Tstartbar = Tstartbar_G, 
  Tstartin = Tstartbar_G, //bc.st_hot_in.T, 
  Tstartout = Tstartbar_G, //bc.st_hot_out.T,   
  redeclare model HeatTransfer = HeatTransfer_G) annotation(
    Placement(transformation(extent = {{-12, 66}, {12, 46}}, rotation = 0)));
  /*
  Thermal.MetalTubeFV metalWall(
  L = exchSurface_F ^ 2 / (fluidVol * pi * 4),  
  Nw = Nw_F, 
  Tstartbar = Tstartbar_M, 
  WallRes = false, 
  lambda = lambda, 
  rext = (metalVol + fluidVol) * 4 / extSurfaceTub / 2, 
  rhomcm = rhomcm, 
  rint = fluidVol * 4 / exchSurface_F / 2) annotation(
    Placement(transformation(extent = {{-10, -24}, {10, -4}})));
  */
  
  PCHEMetalWallFV metalWall(
    L = L,  
    r_c = geo_tube.d / 2, 
    Nw = Nw_F,
    Nt = Nt * 2, 
    Tstartbar = (Tstartbar_G + Tstartbar_F) / 2, 
    Tstart1 =  Tstartbar_G, //bc.st_hot_out.T, 
    TstartN = Tstartbar_F, //bc.st_hot_in.T,   
    WallRes = true, 
    table_th_conductivity = table_k_metalwall,
    // QuasiStatic = metalQuasiStatic , 
    rhomcm = rhomcm 
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
