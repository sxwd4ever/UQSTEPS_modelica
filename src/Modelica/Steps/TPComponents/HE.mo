within Steps.TPComponents;

model HE "Heat Exchanger fluid - gas"
  //extends ThermoPower.PowerPlants.HRSG.Interfaces.HeatExchanger;
  extends Interfaces.HeatExchanger;
  
  import Gas = ThermoPower.Gas;  
  import Water = ThermoPower.Water;
  import Choices = ThermoPower.Choices;
  import Thermal = ThermoPower.Thermal;
  import Options = ThermoPower.Choices.Init.Options;
    
  replaceable model HeatTransfer_F = Thermal.HeatTransferFV.IdealHeatTransfer constrainedby ThermoPower.Thermal.BaseClasses.DistributedHeatTransferFV annotation(
     choicesAllMatching = true);
  replaceable model HeatTransfer_G = Thermal.HeatTransferFV.IdealHeatTransfer constrainedby ThermoPower.Thermal.BaseClasses.DistributedHeatTransferFV annotation(
     choicesAllMatching = true);
  replaceable model HeatExchangerTopology = Thermal.HeatExchangerTopologies.CoCurrentFlow constrainedby ThermoPower.Thermal.BaseClasses.HeatExchangerTopologyData annotation(
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
  constant Real pi = Modelica.Constants.pi;
  final parameter SI.Distance L = 1 "Tube length";
  parameter Choices.FluidPhase.FluidPhases FluidPhaseStart = Choices.FluidPhase.FluidPhases.Liquid "Fluid phase (only for initialization!)" annotation(
    Dialog(tab = "Initialization"));
  Water.Flow1DFV fluidFlow(
    Nt = Nt, 
    N = N_F, 
    Nw = Nw_F, 
    wnom = fluidNomFlowRate, 
    initOpt = if SSInit then Options.steadyState else Options.noInit, 
    redeclare package Medium = FluidMedium, 
    L = exchSurface_F ^ 2 / (fluidVol * pi * 4), 
    A = (fluidVol * 4 / exchSurface_F) ^ 2 / 4 * pi,
    omega = fluidVol * 4 / exchSurface_F * pi, 
    Dhyd = fluidVol * 4 / exchSurface_F, 
    FFtype = FFtype_F, 
    HydraulicCapacitance = HCtype_F, 
    Kfnom = Kfnom_F, dpnom = dpnom_F, 
    rhonom = rhonom_F, Cfnom = Cfnom_F, 
    FluidPhaseStart = FluidPhaseStart,
    pstart = pstart_F, 
    redeclare model HeatTransfer = HeatTransfer_F) 
    annotation(
    Placement(transformation(extent = {{-10, -66}, {10, -46}}, rotation = 0)));
  
  // Gas.Flow1DFV gasFlow(  
  Steps.TPComponents.Flow1DFV gasFlow(
    Nt = 1, 
    N = N_G, 
    Nw = Nw_G, 
    Dhyd = 1, 
    wnom = gasNomFlowRate, 
    initOpt = if SSInit then Options.steadyState else Options.noInit, 
    redeclare package Medium = FlueGasMedium, 
    QuasiStatic = gasQuasiStatic, 
    pstart = pstart_G, 
    L = L, 
    A = gasVol / L, 
    omega = exchSurface_G / L, 
    FFtype = FFtype_G, 
    Kfnom = Kfnom_G, 
    dpnom = dpnom_G, 
    rhonom = rhonom_G, 
    Cfnom = Cfnom_G, 
    Tstartbar = Tstartbar_G, 
    Tstartin    = cfg_G.st_in.T,                                              
    Tstartout   = cfg_G.st_out.T,
    hstartin    = cfg_F.st_in.h,
    hstartout   = cfg_F.st_out.h,
    redeclare model HeatTransfer = HeatTransfer_G) 
    annotation(
    Placement(transformation(extent = {{-12, 66}, {12, 46}}, rotation = 0)));
  
  Thermal.MetalTubeFV metalTube(
    Nw = Nw_F, 
    L = exchSurface_F ^ 2 / (fluidVol * pi * 4), 
    rint = fluidVol * 4 / exchSurface_F / 2, 
    rext = (metalVol + fluidVol) * 4 / extSurfaceTub / 2, 
    rhomcm = rhomcm, 
    lambda = lambda,
    Tstartbar = Tstartbar_M,
    WallRes = false) 
    annotation(
    Placement(transformation(extent = {{-10, -24}, {10, -4}})));

  Thermal.HeatExchangerTopologyFV heatExchangerTopology(
    Nw = Nw_F, 
    redeclare model HeatExchangerTopology = HeatExchangerTopology) 
    annotation(
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
  connect(metalTube.ext, fluidFlow.wall) annotation(
    Line(points = {{0, -17.1}, {0, -51}}, color = {255, 127, 0}, smooth = Smooth.None));
  connect(heatExchangerTopology.side2, metalTube.int) annotation(
    Line(points = {{0, 12.9}, {0, -11}}, color = {255, 127, 0}, smooth = Smooth.None));
  connect(heatExchangerTopology.side1, gasFlow.wall) annotation(
    Line(points = {{0, 19}, {0, 51}}, color = {255, 127, 0}, smooth = Smooth.None));
  annotation(
    Diagram(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}}), graphics));
end HE;
