within Steps.Interfaces;

partial model HeatExchangerG2G "Base class for heat exchanger gas - gas (derived from ThermoPower...HE by H Russell 07-08-20"
  import SI = Modelica.SIunits;
  import Gas = ThermoPower.Gas;  
  import Steps.Model.{EntityConfig, EntityGeoParam, EntityThermoParam, HEBoundaryCondition};  
  
  replaceable package FlueGasMedium = ThermoPower.Media.FlueGas constrainedby Modelica.Media.Interfaces.PartialMedium "Flue gas model";
  replaceable package FluidMedium = ThermoPower.Water.StandardWater constrainedby Modelica.Media.Interfaces.PartialPureSubstance "Fluid model";
  
  parameter HEBoundaryCondition bc;
  
  // copy these configs into independant variables - or else Modelica will create template error during translation. 
  parameter EntityGeoParam geo_hot;
  parameter EntityGeoParam geo_cold;
  parameter EntityGeoParam geo_tube;

  parameter EntityThermoParam thermo_hot;
  parameter EntityThermoParam thermo_cold;
  parameter EntityThermoParam thermo_tube; 
  
  parameter Integer N_G = geo_hot.N_seg + 1 "Number of node of the gas side";
  parameter Integer Nw_G = N_G - 1 "Number of volumes of the gas side wall";
  parameter Integer N_F = geo_cold.N_seg + 1 "Number of node of the fluid side";
  parameter Integer Nw_F = N_F - 1 "Number of volumes of the fluid side wall";
  parameter Integer Nt = geo_hot.N_ch "Number of tubes in parallel";
  //Nominal parameter
  parameter SI.MassFlowRate gasNomFlowRate = bc.st_hot_in.mdot"Nominal flow rate through the gas side";
  parameter SI.MassFlowRate fluidNomFlowRate = bc.st_cold_in.mdot "Nominal flow rate through the fluid side";
  parameter SI.Pressure gasNomPressure = bc.st_hot_in.p "Nominal pressure in the gas side inlet";
  parameter SI.Pressure fluidNomPressure = bc.st_cold_in.p"Nominal pressure in the fluid side inlet";
  //Physical Parameter
  parameter SI.Area exchSurface_G = geo_hot.A_ex"Exchange surface between gas - metal tube";
  parameter SI.Area exchSurface_F = geo_cold.A_ex"Exchange surface between metal tube - fluid";
  parameter SI.Area extSurfaceTub = geo_tube.A_ex"Total external surface of the tubes";
  parameter SI.Volume gasVol = geo_hot.V "Gas volume";
  parameter SI.Volume fluidVol = geo_cold.V "Fluid volume";
  parameter SI.Volume metalVol = geo_tube.V "Volume of the metal part in the tubes";
  parameter Real rhomcm = thermo_tube.rho_mcm "Metal heat capacity per unit volume [J/m^3.K]";
  parameter SI.ThermalConductivity lambda = thermo_tube.lambda "Thermal conductivity of the metal (density by specific heat capacity)";
  //Start value
  parameter SI.Temperature Tstartbar_G = bc.st_hot_in.T "Start value of the average gas temperature" annotation(
    Dialog(tab = "Initialization"));
  parameter SI.Temperature Tstartbar_F = bc.st_cold_in.T"Start value of the average fluid temperature" annotation(
    Dialog(tab = "Initialization"));    
  parameter SI.Pressure pstart_G = gasNomPressure "Pressure start value, gas side" annotation(
    Dialog(tab = "Initialization"));
  parameter SI.Temperature Tstartbar_M = Tstartbar_G - 50 "Start value of the average metal temperature" annotation(
    Dialog(tab = "Initialization"));
  parameter SI.Pressure pstart_F = bc.st_cold_in.p "50e5 Pressure start value, fluid side" annotation(
    Dialog(tab = "Initialization"));
  parameter Boolean SSInit = false "Steady-state initialization" annotation(
    Dialog(tab = "Initialization"));
  Gas.FlangeA gasIn(redeclare package Medium = FlueGasMedium) annotation(
    Placement(transformation(extent = {{-110, -10}, {-90, 10}}, rotation = 0)));
  Gas.FlangeB gasOut(redeclare package Medium = FlueGasMedium) annotation(
    Placement(transformation(extent = {{90, -10}, {110, 10}}, rotation = 0)));
  Gas.FlangeA waterIn(redeclare package Medium = FluidMedium) annotation(
    Placement(transformation(extent = {{-10, 90}, {10, 110}}, rotation = 0)));
  Gas.FlangeB waterOut(redeclare package Medium = FluidMedium) annotation(
    Placement(transformation(extent = {{-10, -110}, {10, -90}}, rotation = 0)));
  annotation(
    Diagram(graphics),
    Icon(graphics = {Rectangle(extent = {{-100, 100}, {100, -100}}, lineColor = {0, 0, 255}, fillColor = {230, 230, 230}, fillPattern = FillPattern.Solid), Line(points = {{0, -80}, {0, -40}, {40, -20}, {-40, 20}, {0, 40}, {0, 80}}, color = {0, 0, 255}, thickness = 0.5), Text(extent = {{-100, -115}, {100, -145}}, lineColor = {85, 170, 255}, textString = "%name")}));
end HeatExchangerG2G;
