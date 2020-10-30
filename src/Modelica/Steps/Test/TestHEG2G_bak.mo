within Steps.Test;

model TestHEG2G_bak
  "backup of Test for HEG2G"     
      package FlueGasMedium = Steps.Media.SCO2;
      package FluidMedium = Steps.Media.SCO2;
      //gas
      parameter Modelica.SIunits.MassFlowRate gasNomFlowRate = 125 "Nominal mass flowrate";
      parameter Modelica.SIunits.Pressure gasNomPressure = 9e6 "Nominal pressure in the gas side inlet";
      parameter Modelica.SIunits.Temperature Tstart_G_In = 883 "Inlet gas temperature start value";
      parameter Modelica.SIunits.Temperature Tstart_G_Out = 643 "Outlet gas temperature start value";
      //fluid
      parameter Modelica.SIunits.MassFlowRate fluidNomFlowRate = 125 "Nominal flow rate through the fluid side";
      parameter Modelica.SIunits.Pressure fluidNomPressure = 9e6 "Nominal pressure in the fluid side inlet";
      parameter Modelica.SIunits.CoefficientOfHeatTransfer gamma_G = 200 "Constant heat transfer coefficient in the gas side";
      parameter Modelica.SIunits.CoefficientOfHeatTransfer gamma_F = 200 "Constant heat transfer coefficient in the fluid side";
      parameter Modelica.SIunits.Temperature Tstart_M_In = Tstart_F_In "Inlet metal wall temperature start value";
      parameter Modelica.SIunits.Temperature Tstart_M_Out = Tstart_F_Out "Outlet metal wall temperature start value";
      parameter Modelica.SIunits.Temperature Tstart_F_In = 633 "Inlet fluid temperature start value";
      parameter Modelica.SIunits.Temperature Tstart_F_Out = 843 "Outlet fluid temperature start value";
      //parameter Modelica.SIunits.SpecificEnthalpy hstart_F_In = FluidMedium.specificEnthalpy_pT(fluidNomPressure, Tstart_F_In) "Nominal specific enthalpy";
      //parameter Modelica.SIunits.SpecificEnthalpy hstart_F_Out = FluidMedium.specificEnthalpy_pT(fluidNomPressure, Tstart_F_Out) "Nominal specific enthalpy";
      //Components
      inner ThermoPower.System system(allowFlowReversal = false) annotation(
        Placement(transformation(extent = {{80, 80}, {100, 100}})));
      ThermoPower.Gas.SourceMassFlow sourceW_water(redeclare package Medium = FluidMedium, T = Tstart_F_In, p0 = fluidNomPressure, use_in_T = false, w0 = fluidNomFlowRate) annotation(
        Placement(transformation(origin = {0, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 270)));
      ThermoPower.Gas.SinkPressure sinkP_water(redeclare package Medium = FluidMedium, p0 = fluidNomPressure, T = Tstart_F_Out) annotation(
        Placement(transformation(origin = {0, -80}, extent = {{-10, -10}, {10, 10}}, rotation = 270)));
      ThermoPower.Gas.SinkPressure sinkP_gas(redeclare package Medium = FlueGasMedium, T = Tstart_G_Out, p0 = gasNomPressure) annotation(
        Placement(transformation(extent = {{60, -10}, {80, 10}}, rotation = 0)));
      ThermoPower.Gas.SourceMassFlow sourceW_gas(redeclare package Medium = FlueGasMedium, T = Tstart_G_In, p0 = gasNomPressure, w0 = gasNomFlowRate) annotation(
        Placement(transformation(extent = {{-70, -10}, {-50, 10}}, rotation = 0)));
      ThermoPower.Gas.SensT T_waterOut(redeclare package Medium = FluidMedium) annotation(
        Placement(transformation(origin = {4, -50}, extent = {{-10, -10}, {10, 10}}, rotation = 270)));
      ThermoPower.Gas.SensT T_gasOut(redeclare package Medium = FlueGasMedium) annotation(
        Placement(transformation(extent = {{30, -6}, {50, 14}}, rotation = 0)));
      ThermoPower.PowerPlants.HRSG.Components.HEG2G hE(
      redeclare package FluidMedium = FluidMedium, 
      redeclare package FlueGasMedium = FlueGasMedium, 
      redeclare replaceable model HeatTransfer_F = ThermoPower.Thermal.HeatTransferFV.ConstantHeatTransferCoefficient(gamma = gamma_F), 
      redeclare replaceable model HeatTransfer_G = ThermoPower.Thermal.HeatTransferFV.ConstantHeatTransferCoefficientTwoGrids(gamma = gamma_G), 
      redeclare model HeatExchangerTopology = ThermoPower.Thermal.HeatExchangerTopologies.CounterCurrentFlow,  
      N_F = 7, 
      N_G = 7, 
      Nw_G = 6, 
      SSInit = true, 
      Tstartbar_G = 823.15, 
      Tstartbar_M = 773.15, 
      exchSurface_F = 225.073, 
      exchSurface_G = 1708.2, 
      extSurfaceTub = 252.286, 
      fluidNomFlowRate = fluidNomFlowRate, 
      fluidNomPressure = fluidNomPressure, 
      fluidVol = 2.234, 
      gasNomFlowRate = gasNomFlowRate, 
      gasNomPressure = gasNomPressure, 
      gasVol = 10, 
      lambda = 20, 
      metalVol = 0.573, 
      pstart_F = fluidNomPressure, 
      rhomcm = 7900 * 578.05) annotation(
        Placement(transformation(extent = {{-20, -20}, {20, 20}}, rotation = 0)));
      //Start value
      parameter Modelica.SIunits.Temperature Tstart_G = (Tstart_G_In + Tstart_G_Out) / 2;
      parameter Modelica.SIunits.Temperature Tstart_M = (Tstart_G_In + Tstart_G_Out + Tstart_F_In + Tstart_F_Out) / 4;
      parameter Boolean SSInit = true "Steady-state initialization";
      
    initial equation
//hstart_F_Out = hE.waterOut.h_outflow;
    equation
      connect(T_gasOut.inlet, hE.gasOut) annotation(
        Line(points = {{34, 0}, {34, 0}, {20, 0}}, color = {159, 159, 223}, thickness = 0.5));
      connect(T_gasOut.outlet, sinkP_gas.flange) annotation(
        Line(points = {{46, 0}, {46, 0}, {60, 0}}, color = {159, 159, 223}, thickness = 0.5));
      connect(sinkP_water.flange, T_waterOut.outlet) annotation(
        Line(points = {{1.83697e-015, -70}, {1.83697e-015, -56}, {-8.88178e-016, -56}}, thickness = 0.5, color = {0, 0, 255}));
      connect(T_waterOut.inlet, hE.waterOut) annotation(
        Line(points = {{8.88178e-016, -44}, {8.88178e-016, -20}, {0, -20}}, thickness = 0.5, color = {0, 0, 255}));
      connect(sourceW_water.flange, hE.waterIn) annotation(
        Line(points = {{-1.83697e-015, 50}, {-1.83697e-015, 20}, {0, 20}}, color = {0, 0, 255}, thickness = 0.5, smooth = Smooth.None));
      connect(sourceW_gas.flange, hE.gasIn) annotation(
        Line(points = {{-50, 0}, {-20, 0}}, color = {159, 159, 223}, thickness = 0.5, smooth = Smooth.None));

annotation(
    Diagram(graphics),
    experiment(StartTime = 0, StopTime = 10, Tolerance = 1e-3, Interval = 0.5),
    __OpenModelica_commandLineOptions = "--matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection -d=initialization,NLSanalyticJacobian,aliasConflicts,bltdump");  
end TestHEG2G_bak;
