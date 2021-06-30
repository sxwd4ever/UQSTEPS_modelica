within Steps.Test.TPComponents;

model TestTP_Compressor_Torque
  "Test for compressor in ThermoPower with varied torque"  
  extends Modelica.Icons.Example;
  package Medium = Modelica.Media.IdealGases.MixtureGases.CombustionAir;

protected
  parameter Real tableEta[6, 4] = [0, 95, 100, 105; 1, 82.5e-2, 81e-2, 80.5e-2; 2, 84e-2, 82.9e-2, 82e-2; 3, 83.2e-2, 82.2e-2, 81.5e-2; 4, 82.5e-2, 81.2e-2, 79e-2; 5, 79.5e-2, 78e-2, 76.5e-2];
  parameter Real tablePhic[6, 4] = [0, 95, 100, 105; 1, 38.3e-3, 43e-3, 46.8e-3; 2, 39.3e-3, 43.8e-3, 47.9e-3; 3, 40.6e-3, 45.2e-3, 48.4e-3; 4, 41.6e-3, 46.1e-3, 48.9e-3; 5, 42.3e-3, 46.6e-3, 49.3e-3];
  parameter Real tablePR[6, 4] = [0, 95, 100, 105; 1, 22.6, 27, 32; 2, 22, 26.6, 30.8; 3, 20.8, 25.5, 29; 4, 19, 24.3, 27.1; 5, 17, 21.5, 24.2];  
public
  // parameter Real N_des = 523.3;
  parameter Real N_des = 2100;
  
  ThermoPower.Gas.SourcePressure SourceP1(
    redeclare package Medium = Medium, 
    p0 = 0.35e5, 
    T = 244.4) 
    annotation(
    Placement(transformation(extent = {{-80, 6}, {-60, 26}}, rotation = 0)));
  ThermoPower.Gas.SinkPressure SinkP1(
    redeclare package Medium = Medium, 
    p0 = 8.3e5, 
    T = 691.4) 
    annotation(
    Placement(transformation(extent = {{40, 6}, {60, 26}}, rotation = 0)));
  Modelica.Mechanics.Rotational.Components.Inertia Inertia1(J = 0.7) annotation(
    Placement(transformation(extent = {{10, -10}, {30, 10}}, rotation = 0)));
  ThermoPower.Gas.Compressor Compressor(
    redeclare package Medium = Medium, 
    pstart_in = 0.35e5, 
    pstart_out = 8.3e5, 
    Tstart_in = 244.4, 
    Tstart_out = 691.4, 
    tablePhic = tablePhic, 
    tableEta = tableEta, 
    tablePR = tablePR, 
    Table = ThermoPower.Choices.TurboMachinery.TableTypes.matrix, 
    explicitIsentropicEnthalpy = false, 
    Ndesign = N_des, 
    Tdes_in = 244.4) 
    annotation(
    Placement(transformation(extent = {{-40, -20}, {0, 20}}, rotation = 0)));
  
  inner ThermoPower.System system annotation(
    Placement(transformation(extent = {{80, 80}, {100, 100}})));
  // hot inlet temperature in K
 
  Modelica.Mechanics.Rotational.Sources.Torque torque_comp(useSupport=false);
 
  Modelica.Blocks.Sources.TimeTable tt_T_motor(
    table = [
      0.28, 84.615; 2.55, 85.714; 20.09, 84.615; 23.21, 124.725;
      26.89, 98.901; 29.15, 100; 90, 94.505; 93.96, 35.165;
      98.49, 72.527; 101.32, 69.231; 159.91, 71.978; 163.3, 109.341;
      167.26, 85.714; 171.23, 86.813; 229.81, 84.615; 233.77, 135.714;
      237.45, 107.692; 241.42, 108.791; 243.96, 100; 246.79, 106.044;
      250.47, 102.747; 299.72, 100.549
    ]) annotation(
    Placement(visible = true, transformation(origin = {76, 44}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
   

initial equation
  Inertia1.w = N_des;
equation
  connect(SourceP1.flange, Compressor.inlet) annotation(
    Line(points = {{-60, 16}, {-36, 16}}, color = {159, 159, 223}, thickness = 0.5));
  connect(Compressor.outlet, SinkP1.flange) annotation(
    Line(points = {{-4, 16}, {40, 16}}, color = {159, 159, 223}, thickness = 0.5));
  connect(Compressor.shaft_b, Inertia1.flange_a) annotation(
    Line(points = {{-8, 0}, {-8, -0.05}, {10, -0.05}, {10, 0}}, color = {0, 0, 0}, thickness = 0.5));
  
  connect(Inertia1.flange_b, torque_comp.flange);

  connect(torque_comp.tau, tt_T_motor.y);
  annotation(
    experiment(StopTime = 100),
    experimentSetupOutput,
    Documentation(info = "<html>
This model test the <tt>Compressor</tt> model with an inertial load. Boundary conditions and data refer to an turbojet engine at 11.000 m.

<p>Simulate for 2 seconds. The compressor slows down.
</html>"));
end TestTP_Compressor_Torque;
