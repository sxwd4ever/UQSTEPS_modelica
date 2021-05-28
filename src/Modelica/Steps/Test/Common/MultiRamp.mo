within Steps.Test.Common;

model MultiRamp "Multiple ramp change as input of transient simulation"
  extends Modelica.Blocks.Interfaces.SO;
  // x ticks
  parameter Modelica.SIunits.Time time_start = 3 "StartTime of 1st ramp";
  parameter Modelica.SIunits.Time interval   = 1 "Interval between end of 1st and beginning of 2nd ramp";
  parameter Real duration                    = 0.15;
  // y ticks
  parameter Real offset     = 1 "Offset of ramps";
  parameter Real ratio_1    = 0.2 "ratio for 1st ramp down and 4th ramp up";
  parameter Real ratio_2    = 0.15 "ratio for 2nd ramp up and 3rd ramp down";
  parameter Real offset_min = offset * (1 - ratio_1);
  parameter Real offset_max = offset * (1 + ratio_2);

  Modelica.Blocks.Sources.Ramp ramp1(
    final height    = offset_min - offset,
    final duration  = duration,
    final startTime = time_start,
    final offset    = 0)
    annotation(
    Placement(visible = true, transformation(extent = {{-82, 66}, {-62, 86}}, rotation = 0)));
  Modelica.Blocks.Sources.Ramp ramp2(
    final height    = offset_max - offset_min,
    final duration  = duration,
    final startTime = time_start + interval,
    final offset    = 0)
    annotation(
    Placement(visible = true, transformation(extent = {{-84, 30}, {-64, 50}}, rotation = 0)));
  Modelica.Blocks.Sources.Ramp ramp3(
    final height    = offset_min - offset_max,
    final duration  = duration,
    final startTime = time_start + interval * 2,
    final offset    = 0)
    annotation(
    Placement(visible = true, transformation(extent = {{-84, -14}, {-64, 6}}, rotation = 0)));
  Modelica.Blocks.Sources.Ramp ramp4(
    final height    = offset - offset_min,
    final duration  = duration,
    final startTime = time_start + interval * 3,
    final offset    = 0)
    annotation(
    Placement(visible = true, transformation(extent = {{-84, -56}, {-64, -36}}, rotation = 0)));
  Modelica.Blocks.Math.Add add1 annotation(
    Placement(visible = true, transformation(extent = {{-42, 48}, {-22, 68}}, rotation = 0)));
  Modelica.Blocks.Math.Add add2 annotation(
    Placement(visible = true, transformation(extent = {{-42, -36}, {-22, -16}}, rotation = 0)));
  Modelica.Blocks.Math.Add add_ramp annotation(
    Placement(visible = true, transformation(extent = {{8, 14}, {28, 34}}, rotation = 0)));
  Modelica.Blocks.Math.Add add_all annotation(
    Placement(visible = true, transformation(extent = {{54, -10}, {74, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.Constant const_offset(k = offset) annotation(
    Placement(visible = true, transformation(origin = {16, -50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
equation
  connect(ramp1.y, add1.u1) annotation(
    Line(points = {{-61, 76}, {-52, 76}, {-52, 64}, {-44, 64}}, color = {0, 0, 127}));
  connect(ramp2.y, add1.u2) annotation(
    Line(points = {{-63, 40}, {-52.5, 40}, {-52.5, 52}, {-44, 52}}, color = {0, 0, 127}));
  connect(add1.y, add_ramp.u1) annotation(
    Line(points = {{-21, 58}, {0, 58}, {0, 30}, {6, 30}}, color = {0, 0, 127}));
  connect(ramp3.y, add2.u1) annotation(
    Line(points = {{-63, -4}, {-48, -4}, {-48, -20}, {-44, -20}}, color = {0, 0, 127}));
  connect(ramp4.y, add2.u2) annotation(
    Line(points = {{-63, -46}, {-48, -46}, {-48, -32}, {-44, -32}}, color = {0, 0, 127}));
  connect(add2.y, add_ramp.u2) annotation(
    Line(points = {{-21, -26}, {0, -26}, {0, 18}, {6, 18}}, color = {0, 0, 127}));
  connect(add_ramp.y, add_all.u1) annotation(
    Line(points = {{29, 24}, {52, 24}, {52, 6}}, color = {0, 0, 127}));
  connect(const_offset.y, add_all.u2) annotation(
    Line(points = {{28, -50}, {52, -50}, {52, -6}, {52, -6}}, color = {0, 0, 127}));
  connect(add_all.y, y) annotation(
    Line(points = {{75, 0}, {110, 0}}, color = {0, 0, 127}));
  annotation(
    Documentation(info = "<HTML>
Block generating the sum of four ramps.
</HTML>"),
    experiment(StartTime = 0, StopTime = 8, Tolerance = 1e-3, Interval = 1),
    Icon(coordinateSystem(initialScale = 0.1), graphics = {Line(points = {{-80, 68}, {-80, -80}}, color = {192, 192, 192}), Polygon(lineColor = {192, 192, 192}, fillColor = {192, 192, 192}, fillPattern = FillPattern.Solid, points = {{-80, 90}, {-88, 68}, {-72, 68}, {-80, 90}}), Line(points = {{-90, -70}, {82, -70}}, color = {192, 192, 192}), Polygon(lineColor = {192, 192, 192}, fillColor = {192, 192, 192}, fillPattern = FillPattern.Solid, points = {{90, -70}, {68, -62}, {68, -78}, {90, -70}}), Line(points = {{-80, -60}, {-50, -60}, {-30, 60}, {10, 60}, {30, -20}, {70, -20}})}));
end MultiRamp;
