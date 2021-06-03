within Steps.TPComponents;

model RampSeq_TwoStages "Two Stages ramp change sequence which can be used as input of transient simulation"
  extends Modelica.Blocks.Interfaces.SO;
  // x ticks
  parameter Modelica.SIunits.Time time_start = 3 "StartTime of 1st ramp";
  parameter Modelica.SIunits.Time interval   = 1 "Interval between end of 1st and beginning of 2nd ramp";
  parameter Real duration_1                  = 0.15 "Duration for the 1st ramp change";
  parameter Real duration_2                  = 0.15 "Duration for the 2nd ramp change";
  // y ticks
  parameter Real offset   = 1 "Offset of ramps";
  parameter Real height_1 = 0.1 "ramp height for 1st Stage";
  parameter Real height_2 = 0.1 "ramp height for 2nd Stage";  

  Modelica.Blocks.Sources.Ramp ramp1(
    final height    = height_1,
    final duration  = duration_1,
    final startTime = time_start,
    final offset    = 0)
    annotation(
    Placement(visible = true, transformation(extent = {{-82, 66}, {-62, 86}}, rotation = 0)));
  Modelica.Blocks.Sources.Ramp ramp2(
    final height    = height_2,
    final duration  = duration_2,
    final startTime = time_start + interval,
    final offset    = 0)
    annotation(
    Placement(visible = true, transformation(extent = {{-84, 30}, {-64, 50}}, rotation = 0)));
  
  Modelica.Blocks.Math.Add add1 annotation(
    Placement(visible = true, transformation(extent = {{-42, 48}, {-22, 68}}, rotation = 0))); 
  Modelica.Blocks.Math.Add add_all annotation(
    Placement(visible = true, transformation(extent = {{54, -10}, {74, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.Constant const_offset(k = offset) annotation(
    Placement(visible = true, transformation(origin = {16, -50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
equation
  connect(ramp1.y, add1.u1) annotation(
    Line(points = {{-61, 76}, {-52, 76}, {-52, 64}, {-44, 64}}, color = {0, 0, 127}));
  connect(ramp2.y, add1.u2) annotation(
    Line(points = {{-63, 40}, {-52.5, 40}, {-52.5, 52}, {-44, 52}}, color = {0, 0, 127}));
  connect(add1.y, add_all.u1) annotation(
    Line(points = {{-21, 58}, {0, 58}, {0, 30}, {6, 30}}, color = {0, 0, 127}));  
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
end RampSeq_TwoStages;
