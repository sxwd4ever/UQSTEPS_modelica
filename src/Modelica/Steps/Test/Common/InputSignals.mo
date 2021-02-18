within Steps.Test.Common;

model InputSignals
  Modelica.Blocks.Sources.IntegerConstant integerConstant(k = 20) annotation(
    Placement(visible = true, transformation(origin = {-84, 52}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.MathInteger.TriggeredAdd triggeredAdd annotation(
    Placement(visible = true, transformation(origin = {-18, 52}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
  Modelica.Blocks.Sources.BooleanPulse booleanPulse(period = 2, width = 10) annotation(
    Placement(visible = true, transformation(origin = {-56, 18}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Interaction.Show.IntegerValue integerValue annotation(
    Placement(visible = true, transformation(origin = {26, 78}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.MathInteger.Sum sum1(nu = 2) annotation(
    Placement(visible = true, transformation(origin = {26, 52}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.IntegerConstant integerConstant1(k = 630) annotation(
    Placement(visible = true, transformation(origin = {0, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Interaction.Show.IntegerValue integerValue1 annotation(
    Placement(visible = true, transformation(origin = {64, 22}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.MathInteger.Sum sum(nu = 2) annotation(
    Placement(visible = true, transformation(origin = {28, -38}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.IntegerConstant integerConstant2(k = 400) annotation(
    Placement(visible = true, transformation(origin = {0, -88}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Interaction.Show.IntegerValue integerValue3 annotation(
    Placement(visible = true, transformation(origin = {64, 4}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.IntegerConstant integerConstant3(k = 20) annotation(
    Placement(visible = true, transformation(origin = {-86, -40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.MathInteger.TriggeredAdd triggeredAdd1 annotation(
    Placement(visible = true, transformation(origin = {-28, -42}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
  Modelica.Blocks.Sources.BooleanPulse booleanPulse1(period = 2, startTime = 1, width = 10) annotation(
    Placement(visible = true, transformation(origin = {-58, -76}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
equation
  connect(integerConstant.y, triggeredAdd.u) annotation(
    Line(points = {{-73, 52}, {-26, 52}}, color = {255, 127, 0}));
  connect(triggeredAdd.y, integerValue.numberPort) annotation(
    Line(points = {{-11, 52}, {2, 52}, {2, 78}, {14, 78}}, color = {255, 127, 0}));
  connect(booleanPulse.y, triggeredAdd.trigger) annotation(
    Line(points = {{-45, 18}, {-22, 18}, {-22, 45}}, color = {255, 0, 255}));
  connect(triggeredAdd.y, sum1.u[1]) annotation(
    Line(points = {{-11, 52}, {16, 52}}, color = {255, 127, 0}));
  connect(integerConstant1.y, sum1.u[2]) annotation(
    Line(points = {{11, 0}, {16, 0}, {16, 52}}, color = {255, 127, 0}));
  connect(sum1.y, integerValue1.numberPort) annotation(
    Line(points = {{38, 52}, {45, 52}, {45, 22}, {52, 22}}, color = {255, 127, 0}));
  connect(integerConstant2.y, sum.u[1]) annotation(
    Line(points = {{11, -88}, {11, -64}, {18, -64}, {18, -38}}, color = {255, 127, 0}));
  connect(sum.y, integerValue3.numberPort) annotation(
    Line(points = {{40, -38}, {52, -38}, {52, 4}}, color = {255, 127, 0}));
  connect(integerConstant3.y, triggeredAdd1.u) annotation(
    Line(points = {{-74, -40}, {-38, -40}, {-38, -42}, {-36, -42}}, color = {255, 127, 0}));
  connect(booleanPulse1.y, triggeredAdd1.trigger) annotation(
    Line(points = {{-46, -76}, {-32, -76}, {-32, -50}, {-32, -50}}, color = {255, 0, 255}));
  connect(triggeredAdd1.y, sum.u[2]) annotation(
    Line(points = {{-20, -42}, {18, -42}, {18, -38}, {18, -38}}, color = {255, 127, 0}));
  annotation(
    experiment(StartTime = 0, StopTime = 10, Interval = 0.1, Tolerance = 1e-6),
    __OpenModelica_commandLineOptions = "--matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection -d=initialization,NLSanalyticJacobian,aliasConflicts");
end InputSignals;
