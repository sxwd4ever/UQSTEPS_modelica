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
  Modelica.Blocks.Sources.Step step(height = -1, offset = 1, startTime = 2) annotation(
    Placement(visible = true, transformation(origin = {-88, 108}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Interaction.Show.RealValue realValue annotation(
    Placement(visible = true, transformation(origin = {-20, 108}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Interaction.Show.RealValue realValue1 annotation(
    Placement(visible = true, transformation(origin = {-22, 82}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Math.MultiSwitch multiSwitch1(nu = 2, expr = {2, 2.5}) annotation(
    Placement(visible = true, transformation(origin = {-96, 80}, extent = {{-10, -10}, {30, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.BooleanStep booleanStep(startTime = 220, startValue = true) annotation(
    Placement(visible = true, transformation(origin = {-148, 98}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.BooleanStep booleanStep1(startTime = 480) annotation(
    Placement(visible = true, transformation(origin = {-150, 58}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.Ramp ramp(duration = 1500, height = 86.374, offset = 49.09959, startTime = 900) annotation(
    Placement(visible = true, transformation(origin = {56, 108}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Interaction.Show.RealValue realValue2 annotation(
    Placement(visible = true, transformation(origin = {114, 108}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.TimeTable timeTable(offset = 0.606862, startTime = 0, table = [865, 0; 960, 0.289268; 1200, 0.492262; 1440, 0.576418; 1680, 0.621235; 1920, 0.64151; 2160, 0.65904; 2400, 0.66217]) annotation(
    Placement(visible = true, transformation(origin = {76, 78}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Interaction.Show.RealValue realValue3 annotation(
    Placement(visible = true, transformation(origin = {118, 78}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.TimeTable timeTable1(table = [0, 48.220577; 910, 49.413678; 1015, 63.370455; 1200, 76.784507; 1440, 92.61866; 1680, 105.822528; 1920, 117.414751; 2160, 128.009885; 2340, 135.194766; 2400, 135.473485]) annotation(
    Placement(visible = true, transformation(origin = {76, 44}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Interaction.Show.RealValue realValue4 annotation(
    Placement(visible = true, transformation(origin = {118, 44}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
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
  connect(step.y, realValue.numberPort) annotation(
    Line(points = {{-76, 108}, {-22, 108}, {-22, 108}, {-32, 108}, {-32, 108}}, color = {0, 0, 127}));
  connect(booleanStep.y, multiSwitch1.u[1]) annotation(
    Line(points = {{-137, 98}, {-122.5, 98}, {-122.5, 80}, {-106, 80}}, color = {255, 0, 255}));
  connect(booleanStep1.y, multiSwitch1.u[2]) annotation(
    Line(points = {{-139, 58}, {-106, 58}, {-106, 80}}, color = {255, 0, 255}));
  connect(multiSwitch1.y, realValue1.numberPort) annotation(
    Line(points = {{-64, 80}, {-34, 80}, {-34, 82}, {-34, 82}}, color = {255, 127, 0}));
  connect(ramp.y, realValue2.numberPort) annotation(
    Line(points = {{67, 108}, {102.5, 108}}, color = {0, 0, 127}));
  connect(timeTable.y, realValue3.numberPort) annotation(
    Line(points = {{87, 78}, {106.5, 78}}, color = {0, 0, 127}));
  connect(timeTable1.y, realValue4.numberPort) annotation(
    Line(points = {{88, 44}, {106, 44}, {106, 44}, {106, 44}}, color = {0, 0, 127}));
  annotation(
    experiment(StartTime = 0, StopTime = 2400, Interval = 10, Tolerance = 1e-6),
    __OpenModelica_commandLineOptions = "--matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection -d=initialization,NLSanalyticJacobian,aliasConflicts");
end InputSignals;
