within Steps.Test;

model TestTP_Mixer
  "StandAlone Component Test For Mixer in ThermoPower"      
   extends Modelica.Icons.Example;
   
  import ThermoPower.Gas;
  import ThermoPower.System; 
  import ThermoPower.Choices; 
//package Medium = Modelica.Media.IdealGases.MixtureGases.CombustionAir;
  package Medium = Media.CO2;// Modelica.Media.IdealGases.MixtureGases.CombustionAir;
  // package Medium = Modelica.Media.IdealGases.SingleGases.CO2;
  
  parameter Real wext=10;
  Components.SSMixer Mixer1(
    redeclare package Medium = Medium,
    gamma=0.8,
S=1,
    V=3,
    pstart=4e6, Tstart=450,
    Tmstart=300) annotation (
      Placement(transformation(extent={{-38,-10},{-18,10}}, rotation=0)));
  //Cm = 0,
  //4e6
  Gas.PressDrop PressDrop1(
    redeclare package Medium = Medium,
    A=0.1,
    dpnom=1e5,
    rhonom=3.5,
    wnom=wext,
    pstart=4e5,
    Tstart=400,
    FFtype=ThermoPower.Choices.PressDrop.FFtypes.OpPoint) annotation (
      Placement(transformation(extent={{0,-10},{22,10}}, rotation=0)));
  Gas.SinkPressure
            SinkP1(
    redeclare package Medium = Medium,
    p0=1e5,
    T=350) annotation (Placement(transformation(extent={{76,-10},{96,10}},
          rotation=0)));
  Gas.SourceMassFlow
              SourceW2(
    redeclare package Medium = Medium,
    w0=15,
    Xnom={1.0}, p0=400000,
    T=350,
    use_in_w0=true) annotation (Placement(transformation(extent={{-76,-40},
            {-56,-20}}, rotation=0)));
  //{0.5,0.5},
  Modelica.Blocks.Sources.Step Step1(
    height=-0.2,
    offset=1.5,
    startTime=15) annotation (Placement(transformation(extent={{20,30},{40,
            50}}, rotation=0)));
  Gas.Valve Valve1(
    redeclare package Medium = Medium,
    wnom=wext,
    CvData=ThermoPower.Choices.Valve.CvTypes.OpPoint,
    pnom=300000,
    dpnom=200000,
    Tstart=400)                                       annotation (Placement(
        transformation(extent={{40,-10},{60,10}}, rotation=0)));
  Modelica.Blocks.Sources.Ramp Ramp1(
    offset=wext,
    height=-1,
    duration=0.1,
    startTime=8) annotation (Placement(transformation(extent={{-100,-20},{-80,
            0}}, rotation=0)));
  Modelica.Blocks.Sources.Ramp Ramp2(
    height=-1,
    offset=5,
    duration=0.1,
    startTime=1) annotation (Placement(transformation(extent={{-100,40},{-80,
            60}}, rotation=0)));
  Gas.SourceMassFlow
              SourceW1(
    redeclare package Medium = Medium,
    p0=400000,
    T=450,
    use_in_w0=true)
           annotation (Placement(transformation(extent={{-74,18},{-54,38}},
          rotation=0)));
  inner System system(initOpt=ThermoPower.Choices.Init.Options.fixedState)
    annotation (Placement(transformation(extent={{80,80},{100,100}})));
equation
  connect(Mixer1.out, PressDrop1.inlet) annotation (Line(
      points={{-18,0},{0,0}},
      color={159,159,223},
      thickness=0.5));
  connect(SourceW2.flange, Mixer1.in2) annotation (Line(
      points={{-56,-30},{-44,-30},{-44,-6},{-36,-6}},
      color={159,159,223},
      thickness=0.5));
  connect(PressDrop1.outlet, Valve1.inlet) annotation (Line(
      points={{22,0},{40,0}},
      color={159,159,223},
      thickness=0.5));
  connect(Valve1.outlet, SinkP1.flange) annotation (Line(
      points={{60,0},{76,0}},
      color={159,159,223},
      thickness=0.5));
  connect(Step1.y, Valve1.theta)
    annotation (Line(points={{41,40},{50,40},{50,7.2}}, color={0,0,127}));
  connect(Ramp1.y, SourceW2.in_w0) annotation (Line(points={{-79,-10},{-72,
          -10},{-72,-25}}, color={0,0,127}));
  connect(SourceW1.flange, Mixer1.in1) annotation (Line(
      points={{-54,28},{-44,28},{-44,6},{-36,6}},
      color={159,159,223},
      thickness=0.5));
  connect(Ramp2.y, SourceW1.in_w0) annotation (Line(points={{-79,50},{-70,
          50},{-70,33}}, color={0,0,127}));
  annotation (
    Diagram(graphics),
    experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-3, Interval = 1),
    __OpenModelica_commandLineOptions = "--matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection -d=initialization,NLSanalyticJacobian,aliasConflicts,bltdump",    
    __OpenModelica_simulationFlags(lv = "LOG_DEBUG,LOG_NLS,LOG_NLS_V,LOG_STATS,LOG_INIT,LOG_STDOUT, -w", outputFormat = "mat", s = "dassl"),
    Documentation(info="<html>
This model tests the <tt>steady state Mixer</tt> model.
</html>"));

end TestTP_Mixer;
