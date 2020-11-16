within Steps.Components;

model ThroughSink 
  "Sink for component test with fixed inlet temperature and pressure"  
  extends ThermoPower.Icons.Gas.SourceW;  
  
  replaceable package Medium = Modelica.Media.Interfaces.PartialMedium annotation(
    choicesAllMatching = true);
  Medium.BaseProperties gas(p(start = p0), T(start = T), Xi(start = Xnom[1:Medium.nXi]));
  parameter Medium.AbsolutePressure p0 = 101325 "Nominal pressure";
  parameter Medium.Temperature T = 300 "Nominal Temperature";
  parameter Medium.MassFraction Xnom[Medium.nX] = Medium.reference_X "Nominal gas composition";
  parameter Medium.MassFlowRate w0 = 0 "Nominal mass flowrate";
  parameter ThermoPower.Units.HydraulicConductance G = 0 "Hydraulic Conductance";
  parameter Boolean allowFlowReversal = system.allowFlowReversal "= true to allow flow reversal, false restricts to design direction" annotation(
    Evaluate = true);
  parameter Boolean use_in_w0 = false "Use connector input for the nominal flow rate" annotation(
    Dialog(group = "External inputs"),
    choices(checkBox = true));
  parameter Boolean use_in_T = false "Use connector input for the temperature" annotation(
    Dialog(group = "External inputs"),
    choices(checkBox = true));
  parameter Boolean use_in_X = false "Use connector input for the composition" annotation(
    Dialog(group = "External inputs"),
    choices(checkBox = true));
  
  final parameter Real mdot_in;
  final parameter Real T_in;
  
  outer ThermoPower.System system "System wide properties";
  Medium.MassFlowRate w "Nominal mass flow rate";
  ThermoPower.Gas.FlangeA flange(redeclare package Medium = Medium, m_flow(min = if allowFlowReversal then -Modelica.Constants.inf else 0)) annotation(
    Placement(transformation(extent = {{-120, -20}, {-80, 20}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealInput in_w0 if use_in_w0 annotation(
    Placement(transformation(origin = {-60, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 270)));
  Modelica.Blocks.Interfaces.RealInput in_T if use_in_T annotation(
    Placement(transformation(origin = {0, 50}, extent = {{-10, 10}, {10, -10}}, rotation = 270)));
  Modelica.Blocks.Interfaces.RealInput in_X[Medium.nX] if use_in_X annotation(
    Placement(transformation(origin = {60, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 270)));
  Modelica.Blocks.Interfaces.RealOutput out_w annotation(
    Placement(visible = true, transformation(origin = {-30, -82}, extent = {{-10, -10}, {10, 10}}, rotation = -90), iconTransformation(origin = {-30, -82}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
  Modelica.Blocks.Interfaces.RealOutput out_T annotation(
    Placement(visible = true, transformation(origin = {36, -82}, extent = {{-10, -10}, {10, 10}}, rotation = -90), iconTransformation(origin = {36, -82}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
protected
  Modelica.Blocks.Interfaces.RealInput in_w0_internal;
  Modelica.Blocks.Interfaces.RealInput in_T_internal;
  Modelica.Blocks.Interfaces.RealInput in_X_internal[Medium.nX];
equation
  if G > 0 then
    flange.m_flow = w + (flange.p - p0) * G;
  else
    flange.m_flow = w;
  end if;
  w = in_w0_internal;
  if not use_in_w0 then
    in_w0_internal = w0 "Flow rate set by parameter";
  end if;
  gas.T = in_T_internal;
  if not use_in_T then
    in_T_internal = T "Temperature set by parameter";
  end if;
  gas.Xi = in_X_internal[1:Medium.nXi];
  if not use_in_X then
    in_X_internal = Xnom "Composition set by parameter";
  end if;
  flange.p = gas.p;
  flange.h_outflow = gas.h;
  flange.Xi_outflow = gas.Xi;
  
  out_w = mdot_in;
  out_T = T_in;
  /*
  mdot_in + flange.m_flow = 0;
  T_in = gas.T;
  */
// Connect protected connectors to public conditional connectors
  connect(in_w0, in_w0_internal);
  connect(in_T, in_T_internal);
  connect(in_X, in_X_internal);
  annotation(
    Icon(graphics),
    Diagram(graphics),
    Documentation(info = "<html>
<p>The actual gas used in the component is determined by the replaceable <tt>GasModel</tt> model. In the case of multiple component, variable composition gases, the nominal gas composition is given by <tt>Xnom</tt>, whose default value is <tt>Medium.reference_X</tt> .
<p>If <tt>G</tt> is set to zero, the flowrate source is ideal; otherwise, the incoming flowrate increases proportionally to the outlet pressure.</p>
<p>If the <tt>in_w0</tt> connector is wired, then the source massflowrate is given by the corresponding signal, otherwise it is fixed to <tt>w0</tt>.</p>
<p>If the <tt>in_T</tt> connector is wired, then the source temperature is given by the corresponding signal, otherwise it is fixed to <tt>T</tt>.</p>
<p>If the <tt>in_X</tt> connector is wired, then the source massfraction is given by the corresponding signal, otherwise it is fixed to <tt>Xnom</tt>.</p>
</html>", revisions = "<html>
<ul>
<li><i>19 Nov 2004</i>
  by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
     Removed <tt>w0fix</tt> and <tt>Tfix</tt> and <tt>Xfix</tt>; the connection of external signals is now detected automatically.</li> <br> Adapted to Modelica.Media
<li><i>1 Oct 2003</i>
  by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
     First release.</li>
</ul>
</html>"));
end ThroughSink;
