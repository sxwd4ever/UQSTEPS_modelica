within Steps.TPComponents;

model SSMixer "Steady state mixer based on ThermoPower's Mixer"
  extends ThermoPower.Icons.Gas.Mixer;
  
  import ThermoPower.Choices;
  import ThermoPower.Gas.FlangeA;
  import ThermoPower.Gas.FlangeB;
  import ThermoPower.Thermal;
  
  replaceable package Medium = Modelica.Media.Interfaces.PartialMedium
    annotation(choicesAllMatching = true);
  Medium.BaseProperties gas(
    p(start=pstart, stateSelect=StateSelect.prefer),
    T(start=Tstart, stateSelect=StateSelect.prefer),
    Xi(start=Xstart[1:Medium.nXi], each stateSelect=StateSelect.prefer));
  parameter SI.Volume V "Inner volume";
  parameter SI.Area S=0 "Inner surface";
  parameter SI.CoefficientOfHeatTransfer gamma=0 "Heat Transfer Coefficient"
    annotation (Evaluate=true);
  parameter SI.HeatCapacity Cm=0 "Metal heat capacity" annotation (Evaluate=true);
  parameter Boolean allowFlowReversal=system.allowFlowReversal
    "= true to allow flow reversal, false restricts to design direction"
    annotation(Evaluate=true);
  outer ThermoPower.System system "System wide properties";
  parameter Medium.AbsolutePressure pstart=1e5 "Pressure start value"
    annotation (Dialog(tab="Initialisation"));
  parameter Medium.Temperature Tstart=300 "Temperature start value"
    annotation (Dialog(tab="Initialisation"));
  parameter Medium.MassFraction Xstart[Medium.nX]=Medium.reference_X
    "Start gas composition" annotation (Dialog(tab="Initialisation"));
  parameter Medium.Temperature Tmstart=Tstart "Metal wall start temperature"
    annotation (Dialog(tab="Initialisation"));
  parameter Choices.Init.Options initOpt=system.initOpt
    "Initialisation option"
    annotation (Dialog(tab="Initialisation"));
  parameter Boolean noInitialPressure=false
    "Remove initial equation on pressure"
    annotation (Dialog(tab="Initialisation"),choices(checkBox=true));
  parameter Boolean noInitialTemperature=false
    "Remove initial equation on temperature"
    annotation (Dialog(tab="Initialisation"),choices(checkBox=true));

  SI.Mass M "Gas total mass";
  SI.InternalEnergy E "Gas total energy";
  SI.Temperature Tm(start=Tmstart) "Wall temperature";
  Medium.SpecificEnthalpy hi1 "Inlet 1 specific enthalpy";
  Medium.SpecificEnthalpy hi2 "Inlet 2 specific enthalpy";
  Medium.SpecificEnthalpy ho "Outlet specific enthalpy";
  Medium.MassFraction Xi1[Medium.nXi] "Inlet 1 composition";
  Medium.MassFraction Xi2[Medium.nXi] "Inlet 2 composition";
  Medium.MassFraction Xo[Medium.nXi] "Outlet composition";
  SI.Time Tr "Residence time";

  FlangeA in1(redeclare package Medium = Medium, m_flow(min=if
          allowFlowReversal then -Modelica.Constants.inf else 0)) annotation (
     Placement(transformation(extent={{-100,40},{-60,80}}, rotation=0)));
  FlangeB out(redeclare package Medium = Medium, m_flow(max=if
          allowFlowReversal then +Modelica.Constants.inf else 0)) annotation (
     Placement(transformation(extent={{80,-20},{120,20}}, rotation=0)));
  FlangeA in2(redeclare package Medium = Medium, m_flow(min=if
          allowFlowReversal then -Modelica.Constants.inf else 0)) annotation (
     Placement(transformation(extent={{-100,-80},{-60,-40}}, rotation=0)));

  replaceable Thermal.HT thermalPort annotation (Placement(transformation(
          extent={{-38,60},{42,80}}, rotation=0)));
equation
  M = gas.d*V "Gas mass";
  E = M*gas.u "Gas internal energy";
  // der(M) = in1.m_flow + in2.m_flow + out.m_flow "Mass balance";
  0 = in1.m_flow + in2.m_flow + out.m_flow "Mass balance"; // steady state mass balance
  /* der(E) = in1.m_flow*hi1 + in2.m_flow*hi2 + out.m_flow*ho - gamma*S*(gas.T
     - Tm) + thermalPort.Q_flow "Energy balance"; */

  0 = in1.m_flow*hi1 + in2.m_flow*hi2 + out.m_flow*ho - gamma*S*(gas.T
     - Tm) + thermalPort.Q_flow "state-state Energy balance";  

  for j in 1:Medium.nXi loop
    0 = in1.m_flow*(Xi1[j] - gas.Xi[j]) + in2.m_flow*(Xi2[j]
       - gas.Xi[j]) + out.m_flow*(Xo[j] - gas.Xi[j])
      "steady state Independent component mass balance";
    /*
    M*der(gas.Xi[j]) = in1.m_flow*(Xi1[j] - gas.Xi[j]) + in2.m_flow*(Xi2[j]
       - gas.Xi[j]) + out.m_flow*(Xo[j] - gas.Xi[j])
      "Independent component mass balance";
    */
  end for;

  if Cm > 0 and gamma > 0 then
    Cm*der(Tm) = gamma*S*(gas.T - Tm) "Metal wall energy balance";
  else
    Tm = gas.T;
  end if;

  // Boundary conditions
  hi1 = homotopy(if not allowFlowReversal then inStream(in1.h_outflow) else
    actualStream(in1.h_outflow), inStream(in1.h_outflow));
  Xi1 = homotopy(if not allowFlowReversal then inStream(in1.Xi_outflow) else
    actualStream(in1.Xi_outflow), inStream(in1.Xi_outflow));
  hi2 = homotopy(if not allowFlowReversal then inStream(in2.h_outflow) else
    actualStream(in2.h_outflow), inStream(in2.h_outflow));
  Xi2 = homotopy(if not allowFlowReversal then inStream(in2.Xi_outflow) else
    actualStream(in2.Xi_outflow), inStream(in2.Xi_outflow));
  ho = homotopy(if not allowFlowReversal then gas.h else actualStream(out.h_outflow),
    gas.h);
  Xo = homotopy(if not allowFlowReversal then gas.Xi else actualStream(out.Xi_outflow),
    gas.Xi);
  in1.p = gas.p;
  in1.h_outflow = gas.h;
  in1.Xi_outflow = gas.Xi;
  in2.p = gas.p;
  in2.h_outflow = gas.h;
  in2.Xi_outflow = gas.Xi;
  out.p = gas.p;
  out.h_outflow = gas.h;
  out.Xi_outflow = gas.Xi;
  thermalPort.T = gas.T;

  Tr = noEvent(M/max(abs(-out.m_flow), Modelica.Constants.eps))
    "Residence time";
initial equation
  // Initial conditions
  if initOpt == Choices.Init.Options.noInit then
    // do nothing
  elseif initOpt == Choices.Init.Options.fixedState then
    if not noInitialPressure then
      gas.p = pstart;
    end if;
    if not noInitialTemperature then
      gas.T = Tstart;
    end if;
    gas.Xi = Xstart[1:Medium.nXi];
    if (Cm > 0 and gamma > 0) then
      Tm  = Tmstart;
    end if;
  elseif initOpt == Choices.Init.Options.steadyState then
    if not noInitialPressure then
      der(gas.p) = 0;
    end if;
    if not noInitialTemperature then
      der(gas.T) = 0;
    end if;
    der(gas.Xi) = zeros(Medium.nXi);
    if (Cm > 0 and gamma > 0) then
      der(Tm) = 0;
    end if;
  elseif initOpt == Choices.Init.Options.steadyStateNoP then
    if not noInitialTemperature then
      der(gas.T) = 0;
    end if;
    der(gas.Xi) = zeros(Medium.nXi);
    if (Cm > 0 and gamma > 0) then
      der(Tm) = 0;
    end if;
  else
    assert(false, "Unsupported initialisation option");
  end if;

  annotation (
    Documentation(info="<html>
<p>This model describes a constant volume mixer with metal walls. The metal wall temperature and the heat transfer coefficient between the wall and the fluid are uniform. The wall is thermally insulated from the outside.</p>
<p><b>Modelling options</b></p>
<p>The actual gas used in the component is determined by the replaceable <tt>Medium</tt> package. In the case of multiple component, variable composition gases, the start composition is given by <tt>Xstart</tt>, whose default value is <tt>Medium.reference_X</tt>.
</html>", revisions="<html>
<ul>
<li><i>30 May 2005</i>
  by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
     Initialisation support added.</li>
<li><i>19 Nov 2004</i>
  by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
     Adapted to Modelica.Media.</li>
<li><i>5 Mar 2004</i>
  by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
     First release.</li>
</ul>
</html>"),
    Icon(graphics),
    Diagram(graphics));

end SSMixer;
