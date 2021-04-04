within Steps.TPComponents;

    partial model Flow1DBase "Basic interface for 1-dimensional water/steam gas flow models"
      extends Icons.Gas.Tube;
      import ThermoPower.Choices.Flow1D.FFtypes;
      import ThermoPower.Choices.Flow1D.HCtypes;
      replaceable package Medium = Modelica.Media.Interfaces.PartialMedium annotation(
        choicesAllMatching = true);
      parameter Integer N(min = 2) = 2 "Number of nodes for thermal variables";
      parameter Integer Nw = N - 1 "Number of volumes on the wall interface";
      parameter Integer Nt = 1 "Number of tubes in parallel";
      parameter SI.Distance L "Tube length";
      parameter SI.Position H = 0 "Elevation of outlet over inlet";
      parameter SI.Area A "Cross-sectional area (single tube)";
      parameter SI.Length omega "Perimeter of heat transfer surface (single tube)";
      parameter SI.Length Dhyd "Hydraulic Diameter (single tube)";
      parameter Medium.MassFlowRate wnom "Nominal mass flowrate (total)";
      parameter ThermoPower.Choices.Flow1D.FFtypes FFtype = ThermoPower.Choices.Flow1D.FFtypes.NoFriction "Friction Factor Type" annotation(
        Evaluate = true);
      parameter SI.PressureDifference dpnom = 0 "Nominal pressure drop";
      parameter Real Kfnom = 0 "Nominal hydraulic resistance coefficient" annotation(
        Dialog(enable = FFtype == ThermoPower.Choices.Flow1D.FFtypes.Kfnom));
      parameter Medium.Density rhonom = 0 "Nominal inlet density" annotation(
        Dialog(enable = FFtype == ThermoPower.Choices.Flow1D.FFtypes.OpPoint));
      parameter SI.PerUnit Cfnom = 0 "Nominal Fanning friction factor" annotation(
        Dialog(enable = FFtype == ThermoPower.Choices.Flow1D.FFtypes.Cfnom));
      parameter SI.PerUnit e = 0 "Relative roughness (ratio roughness/diameter)";
      parameter Real Kfc = 1 "Friction factor correction coefficient";
      parameter Boolean DynamicMomentum = false "Inertial phenomena accounted for" annotation(
        Evaluate = true);
      parameter Boolean UniformComposition = true "Uniform gas composition is assumed" annotation(
        Evaluate = true);
      parameter Boolean QuasiStatic = false "Quasi-static model (mass, energy and momentum static balances" annotation(
        Evaluate = true);
      parameter HCtypes HydraulicCapacitance = HCtypes.Downstream "1: Upstream, 2: Downstream";
      parameter Boolean avoidInletEnthalpyDerivative = true "Avoid inlet enthalpy derivative";
      parameter Boolean allowFlowReversal = system.allowFlowReversal "= true to allow flow reversal, false restricts to design direction" annotation(
        Evaluate = true);
      outer ThermoPower.System system "System wide properties";
      parameter Choices.FluidPhase.FluidPhases FluidPhaseStart=Choices.FluidPhase.FluidPhases.Liquid
        "Fluid phase (only for initialization!)"
        annotation (Dialog(tab="Initialisation"));
      parameter Medium.AbsolutePressure pstart = 1e5 "Pressure start value" annotation(
        Dialog(tab = "Initialisation"));
      parameter Medium.SpecificEnthalpy hstartin=if FluidPhaseStart == Choices.FluidPhase.FluidPhases.Liquid
           then 1e5 else if FluidPhaseStart == Choices.FluidPhase.FluidPhases.Steam
           then 3e6 else 1e6 "Inlet enthalpy start value"
        annotation (Dialog(tab="Initialisation"));
      parameter Medium.SpecificEnthalpy hstartout=if FluidPhaseStart == Choices.FluidPhase.FluidPhases.Liquid
           then 1e5 else if FluidPhaseStart == Choices.FluidPhase.FluidPhases.Steam
           then 3e6 else 1e6 "Outlet enthalpy start value"
        annotation (Dialog(tab="Initialisation"));
      parameter Medium.SpecificEnthalpy hstart[N]=linspace(
              hstartin,
              hstartout,
              N) "Start value of enthalpy vector (initialized by default)"
        annotation (Dialog(tab="Initialisation")); 
      parameter Medium.Temperature Tstartbar = 300 "Avarage temperature start value" annotation(
        Dialog(tab = "Initialisation"));
      parameter Medium.Temperature Tstartin = Tstartbar "Inlet temperature start value" annotation(
        Dialog(tab = "Initialisation"));
      parameter Medium.Temperature Tstartout = Tstartbar "Outlet temperature start value" annotation(
        Dialog(tab = "Initialisation"));
      parameter Medium.Temperature Tstart[N] = linspace(Tstartin, Tstartout, N) "Start value of temperature vector (initialized by default)" annotation(
        Dialog(tab = "Initialisation"));
      final parameter SI.Velocity unom = 10 "Nominal velocity for simplified equation";
      parameter Real wnf = 0.01 "Fraction of nominal flow rate at which linear friction equals turbulent friction";
      parameter Medium.MassFraction Xstart[nX] = Medium.reference_X "Start gas composition" annotation(
        Dialog(tab = "Initialisation"));
      parameter Choices.Init.Options initOpt = system.initOpt "Initialisation option" annotation(
        Dialog(tab = "Initialisation"));
      parameter Boolean noInitialPressure = false "Remove initial equation on pressure" annotation(
        Dialog(tab = "Initialisation"),
        choices(checkBox = true));
      function squareReg = ThermoPower.Functions.squareReg;
    protected
      parameter Integer nXi = Medium.nXi "number of independent mass fractions";
      parameter Integer nX = Medium.nX "total number of mass fractions";
    public
      FlangeA infl(redeclare package Medium = Medium, m_flow(start = wnom, min = if allowFlowReversal then -Modelica.Constants.inf else 0)) annotation(
        Placement(transformation(extent = {{-120, -20}, {-80, 20}}, rotation = 0)));
      FlangeB outfl(redeclare package Medium = Medium, m_flow(start = -wnom, max = if allowFlowReversal then +Modelica.Constants.inf else 0)) annotation(
        Placement(transformation(extent = {{80, -20}, {120, 20}}, rotation = 0)));
      SI.Power Q "Total heat flow through the lateral boundary (all Nt tubes)";
      SI.Time Tr "Residence time";
      final parameter SI.PerUnit dzdx=H/L "Slope" annotation (Evaluate=true);
      final parameter SI.Length l=L/(N - 1) "Length of a single volume";
      final parameter SI.Volume V = Nt*A*L "Total volume (all Nt tubes)";    
    
    initial equation
      assert(wnom > 0, "Please set a positive value for wnom");
      assert(FFtype == FFtypes.NoFriction or dpnom > 0, "dpnom=0 not valid, it is also used in the homotopy trasformation during the inizialization");
      assert(not (FFtype == FFtypes.Kfnom and not Kfnom > 0), "Kfnom = 0 not valid, please set a positive value");
      assert(not (FFtype == FFtypes.OpPoint and not rhonom > 0), "rhonom = 0 not valid, please set a positive value");
      assert(not (FFtype == FFtypes.Cfnom and not Cfnom > 0), "Cfnom = 0 not valid, please set a positive value");
      assert(not (FFtype == FFtypes.Colebrook and not Dhyd > 0), "Dhyd = 0 not valid, please set a positive value");
      assert(not (FFtype == FFtypes.Colebrook and not e > 0), "e = 0 not valid, please set a positive value");
      annotation(
        Dialog(enable = FFtype == ThermoPower.Choices.Flow1D.FFtypes.Colebrook),
        Documentation(info = "<HTML>
Basic interface of the <tt>Flow1D</tt> models, containing the common parameters and connectors.
</HTML>
        ", revisions = "<html>
<ul>
<li><i>7 Apr 2014</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       Added base class.</li>

</ul>
</html>"),
        Diagram(graphics),
        Icon(graphics));
    end Flow1DBase;