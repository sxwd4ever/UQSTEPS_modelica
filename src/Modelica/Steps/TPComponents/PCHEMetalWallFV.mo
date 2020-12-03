within Steps.TPComponents;

model PCHEMetalWallFV "simplified PCHE's fin-like metal wall between parallel micro channels"
  extends ThermoPower.Icons.MetalWall;
  import ThermoPower.Choices;
  import ThermoPower.Functions;
  
  parameter Integer Nw = 1 "Number of volumes on the wall ports";
  parameter Integer Nt = 1 "Number of tubes in parallel";
  parameter SI.Length L "Tube length";
  parameter SI.Length r_c "channel raidus";
  parameter SI.Length w_ch = 2.5e-3 "Wall width, mm. see table 2 in Meshram [2016]";
  parameter SI.Length h_ch = 3.2e-3 "height of computational domain, mm. see table 2 in Meshram [2016]";
  parameter SI.Length dx = 0.4e-3 "(equivalent) delta x between channels";
  parameter Real rhomcm "Metal heat capacity per unit volume [J/m^3.K]";
  parameter SI.ThermalConductivity lambda "Thermal conductivity";
  parameter Boolean WallRes = true "Wall thermal resistance accounted for";
  parameter SI.Temperature Tstartbar = 300 "Avarage temperature" annotation(
    Dialog(tab = "Initialisation"));
  parameter SI.Temperature Tstart1 = Tstartbar "Temperature start value - first volume" annotation(
    Dialog(tab = "Initialisation"));
  parameter SI.Temperature TstartN = Tstartbar "Temperature start value - last volume" annotation(
    Dialog(tab = "Initialisation"));
  parameter SI.Temperature Tvolstart[Nw] = Functions.linspaceExt(Tstart1, TstartN, Nw) annotation(
    Dialog(tab = "Initialisation"));
  parameter Choices.Init.Options initOpt = system.initOpt "Initialisation option" annotation(
    Dialog(tab = "Initialisation"));
  constant Real pi = Modelica.Constants.pi;
  // parameter SI.Area Am = (rext ^ 2 - rint ^ 2) * pi "Area of the metal tube cross-section";
  // parameter SI.Area Am =  w_ch * h_ch - pi * r_c * r_c * 4 "Area of the metal tube cross-section";
  parameter SI.Area Am =  w_ch * h_ch - pi * r_c * r_c "Area of the metal tube cross-section";
  parameter Boolean QuasiStatic = true "=True: Dynamic behavior or heat storage will NOT be simulated for metal wall";  
  final parameter SI.HeatCapacity Cm = Nt * L * Am * rhomcm "Total heat capacity";
  outer ThermoPower.System system "System wide properties";
  // Real Q_res[Nw];
  SI.Temperature Tvol[Nw](start = Tvolstart) "Volume temperatures";
  ThermoPower.Thermal.DHTVolumes int(final N = Nw, T(start = Tvolstart)) "Internal surface" annotation(
    Placement(transformation(extent = {{-40, 20}, {40, 40}}, rotation = 0)));
  ThermoPower.Thermal.DHTVolumes ext(final N = Nw, T(start = Tvolstart)) "External surface" annotation(
    Placement(transformation(extent = {{-40, -42}, {40, -20}}, rotation = 0)));

equation
   assert(Am > 0, "Area of the metal wall cross-section must be positive");  
  // L / Nw * Nt * rhomcm * Am * der(Tvol) = int.Q + ext.Q "Energy balance";
  // For now, we consider static-state simulation and ignore the heat storage in metal tube for PCHE;
  // fill(0.0, Nw) = int.Q + ext.Q "Energy balance";
  if QuasiStatic then   
    fill(0.0, Nw) = int.Q + ext.Q "Energy balance";
  else
    L / Nw * Nt * rhomcm * Am * der(Tvol) = int.Q + ext.Q "Energy balance";
  end if;
  if WallRes then
// Thermal resistance of the tube walls accounted for
    int.Q = lambda * pi * r_c * L / Nw * (int.T - Tvol) / dx * Nt "Heat conduction through the internal half-thickness";
    ext.Q = lambda * 2 * r_c * L / Nw * (ext.T - Tvol) / dx * Nt "Heat conduction through the external half-thickness";
  else
// No temperature gradients across the thickness
    ext.T = Tvol;
    int.T = Tvol;
  end if;
initial equation
  if initOpt == Choices.Init.Options.noInit then
// do nothing
  elseif initOpt == Choices.Init.Options.fixedState then
    Tvol = Tvolstart;
  elseif initOpt == Choices.Init.Options.steadyState then
    der(Tvol) = zeros(Nw);
  elseif initOpt == Choices.Init.Options.steadyStateNoT then
// do nothing
  else
    assert(false, "Unsupported initialisation option");
  end if;
  annotation(
    Icon(graphics = {Text(extent = {{-100, 60}, {-40, 20}}, lineColor = {0, 0, 0}, fillColor = {128, 128, 128}, fillPattern = FillPattern.Forward, textString = "Int"), Text(extent = {{-100, -20}, {-40, -60}}, lineColor = {0, 0, 0}, fillColor = {128, 128, 128}, fillPattern = FillPattern.Forward, textString = "Ext"), Text(extent = {{-138, -60}, {142, -100}}, lineColor = {191, 95, 0}, textString = "%name")}),
    Documentation(info = "<HTML>
<p>This is the model of a cylindrical tube of solid material.
<p>The heat capacity (which is lumped at the center of the tube thickness) is accounted for, as well as the thermal resistance due to the finite heat conduction coefficient. Longitudinal heat conduction is neglected.
<p><b>Modelling options</b></p>
<p>The following options are available:
<ul>
<li><tt>WallRes = false</tt>: the thermal resistance of the tube wall is neglected.
<li><tt>WallRes = true</tt>: the thermal resistance of the tube wall is accounted for.
</ul>
</HTML>", revisions = "<html>
<ul>
<li><i>30 May 2005</i>
  by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
     Initialisation support added.</li>
<li><i>1 Oct 2003</i>
  by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
     First release.</li>
</ul>
</html>
    "),
    Diagram(graphics));

end PCHEMetalWallFV;
