within Steps.TPComponents;

model CompressorFixedPController "A controller to maintain the compressor's outlet pressure by setting the rotational speed"
   
  import ThermoPower.Choices.TurboMachinery.TableTypes;
  replaceable package Medium = Modelica.Media.Interfaces.PartialMedium;
  
  parameter Medium.Temperature Tdes_in "inlet design temperature";
  parameter SI.AngularVelocity Ndesign "Design velocity";
  parameter Real tablePhic[:, :] = fill(0, 0, 2) "Table for phic(N_T,beta)";
  parameter Real tablePR[:, :] = fill(0, 0, 2) "Table for eta(N_T,beta)";
  parameter String fileName = "noName" "File where matrix is stored";
  parameter TableTypes Table = TableTypes.matrix "Selection of the way of definition of table matrix";
  
  Modelica.Blocks.Tables.CombiTable2D PressRatio(tableOnFile = if Table == TableTypes.matrix then false else true, table = tablePR, tableName = if Table == TableTypes.matrix then "NoName" else "tabPR", fileName = if Table == TableTypes.matrix then "NoName" else fileName, smoothness = Modelica.Blocks.Types.Smoothness.ContinuousDerivative) annotation(
    Placement(transformation(extent = {{-12, 0}, {8, 20}}, rotation = 0)));

  Modelica.Blocks.Tables.CombiTable2D Phic(tableOnFile = if Table == TableTypes.matrix then false else true, table = tablePhic, tableName = if Table == TableTypes.matrix then "NoName" else "tabPhic", fileName = if Table == TableTypes.matrix then "NoName" else fileName, smoothness = Modelica.Blocks.Types.Smoothness.ContinuousDerivative) 
annotation(
    Placement(transformation(extent = {{-12, 30}, {8, 50}}, rotation = 0)));

  parameter SI.Pressure p0 = 20e6 "parameter fixed outlet pressure";
  parameter Boolean use_in_p0 = false "Use connector input for the pressure";
  Modelica.Blocks.Interfaces.RealInput in_p0 if use_in_p0;

  Modelica.Blocks.Interfaces.RealInput p_inlet "inlet pressure of the compressor";
  Modelica.Blocks.Interfaces.RealInput T_inlet "inlet temperature of the compressor";
  Modelica.Blocks.Interfaces.RealInput mdot_inlet "inlet mass flow";  
  
  Real N_T "Referred speed ";  
  Real N_T_design "Referred design velocity";
  Real phic "Flow number ";
  Real beta(start = integer(size(tablePhic, 1) / 2)) "Number of beta line";
  SI.PerUnit PR "pressure ratio";
  Modelica.Blocks.Interfaces.RealOutput omega "shaft angular velocity - calculated output";
  
protected
  Modelica.Blocks.Interfaces.RealInput in_p0_internal;  
  
equation
  N_T_design = Ndesign / sqrt(Tdes_in) "Referred design velocity";
  N_T = 100 * omega / (sqrt(T_inlet) * N_T_design) "Referred speed definition, as percentage of design velocity";
  phic = mdot_inlet * sqrt(T_inlet) / p_inlet "Flow number definition";

// for constant pressure output  
  PR = in_p0_internal / p_inlet;

// phic = Phic(beta, N_T)
  Phic.u1 = beta;
  Phic.u2 = N_T;
  phic = Phic.y;

// PR = PressRatio(beta, N_T)  
  PressRatio.u1 = beta;
  PressRatio.u2 = N_T;
  PR = PressRatio.y;  
  
  // set the value of speed component
  // need to connect speed.flange with shaft_a/b outside this component to set compressor's speed
  // can not do it in this component (as in the commented line below)
  // since it will create under-determined situation. 
  
  if not use_in_p0 then
    in_p0_internal = p0 "Pressure set by parameter";
  end if;
    
  // connect(speed.w_ref, omega);  
  // connect(speed.flange, shaft_a);
  
// Connect protected connectors to public conditional connectors
  connect(in_p0, in_p0_internal);
    
  annotation(
    Documentation(info = "<html>
This model adds the performance characteristics to the Compressor_Base model, by means of 2D interpolation tables.</p>
<p>The perfomance characteristics are specified by two characteristic equations: the first relates the flow number <tt>phic</tt>, the pressure ratio <tt>PR</tt> and the referred speed <tt>N_T</tt>; the second relates the efficiency <tt>eta</tt>, the flow number <tt>phic</tt>, and the referred speed <tt>N_T</tt> [1]. To avoid singularities, the two characteristic equations are expressed in parametric form by adding a further variable <tt>beta</tt> (method of beta lines [2]).
<p>The performance maps are thus tabulated into three differents tables, <tt>tablePhic</tt>,  <tt>tablePR</tt> and <tt>tableEta</tt>, which express <tt>phic</tt>, <tt>PR</tt> and <tt>eta</tt> as a function of <tt>N_T</tt> and <tt>beta</tt>, respectively, where <tt>N_T</tt> is the first row while <tt>beta</tt> is the first column. The referred speed <tt>N_T</tt> is defined as a percentage of the design referred speed and <tt>beta</tt> are arbitrary lines, usually drawn parallel to the surge-line on the performance maps.
<p><tt>Modelica.Blocks.Tables.CombiTable2D</tt> interpolates the tables to obtain values of referred flow, pressure ratio and efficiency at given levels of referred speed and beta.
<p><b>Modelling options</b></p>
<p>The following options are available to determine how the table is defined:
<ul><li><tt>Table = 0</tt>: the table is explicitly supplied as matrix parameter.
<li><tt>Table = 1</tt>: the table is read from a file; the string <tt>fileName</tt> contains the path to the files where the tables are stored, either in ASCII or Matlab binary format.
</ul>
<p><b>References:</b></p>
<ol>
<li>S. L. Dixon: <i>Fluid mechanics, thermodynamics of turbomachinery</i>, Oxford, Pergamon press, 1966, pp. 213.
<li>P. P. Walsh, P. Fletcher: <i>Gas Turbine Performance</i>, 2nd ed., Oxford, Blackwell, 2004, pp. 646.
</ol>
</html>", revisions = "<html>
<ul>
<li><i>13 Apr 2005</i>
  by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
     New method for calculating performance parameters using tables.</li>
</li>
<li><i>14 Jan 2005</i>
  by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
     Adapted to Modelica.Media.</li>
<br> Compressor model restructured using inheritance.
</li>
<li><i>5 Mar 2004</i>
  by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
     First release.</li>
</ul>
</html>"));
end CompressorFixedPController;
