within Steps.TPComponents;

model CompressorPham "Gas compressor using alternative map methodology based on equivalent conditions."
  extends BaseClasses.CompressorBase;
  import ThermoPower.Choices.TurboMachinery.TableTypes;
  parameter SI.AngularVelocity Ndesign "Design velocity";
  parameter SI.Pressure Pdes_in "Design inlet pressure";
  parameter Real ns_des;
  parameter Real Zdes_in;
  parameter Real R_CO2;
  
  parameter Real table_mdote[:, :]=fill(
        0,
        0,
        2) "Table for equivalent mass flow mdote (N_T,beta)";
  parameter Real table_eta[:, :]=fill(
        0,
        0,
        2) "Table for eta(N_T,beta)";
  parameter Real table_dhe[:, :]=fill(
        0,
        0,
        2) "Table for equivelent enthalpy rise dhe (N_T,beta)";
  parameter String fileName="noName" "File where matrix is stored";
  parameter TableTypes Table = TableTypes.matrix
    "Selection of the way of definition of table matrix";
  Modelica.Blocks.Tables.CombiTable2D Eta(
    tableOnFile = if (Table == TableTypes.matrix) then false else true,
    table       = table_eta,
    tableName   = if (Table == TableTypes.matrix) then "NoName" else "tabEta",
    fileName    = if (Table == TableTypes.matrix) then "NoName" else fileName,
    smoothness  = Modelica.Blocks.Types.Smoothness.ContinuousDerivative)
    annotation (Placement(transformation(extent={{-12,60},{8,80}}, rotation=0)));
  Modelica.Blocks.Tables.CombiTable2D Dhe(
    tableOnFile = if (Table == TableTypes.matrix) then false else true,
    table       = table_dhe,
    tableName   = if (Table == TableTypes.matrix) then "NoName" else "tabPR",
    fileName    = if (Table == TableTypes.matrix) then "NoName" else fileName,
    smoothness  = Modelica.Blocks.Types.Smoothness.ContinuousDerivative)
    annotation (Placement(transformation(extent={{-12,0},{8,20}}, rotation=0)));
  Modelica.Blocks.Tables.CombiTable2D Mdote(
    tableOnFile = if (Table == TableTypes.matrix) then false else true,
    table       = table_mdote,
    tableName   = if (Table == TableTypes.matrix) then "NoName" else "tabMdote",
    fileName    = if (Table == TableTypes.matrix) then "NoName" else fileName,
    smoothness  = Modelica.Blocks.Types.Smoothness.ContinuousDerivative)
    annotation (Placement(transformation(extent={{-12,30},{8,50}}, rotation=0)));
  Real N_T "Referred speed ";
  Real mdote "Flow number ";  
  Real beta(start=integer(size(table_mdote, 1)/2)) "Number of beta line";  
  Real Z_in;
  Real ns_in;
  Real dh;

  //Medium.SpecificEnthalpy dh "Enthalpy rise";
  Medium.ThermodynamicState state_iso;

equation
   
  // Z 
  Z_in = gas_in.p / (gas_in.d*R_CO2*gas_in.T);

  //ns
  ns_in = (gas_in.state.cp/gas_in.state.cv)*(1/gas_in.p)*(1/gas_in.state.kappa);

  // Equivalent conditions (mdote, N_T = Neq/Ndes)
  mdote = w * ( sqrt(ns_in*Z_in*gas_in.T)/(ns_in*gas_in.p) ) * ( (ns_des*Pdes_in)/sqrt(ns_des*Zdes_in*Tdes_in) ) "Equivalent mass flow definition";
  N_T = 100* (1/Ndesign) * omega * ( 1/sqrt(ns_in*Z_in*gas_in.T) ) * sqrt(ns_des*Zdes_in*Tdes_in) "Referred speed definition, as percentage of design velocity";

  // Mdote = Mdote(beta, N_T)
  Mdote.u1 = beta;
  Mdote.u2 = N_T;
  mdote    = Mdote.y;
  
  // eta = Eta(beta, N_T)
  Eta.u1 = beta;
  Eta.u2 = N_T;
  eta    = Eta.y;
  
  // PR = Dhe(beta, N_T)
  Dhe.u1 = beta;
  Dhe.u2 = N_T;  

  dh = Dhe.y * ( ns_in*Z_in*gas_in.T ) * ( 1/(ns_des*Zdes_in*Tdes_in));

  state_iso = Medium.setState_hs(gas_in.h+dh*eta, gas_in.state.s);
  
  PR = state_iso.p/gas_in.p;
  
  annotation (Documentation(info="<html>
This model adds the performance characteristics to the Compressor_Base model, by means of 2D interpolation tables. This model is derived from the ThermoPower.Gas.Compressor class. It implements the methodology of </p>
<p>The perfomance characteristics are specified by two characteristic equations: the first relates the equivalent mass flow rate <tt>mdote</tt>, the enthalpy rise <tt>dh</tt> and the referred speed <tt>N_T</tt>; the second relates the efficiency <tt>eta</tt>, the equivalent mass flow <tt>mdote</tt>, and the referred speed <tt>N_T</tt> [1]. To avoid singularities, the two characteristic equations are expressed in parametric form by adding a further variable <tt>beta</tt> (method of beta lines [2]).
<p>The performance maps are thus tabulated into three differents tables, <tt>table_mdote</tt>,  <tt>table_dhe</tt> and <tt>table_eta</tt>, which express <tt>mdote</tt>, <tt>dh</tt> and <tt>eta</tt> as a function of <tt>N_T</tt> and <tt>beta</tt>, respectively, where <tt>N_T</tt> is the first row while <tt>beta</tt> is the first column. The referred speed <tt>N_T</tt> is defined as a percentage of the design referred speed and <tt>beta</tt> are arbitrary lines, usually drawn parallel to the surge-line on the performance maps.
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
<li>H. S. Pham et al: <i>An approach for establishing the performance maps of the sc-CO2 compressor: Development and qualification by means of CFD simulations</i>, International Journal of Heat and Fluid Flow, vol. 61, pp. 379–394, Oct. 2016, doi: 10.1016/j.ijheatfluidflow.2016.05.017.

</ol>
</html>", revisions="<html>
<ul>
<li><i>19 May 2021</i>
  by <a href=\"mailto:h.russell@uq.edu.au\">Hugh Russell</a>:<br>
     Pham et al (2016) equivalent condition methodology for sCO2 compressor off-design performance.</li>
</li>
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
end CompressorPham;