within Steps.TPComponents;

model GasFlow1DFV "1-dimensional fluid flow model for gas (finite volumes)"
  extends ThermoPower.Gas.BaseClasses.Flow1DBase;
  import MyUtil = Steps.Utilities.Util;
  import Thermal = ThermoPower.Thermal;
  import Choices = ThermoPower.Choices;
  
  Thermal.DHTVolumes wall(final N = Nw) annotation(
    Placement(transformation(extent = {{-60, 40}, {60, 60}}, rotation = 0)));
  replaceable model HeatTransfer = Thermal.HeatTransferFV.IdealHeatTransfer constrainedby ThermoPower.Thermal.BaseClasses.DistributedHeatTransferFV annotation(
     choicesAllMatching = true);
  HeatTransfer heatTransfer(
    redeclare package Medium = Medium, 
    final Nf = N, 
    final Nw = Nw, 
    final Nt = Nt, 
    final L = L, 
    final A = A, 
    final Dhyd = Dhyd, 
    final omega = omega, 
    final wnom = wnom / Nt, 
    final w = w * ones(N), 
    final fluidState = gas.state) "Instantiated heat transfer model";
    
  parameter SI.PerUnit wnm = 1e-2 "Maximum fraction of the nominal flow rate allowed as reverse flow";
  Medium.BaseProperties gas[N] "Gas nodal properties";
  SI.Pressure Dpfric "Pressure drop due to friction";
  SI.Length omega_hyd "Wet perimeter (single tube)";
  Real Kf "Friction factor";
  Real Kfl "Linear friction factor";
  Real dwdt "Time derivative of mass flow rate";
  SI.PerUnit Cf "Fanning friction factor";
  Medium.MassFlowRate w(start = wnom / Nt) "Mass flowrate (single tube)";
  Medium.Temperature Ttilde[N - 1](start = Tstart[2:N]) "Temperature state variables";
  Medium.Temperature T[N](start = Tstart, each stateSelect = StateSelect.prefer) "Node temperatures";
  Medium.SpecificEnthalpy h[N] "Node specific enthalpies";
  // Medium.MassFraction X[if UniformComposition or Medium.fixedX then 1 else N, nX] "Node mass fraction - all compositions";    
  Medium.Temperature Tin(start = Tstartin);
  Medium.MassFraction Xtilde[if UniformComposition or Medium.fixedX then 1 else N - 1, nX](start = ones(size(Xtilde, 1), size(Xtilde, 2)) * diagonal(Xstart[1:nX]), each stateSelect = StateSelect.prefer) "Composition state variables";
  Medium.MassFlowRate wbar[N - 1](each start = wnom / Nt);
  SI.Power Q_single[N - 1] = heatTransfer.Qvol / Nt "Heat flows entering the volumes from the lateral boundary (single tube)";
  SI.Velocity u[N] "Fluid velocity";
  Medium.AbsolutePressure p(start = pstart, stateSelect = StateSelect.prefer);
  SI.Time Tr "Residence time";
  SI.Mass M "Gas Mass (single tube)";
  SI.Mass Mtot "Gas Mass (total)";
  SI.Power Q "Total heat flow through the wall (all Nt tubes)";
protected
  parameter SI.Length l = L / (N - 1) "Length of a single volume";
  Medium.Density rhobar[N - 1] "Fluid average density";
  SI.SpecificVolume vbar[N - 1] "Fluid average specific volume";
  Medium.DerDensityByPressure drbdp[N - 1] "Derivative of average density by pressure";
  Medium.DerDensityByTemperature drbdT1[N - 1] "Derivative of average density by left temperature";
  Medium.DerDensityByTemperature drbdT2[N - 1] "Derivative of average density by right temperature";
  Real drbdX1[N - 1, nX](each unit = "kg/m3") "Derivative of average density by left composition";
  Real drbdX2[N - 1, nX](each unit = "kg/m3") "Derivative of average density by right composition";
  Medium.SpecificHeatCapacity cvbar[N - 1] "Average cv";
  SI.MassFlowRate dMdt[N - 1] "Derivative of mass in a finite volume";
  Medium.SpecificHeatCapacity cv[N];
  Medium.DerDensityByTemperature dddT[N] "Derivative of density by temperature";
  Medium.DerDensityByPressure dddp[N] "Derivative of density by pressure";
  Real dddX[N, nX](each unit = "kg/m3") "Derivative of density by composition";
equation
  assert(FFtype == ThermoPower.Choices.Flow1D.FFtypes.NoFriction or dpnom > 0, "dpnom=0 not supported, it is also used in the homotopy trasformation during the inizialization");
//All equations are referred to a single tube
// Friction factor selection
  omega_hyd = 4 * A / Dhyd;
  if FFtype == ThermoPower.Choices.Flow1D.FFtypes.Kfnom then
    Kf = Kfnom * Kfc;
    Cf = 2 * Kf * A ^ 3 / (omega_hyd * L);
  elseif FFtype == ThermoPower.Choices.Flow1D.FFtypes.OpPoint then
    Kf = dpnom * rhonom / (wnom / Nt) ^ 2 * Kfc;
    Cf = 2 * Kf * A ^ 3 / (omega_hyd * L);
  elseif FFtype == ThermoPower.Choices.Flow1D.FFtypes.Cfnom then
    Kf = Cfnom * omega_hyd * L / (2 * A ^ 3) * Kfc;
    Cf = Cfnom * Kfc;
  elseif FFtype == ThermoPower.Choices.Flow1D.FFtypes.Colebrook then
    Cf = f_colebrook(w, Dhyd / A, e, Medium.dynamicViscosity(gas[integer(N / 2)].state)) * Kfc;
    Kf = Cf * omega_hyd * L / (2 * A ^ 3);
  elseif FFtype == ThermoPower.Choices.Flow1D.FFtypes.NoFriction then
    Cf = 0;
    Kf = 0;
  else
    assert(false, "Unsupported FFtype");
    Cf = 0;
    Kf = 0;
  end if;
  assert(Kf >= 0, "Negative friction coefficient");
  Kfl = wnom / Nt * wnf * Kf "Linear friction factor";
// Dynamic momentum term
  dwdt = if DynamicMomentum and not QuasiStatic then der(w) else 0;
  sum(dMdt) = (infl.m_flow + outfl.m_flow) / Nt "Mass balance";
  L / A * dwdt + outfl.p - infl.p + Dpfric = 0 "Momentum balance";
  Dpfric = if FFtype == ThermoPower.Choices.Flow1D.FFtypes.NoFriction then 0 else homotopy(smooth(1, Kf * squareReg(w, wnom / Nt * wnf)) * sum(vbar) / (N - 1), dpnom / (wnom / Nt) * w) "Pressure drop due to friction";
  for j in 1:N - 1 loop
    if not QuasiStatic then
// Dynamic mass and energy balances
//      MyUtil.myAssert(debug = false, val_test = Ttilde[j], min = 274, max = 1e6, name_val = "Ttilde[j]", val_ref = {j, A, l, rhobar[j], cvbar[j], wbar[j], gas[j + 1].h,  gas[j].h, Q_single[j]}, name_val_ref = {"j", "A", "l", "rhobar[j]", "cvbar[j]", "wbar[j]", "gas[j + 1].h", "gas[j].h", "Q_single[j]"});
      A * l * rhobar[j] * cvbar[j] * der(Ttilde[j]) + wbar[j] * (gas[j + 1].h - gas[j].h) = Q_single[j] "Energy balance";
      // dMdt[j] = A * l * (drbdp[j] * der(p) + drbdT1[j] * der(T[j]) + drbdT2[j] * der(T[j + 1]) + vector(drbdX1[j, :]) * vector(X[j]) + vector(drbdX2[j, :]) * vector(X[j+1])) "Mass balance";
      dMdt[j] = A * l * (drbdp[j] * der(p) + (drbdT1[j]  + drbdT2[j]) / 2 * der(Ttilde[j]) + vector(dddX[j, :]) * vector(Xtilde[j])) "Mass balance";
       
      
/*
    dMdt[j] = A*l*(drbdT[j]*der(Ttilde[j]) + drbdp[j]*der(p) + vector(drbdX[j, :])*
    vector(der(Xtilde[if UniformComposition then 1 else j, :])))
    "Mass balance";
*/
// Average volume quantities
      if avoidInletEnthalpyDerivative and j == 1 then
// first volume properties computed by the volume outlet properties
        rhobar[j] = gas[j + 1].d;
        drbdp[j] = dddp[j + 1];
        drbdT1[j] = 0;
        drbdT2[j] = dddT[j + 1];
        drbdX1[j, :] = zeros(size(Xtilde, 2));
        drbdX2[j, :] = dddX[j + 1, :];
      else
// volume properties computed by averaging
        rhobar[j] = (gas[j].d + gas[j + 1].d) / 2;
        drbdp[j] = (dddp[j] + dddp[j + 1]) / 2;
        drbdT1[j] = dddT[j] / 2;
        drbdT2[j] = dddT[j + 1] / 2;
        drbdX1[j, :] = dddX[j, :] / 2;
        drbdX2[j, :] = dddX[j + 1, :] / 2;
      end if;
      vbar[j] = 1 / rhobar[j];
      wbar[j] = homotopy(infl.m_flow / Nt - sum(dMdt[1:j - 1]) - dMdt[j] / 2, wnom / Nt);
      cvbar[j] = (cv[j] + cv[j + 1]) / 2;
    else
// Static mass and energy balances
      wbar[j] * (gas[j + 1].h - gas[j].h) = Q_single[j] "Energy balance";
      dMdt[j] = 0 "Mass balance";
// Dummy values for unused average quantities
      rhobar[j] = 0;
      drbdp[j] = 0;
      drbdT1[j] = 0;
      drbdT2[j] = 0;
      drbdX1[j, :] = zeros(nX);
      drbdX2[j, :] = zeros(nX);
      vbar[j] = 0;
      wbar[j] = infl.m_flow / Nt;
      cvbar[j] = 0;
    end if;
  end for;
  Q = heatTransfer.Q "Total heat flow through the lateral boundary";
  if Medium.fixedX then
    Xtilde = fill(Medium.reference_X, 1);
  elseif QuasiStatic then
    Xtilde = fill(gas[1].X, size(Xtilde, 1)) "Gas composition equal to actual inlet";
  elseif UniformComposition then
    der(Xtilde[1, :]) = homotopy(1 / L * sum(u) / N * (gas[1].X - gas[N].X), 1 / L * unom * (gas[1].X - gas[N].X)) "Partial mass balance for the whole pipe";
  else
    for j in 1:N - 1 loop
      der(Xtilde[j, :]) = homotopy((u[j + 1] + u[j]) / (2 * l) * (gas[j].X - gas[j + 1].X), 1 / L * unom * (gas[j].X - gas[j + 1].X)) "Partial mass balance for single volume";
    end for;
  end if;
  for j in 1:N loop
    u[j] = w / (gas[j].d * A) "Gas velocity";
    gas[j].p = p;
    gas[j].T = T[j];
    gas[j].h = h[j];
   
  end for;
// Fluid property computations
  for j in 1:N loop
    if not QuasiStatic then
      cv[j] = Medium.heatCapacity_cv(gas[j].state);
      dddT[j] = Medium.density_derT_p(gas[j].state);
      dddp[j] = Medium.density_derp_T(gas[j].state);
// cv[j] = 0;
// dddT[j] = 0;
// dddp[j] = 0;
      if nX > 0 then
        dddX[j, :] = Medium.density_derX(gas[j].state);
      end if;
    else
// Dummy values (not needed by dynamic equations)
      cv[j] = 0;
      dddT[j] = 0;
      dddp[j] = 0;
      dddX[j, :] = zeros(nX);
    end if;
  end for;
// Selection of representative pressure and flow rate variables
  if HydraulicCapacitance == ThermoPower.Choices.Flow1D.HCtypes.Upstream then
    p = infl.p;
    w = -outfl.m_flow / Nt;
  else
    p = outfl.p;
    w = infl.m_flow / Nt;
  end if;
// Boundary conditions
  infl.h_outflow = gas[1].h;
  outfl.h_outflow = gas[N].h;
  infl.Xi_outflow = gas[1].Xi;
  outfl.Xi_outflow = gas[N].Xi;
  gas[1].h = inStream(infl.h_outflow);
  gas[2:N].T = Ttilde;
  gas[1].Xi = inStream(infl.Xi_outflow);
  
  for j in 2:N loop
    gas[j].Xi = Xtilde[if UniformComposition then 1 else j - 1, 1:nXi];
  end for;  
  
  connect(wall, heatTransfer.wall);
  Tin = gas[1].T;
  M = sum(rhobar) * A * l "Fluid mass (single tube)";
  Mtot = M * Nt "Fluid mass (total)";
  Tr = noEvent(M / max(infl.m_flow / Nt, Modelica.Constants.eps)) "Residence time";
  assert(infl.m_flow > (-wnom * wnm), "Reverse flow not allowed, maybe you connected the component with wrong orientation");
initial equation
  if initOpt == Choices.Init.Options.noInit or QuasiStatic then
// do nothing
  elseif initOpt == Choices.Init.Options.fixedState then
    if not noInitialPressure then
      p = pstart;
    end if;
    Ttilde = Tstart[2:N];
    if not Medium.fixedX then
      Xtilde = ones(size(Xtilde, 1), size(Xtilde, 2)) * diagonal(Xstart[1:nX]);
    end if;
  elseif initOpt == Choices.Init.Options.steadyState then
    if not Medium.singleState and not noInitialPressure then
      der(p) = 0;
    end if;
    der(Ttilde) = zeros(N - 1);    
    if not Medium.fixedX then
      der(Xtilde) = zeros(size(Xtilde, 1), size(Xtilde, 2));
    end if;
  elseif initOpt == Choices.Init.Options.steadyStateNoP then
    der(Ttilde) = zeros(N - 1);
    if not Medium.fixedX then
      der(Xtilde) = zeros(size(Xtilde, 1), size(Xtilde, 2));
    end if;
  else
    assert(false, "Unsupported initialisation option");
  end if;
  annotation(
    Icon(graphics = {Text(extent = {{-100, -60}, {100, -100}}, textString = "%name")}),
    Diagram(graphics),
    Documentation(info = "<html>
<p>This model describes the flow of a gas in a rigid tube. The basic modelling assumptions are:
<ul>
<li>Uniform velocity is assumed on the cross section, leading to a 1-D distributed parameter model.
<li>Turbulent friction is always assumed; a small linear term is added to avoid numerical singularities at zero flowrate. The friction effects are not accurately computed in the laminar and transitional flow regimes, which however should not be an issue in most power generation applications.
<li>The model is based on dynamic mass, momentum, and energy balances. The dynamic momentum term can be switched off, to avoid the fast oscillations that can arise from its coupling with the mass balance (sound wave dynamics).
<li>The longitudinal heat diffusion term is neglected.
<li>The energy balance equation is written by assuming a uniform pressure distribution; the pressure drop is lumped either at the inlet or at the outlet.
<li>The fluid flow can exchange thermal power through the lateral tube boundary, by means of the <tt>wall</tt> connector, that actually represents the wall surface with its temperature. The heat flow is computed by an instance of the replaceable HeatTransfer model; various heat transfer models are available in the ThermoPower.Thermal.HeatTransferFV package.
</ul>
<p>The mass, momentum and energy balance equation are discretised with the finite volume method. The state variables are one pressure, one flowrate (optional), N-1 temperatures, and either one or N-1 gas composition vectors.
<p>The turbulent friction factor can be either assumed as a constant, or computed by Colebrook's equation. In the former case, the friction factor can be supplied directly, or given implicitly by a specified operating point. In any case, the multiplicative correction coefficient <tt>Kfc</tt> can be used to modify the friction coefficient, e.g. to fit experimental data.
<p>A small linear pressure drop is added to avoid numerical singularities at low or zero flowrate. The <tt>wnom</tt> parameter must be always specified: the additional linear pressure drop is such that it is equal to the turbulent pressure drop when the flowrate is equal to <tt>wnf*wnom</tt> (the default value is 1% of the nominal flowrate). Increase <tt>wnf</tt> if numerical problems occur in tubes with very low pressure drops.
<p>Flow reversal is not supported by this model; if you need flow reversal, please consider using the Flow1DFEM model.
<p><b>Modelling options</b></p>
<p>The actual gas used in the component is determined by the replaceable <tt>Medium</tt> package.In the case of multiple component, variable composition gases, the start composition is given by <tt>Xstart</tt>, whose default value is <tt>Medium.reference_X</tt>.
<p>Thermal variables (enthalpy, temperature, density) are computed in <tt>N</tt> equally spaced nodes, including the inlet (node 1) and the outlet (node N); <tt>N</tt> must be greater than or equal to 2.
<p>if <tt>UniformComposition</tt> is true, then a uniform compostion is assumed for the gas through the entire tube length; otherwise, the gas compostion is computed in <tt>N</tt> equally spaced nodes, as in the case of thermal variables.
<p>The following options are available to specify the friction coefficient:
<ul><li><tt>FFtype = FFtypes.Kfnom</tt>: the hydraulic friction coefficient <tt>Kf</tt> is set directly to <tt>Kfnom</tt>.
<li><tt>FFtype = FFtypes.OpPoint</tt>: the hydraulic friction coefficient is specified by a nominal operating point (<tt>wnom</tt>,<tt>dpnom</tt>, <tt>rhonom</tt>).
<li><tt>FFtype = FFtypes.Cfnom</tt>: the friction coefficient is computed by giving the (constant) value of the Fanning friction factor <tt>Cfnom</tt>.
<li><tt>FFtype = FFtypes.Colebrook</tt>: the Fanning friction factor is computed by Colebrook's equation (assuming Re > 2100, e.g. turbulent flow).
<li><tt>FFtype = FFtypes.NoFriction</tt>: no friction is assumed across the pipe.</ul>
<p>If <tt>QuasiStatic</tt> is set to true, the dynamic terms are neglected in the mass, momentum, and energy balances, i.e., quasi-static behaviour is modelled. It is also possible to neglect only the dynamic momentum term by setting <tt>DynamicMomentum = false</tt>.
<p>If <tt>HydraulicCapacitance = 2</tt> (default option) then the mass buildup term depending on the pressure is lumped at the outlet, while the optional momentum buildup term depending on the flowrate is lumped at the inlet; therefore, the state variables are the outlet pressure and the inlet flowrate. If <tt>HydraulicCapacitance = 1</tt> the reverse takes place.
<p>Start values for the pressure and flowrate state variables are specified by <tt>pstart</tt>, <tt>wstart</tt>. The start values for the node temperatures are linearly distributed from <tt>Tstartin</tt> at the inlet to <tt>Tstartout</tt> at the outlet. The (uniform) start value of the gas composition is specified by <tt>Xstart</tt>.
<p>A bank of <tt>Nt</tt> identical tubes working in parallel can be modelled by setting <tt>Nt > 1</tt>. The geometric parameters always refer to a <i>single</i> tube.
<p>This models makes the temperature and external heat flow distributions available to connected components through the <tt>wall</tt> connector. If other variables (e.g. the heat transfer coefficient) are needed by external components to compute the actual heat flow, the <tt>wall</tt> connector can be replaced by an extended version of the <tt>DHT</tt> connector.
</html>", revisions = "<html>
<ul>
<li><i>30 May 2005</i>
  by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
     Initialisation support added.</li>
<li><i>24 Mar 2005</i>
  by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
     <tt>QuasiStatic</tt> added.<br>
     <tt>FFtypes</tt> package and <tt>NoFriction</tt> option added.</li>
<li><i>19 Nov 2004</i>
  by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
     Adapted to Modelica.Media.</li>
<li><i>5 Mar 2004</i>
  by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
     First release.</li>
</ul>
</html>"));
end GasFlow1DFV;
