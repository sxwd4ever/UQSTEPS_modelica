within Steps.TPComponents;

model Flow1DFV 
  "1-dimensional fluid flow model for water/steam (finite volumes)"
    extends Flow1DBase;
    import ThermoPower.Choices.Flow1D.FFtypes;
    import ThermoPower.Choices.Flow1D.HCtypes;

    ThermoPower.Thermal.DHTVolumes wall(final N=Nw)
      annotation (Placement(transformation(extent={{-40,40},{40,60}},
            rotation=0)));
    
    replaceable model HeatTransfer = ThermoPower.Thermal.HeatTransferFV.IdealHeatTransfer
      constrainedby ThermoPower.Thermal.BaseClasses.DistributedHeatTransferFV
      annotation (choicesAllMatching=true);
    HeatTransfer heatTransfer(
      redeclare package Medium = Medium,
      final Nf=N,
      final Nw=Nw,
      final Nt=Nt,
      final L=L,
      final A=A,
      final Dhyd=Dhyd,
      final omega=omega,
      final wnom=wnom/Nt,
      final w=w*ones(N),
      final fluidState=gas) "Instantiated heat transfer model";
    
    parameter SI.PerUnit wnm = 1e-3 "Maximum fraction of the nominal flow rate allowed as reverse flow";
    parameter Boolean fixedMassFlowSimplified = false "Fix flow rate = wnom for simplified homotopy model"
        annotation (Dialog(tab="Initialisation"));

    Medium.ThermodynamicState gas[N]
      "Thermodynamic state of the fluid at the nodes";
    SI.Pressure Dpfric "Pressure drop due to friction (total)";
    SI.Length omega_hyd "Wet perimeter (single tube)";
    Medium.MassFlowRate win "Flow rate at the inlet (single tube)";
    Medium.MassFlowRate wout "Flow rate at the outlet (single tube)";    
    Real Kf "Hydraulic friction coefficient";
    Real Kfl "Linear friction factor";
    Real dwdt "Dynamic momentum term";
    SI.PerUnit Cf "Fanning friction factor";
    Medium.MassFlowRate w(start=wnom/Nt) "Mass flowrate (single tube)";
    // Medium.SpecificEnthalpy htilde[N - 1](start=hstart[2:N],each stateSelect=StateSelect.prefer)
    //   "Enthalpy state variables";
    Medium.Temperature Ttilde[N - 1](start=ones(N - 1)*Tstartin + (1:(N - 1))/
        (N - 1)*(Tstartout - Tstartin), each stateSelect=StateSelect.prefer)
    "Temperature state variables";      
    Medium.Temperature T[N] "Node temperatures";
    Medium.SpecificEnthalpy h[N] "Node specific enthalpies";
    // Medium.Temperature Tin(start=Tstartin);    
    Medium.MassFlowRate wbar[N - 1](each start=wnom/Nt)
      "Average flow rate through volumes (single tube)";
    SI.Power Q_single[N-1] = heatTransfer.Qvol/Nt
      "Heat flows entering the volumes from the lateral boundary (single tube)";
  //   MassFlowRate wstar[N];
    SI.Velocity u[N] "Fluid velocity";
    Medium.AbsolutePressure p(start=pstart,stateSelect=StateSelect.prefer) "Fluid pressure for property calculations";
    // Tr was defined in BaseClass Flow1DBase
    Medium.Density rho[N] "Fluid nodal density";
    SI.Mass M "Fluid mass (single tube)";
    SI.Mass Mtot "Fluid mass (total)";
    // SI.Power Q "Total heat flow through the wall (all Nt tubes)" - Defined in baseClass;    
/*
  variables in Water.Flow1DFV
    SI.Pressure Dpfric1
      "Pressure drop due to friction (from inlet to capacitance)";
    SI.Pressure Dpfric2
      "Pressure drop due to friction (from capacitance to outlet)";
    SI.Pressure Dpstat "Pressure drop due to static head";
    Medium.MassFlowRate win "Flow rate at the inlet (single tube)";
    Medium.MassFlowRate wout "Flow rate at the outlet (single tube)";  
*/    
  
 protected
    Medium.Density rhobar[N - 1] "Fluid average density";
    SI.SpecificVolume vbar[N - 1] "Fluid average specific volume";
    //HeatFlux phibar[N - 1] "Average heat flux";
    SI.DerDensityByEnthalpy drdh[N] "Derivative of density by enthalpy";
    SI.DerDensityByEnthalpy drbdh[N - 1]
      "Derivative of average density by enthalpy";
    SI.DerDensityByPressure drdp[N] "Derivative of density by pressure";
    SI.DerDensityByPressure drbdp[N - 1]
      "Derivative of average density by pressure";     

    Medium.DerDensityByTemperature drbdT1[N - 1]
      "Derivative of average density by left temperature";
    Medium.DerDensityByTemperature drbdT2[N - 1]
      "Derivative of average density by right temperature";
    /*
    Real drbdX1[N - 1, nX](each unit="kg/m3")
      "Derivative of average density by left composition";
    Real drbdX2[N - 1, nX](each unit="kg/m3")
      "Derivative of average density by right composition";
    */
    Medium.SpecificHeatCapacity cvbar[N - 1] "Average cv";
    SI.MassFlowRate dMdt[N - 1] "Derivative of mass in a finite volume";
    Medium.SpecificHeatCapacity cv[N];
    Medium.DerDensityByTemperature dddT[N]
      "Derivative of density by temperature";
    Medium.DerDensityByPressure dddp[N] "Derivative of density by pressure";
    // Real dddX[N, nX](each unit="kg/m3") "Derivative of density by composition";           
  equation
    assert(FFtype == ThermoPower.Choices.Flow1D.FFtypes.NoFriction or dpnom > 0,
      "dpnom=0 not supported, it is also used in the homotopy trasformation during the inizialization");  
    //All equations are referred to a single tube
    // Friction factor selection
    omega_hyd = 4*A/Dhyd;
    if FFtype == FFtypes.Kfnom then
      Kf = Kfnom*Kfc;
      Cf = 2*Kf*A^3/(omega_hyd*L);
    elseif FFtype == FFtypes.OpPoint then
      Kf = dpnom*rhonom/(wnom/Nt)^2*Kfc;
      Cf = 2*Kf*A^3/(omega_hyd*L);
    elseif FFtype == FFtypes.Cfnom then
      Kf = Cfnom*omega_hyd*L/(2*A^3)*Kfc;
      Cf = Cfnom*Kfc;
    elseif FFtype == FFtypes.Colebrook then
      Cf = f_colebrook(
          w,
          Dhyd/A,
          e,
          Medium.dynamicViscosity(gas[integer(N/2)]))*Kfc;
    elseif FFtype == ThermoPower.Choices.Flow1D.FFtypes.NoFriction then
      Cf = 0;
      Kf = 0;
    else
      assert(false, "Unsupported FFtype");
      Cf = 0;
      Kf = 0;
    end if;
    // Kf = Cf*omega_hyd*L/(2*A^3)   "Relationship between friction coefficient and Fanning friction factor";
    assert(Kf >= 0, "Negative friction coefficient");
    Kfl = wnom/Nt*wnf*Kf "Linear friction factor";
    
    // Dynamic momentum term
    if DynamicMomentum then
      dwdt = der(w);
    else
      dwdt = 0;
    end if;

    sum(dMdt) = (infl.m_flow + outfl.m_flow)/Nt "Mass balance";
    L/A*dwdt + (outfl.p - infl.p) + Dpfric = 0 "Momentum balance";
    Dpfric = (if FFtype == ThermoPower.Choices.Flow1D.FFtypes.NoFriction then 0
            else homotopy((smooth(1, Kf*squareReg(w, wnom/Nt*wnf))*sum(vbar)/(N - 1)),
                            dpnom/(wnom/Nt)*w))
    "Pressure drop due to friction";
    
    for j in 1:N-1 loop
      if not QuasiStatic then
        // "Energy balance"
        
        if Medium.singleState then
          
          A*l*rhobar[j]*cvbar[j]*der(Ttilde[j]) + wbar[j]*(h[j + 1] - h[j]) = Q_single[j]
            "(pressure effects neglected)";
            
        else
          A*l*rhobar[j]*cvbar[j]*der(Ttilde[j]) + wbar[j]*(h[j + 1] - h[j]) - A*l*der(p) =
            Q_single[j] "Energy balance";  
            
        end if;
        
        // dMdt[j] = A*l*(drbdh[j]*der(htilde[j]) + drbdp[j]*der(p))        
        // Mass balance Equation in Gas.Flow1DFV. It is gave up for using gas.T as derivatived variable
        // which may cause unbalanced equations. Use independent variables such as Ttilde instead. 
        // dMdt[j] = A*l*(drbdp[j]*der(p) + drbdT1[j]*der(gas[j].T) + drbdT2[j]*der(gas[j + 1].T)) "Mass balance";
        
        dMdt[j] = A*l*(drbdT1[j]*der(Ttilde[j-1]) + drbdT2[j]*der(Ttilde[j]) + drbdp[j]*der(p))
          "Mass derivative for each volume";
          
        if avoidInletEnthalpyDerivative and j == 1 then
          // first volume properties computed by the volume outlet properties
          rhobar[j] = gas[j + 1].d;
          drbdp[j] = dddp[j + 1];
          drbdT1[j] = 0;
          drbdT2[j] = dddT[j + 1];
          // drbdX1[j, :] = zeros(size(Xtilde, 2));
          // drbdX2[j, :] = dddX[j + 1, :];
        else
          // volume properties computed by averaging
          rhobar[j] = (gas[j].d + gas[j + 1].d)/2;
          drbdp[j] = (dddp[j] + dddp[j + 1])/2;
          drbdT1[j] = dddT[j]/2;
          drbdT2[j] = dddT[j + 1]/2;
          // drbdX1[j, :] = dddX[j, :]/2;
          // drbdX2[j, :] = dddX[j + 1, :]/2;
        end if;
                
        // Average volume quantities
        // rhobar[j] = (rho[j] + rho[j + 1])/2;
        // drbdp[j] = (drdp[j] + drdp[j + 1])/2;
        drbdh[j] = (drdh[j] + drdh[j + 1])/2;
        vbar[j] = 1/rhobar[j];
        if fixedMassFlowSimplified then
          wbar[j] = homotopy(infl.m_flow/Nt - sum(dMdt[1:j - 1]) - dMdt[j]/2, wnom/Nt);
        else
          wbar[j] = infl.m_flow/Nt - sum(dMdt[1:j - 1]) - dMdt[j]/2;
        end if;
        cvbar[j] = (cv[j] + cv[j + 1])/2;
      else
        // Static mass and energy balances
        wbar[j]*(gas[j + 1].h - gas[j].h) = Q_single[j] "Energy balance";
        dMdt[j] = 0 "Mass balance";
        // Dummy values for unused average quantities
        rhobar[j] = 0;
        drbdp[j] = 0;
        drbdT1[j] = 0;
        drbdT2[j] = 0;
        // drbdX1[j, :] = zeros(nX);
        // drbdX2[j, :] = zeros(nX);
        vbar[j] = 0;
        wbar[j] = infl.m_flow/Nt;
        cvbar[j] = 0;      
      end if;
    end for;

    // for j in 1:N loop
    //   wstar[j] = homotopy(infl.m_flow/Nt - sum(dMdt[1:j - 1]) - dMdt[j]/2, wnom/Nt);
    // end for;

    // Fluid property calculations
    for j in 1:N loop
      gas[j] = Medium.setState_ph(p, h[j]);
      T[j] = Medium.temperature(gas[j]);
      rho[j] = Medium.density(gas[j]);
      drdp[j] = if Medium.singleState then 0 else Medium.density_derp_h(gas[j]);
      drdh[j] = Medium.density_derh_p(gas[j]);
      u[j] = w/(rho[j]*A);
    end for;
    
    // Fluid property computations
    for j in 1:N loop
      if not QuasiStatic then
        cv[j] = Medium.heatCapacity_cv(gas[j]);
        dddT[j] = Medium.density_derT_p(gas[j]);
        dddp[j] = Medium.density_derp_T(gas[j]);
        /*
        if nX > 0 then
          dddX[j, :] = Medium.density_derX(gas[j].state);
        end if;
        */
      else
        // Dummy values (not needed by dynamic equations)
        cv[j] = 0;
        dddT[j] = 0;
        dddp[j] = 0;
        // dddX[j, :] = zeros(nX);
      end if;
    end for;    

    // Boundary conditions
    win = infl.m_flow/Nt;
    wout = -outfl.m_flow/Nt;
    assert(HydraulicCapacitance == HCtypes.Upstream or
             HydraulicCapacitance == HCtypes.Downstream,
             "Unsupported HydraulicCapacitance option");
    /*
    if HydraulicCapacitance == HCtypes.Middle then
      p = infl.p - Dpfric1 - Dpstat/2;
      w = win;
    else
    */
    if HydraulicCapacitance == HCtypes.Upstream then
      p = infl.p;
      w = -outfl.m_flow/Nt;
    else  // if HydraulicCapacitance == HCtypes.Downstream then
      p = outfl.p;
      w = win;
    end if;
    infl.h_outflow = h[1];
    outfl.h_outflow = h[N]; //htilde[N - 1];

    h[1] = inStream(infl.h_outflow);
    // h[2:N] = htilde;
    T[2:N] = Ttilde;

    connect(wall,heatTransfer.wall);

    Q = heatTransfer.Q "Total heat flow through lateral boundary";
    M = sum(rhobar)*A*l "Fluid mass (single tube)";
    Mtot = M*Nt "Fluid mass (total)";
    Tr = noEvent(M/max(win, Modelica.Constants.eps)) "Residence time";

    assert(w > -wnom*wnm, "Reverse flow not allowed, maybe you connected the component with wrong orientation");
  initial equation
    if initOpt == ThermoPower.Choices.Init.Options.noInit then
      // do nothing
    elseif initOpt == ThermoPower.Choices.Init.Options.fixedState then
      if not noInitialPressure then
        p = pstart;
      end if;
      // htilde = hstart[2:N];
      Ttilde = Tstart[2:N];
    elseif initOpt == ThermoPower.Choices.Init.Options.steadyState then
      // der(htilde) = zeros(N - 1);
      der(Ttilde) = zeros(N - 1);
      if (not Medium.singleState) and not noInitialPressure then
        der(p) = 0;
      end if;
    elseif initOpt == ThermoPower.Choices.Init.Options.steadyStateNoP then
      // der(htilde) = zeros(N - 1);
      der(Ttilde) = zeros(N - 1);
      assert(false, "initOpt = steadyStateNoP deprecated, use steadyState and noInitialPressure",AssertionLevel.warning);
    elseif initOpt == ThermoPower.Choices.Init.Options.steadyStateNoT and not Medium.singleState then
      der(p) = 0;
    else
      assert(false, "Unsupported initialisation option");
    end if;
end Flow1DFV;
