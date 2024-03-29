within Steps.Test.Components;

model TestPCHECImpl "PCHE c-based implementation"
  import Modelica.SIunits.Conversions.{from_bar,from_deg,from_degC};
  import Steps.Utilities.CoolProp.PropsSI;
  import Steps.Components.{PCHEBoundaryCondition,PCHEGeoParam,ThermoState,SimParam};
  import Steps.Cycle.OffDPBParamSet;
  parameter OffDPBParamSet param(geo_HTR.N_seg = 50);
  parameter Boolean SourceFixed_hot = true;
  parameter Boolean SourceFixed_cold = true;
  // select the group of parameter set as input here
  parameter PCHEBoundaryCondition bc = param.bc_HTR;
  parameter PCHEGeoParam geo = param.geo_HTR;
  parameter SimParam sim_param = param.sim_param_def;
  Components.Source source_hot(p_outlet = bc.st_hot_in.p, T_outlet = bc.st_hot_in.T, mdot_init = bc.st_hot_in.mdot, fix_state = SourceFixed_hot);
  Components.Source source_cold(p_outlet = bc.st_cold_in.p, T_outlet = bc.st_cold_in.T, mdot_init = bc.st_cold_in.mdot, fix_state = SourceFixed_cold);
  // True if set its state as boundary condition
  Components.Sink sink_hot(p_inlet = bc.st_hot_out.p, T_inlet = bc.st_hot_out.T, mdot_init = bc.st_hot_out.mdot, fix_state = not SourceFixed_hot);
  Components.Sink sink_cold(p_inlet = bc.st_cold_out.p, T_inlet = bc.st_cold_out.T, mdot_init = bc.st_cold_out.mdot, fix_state = not SourceFixed_cold);
  Components.PCHECImpl pche(geo = geo, bc = bc, sim_param = sim_param);
  // step_rel will affect result's error and simulation speed
  // step_rel will affect result's error and simulation speed
equation
  connect(source_hot.outlet, pche.inlet_hot);
  connect(pche.outlet_hot, sink_hot.inlet);
  connect(source_cold.outlet, pche.inlet_cold);
  connect(pche.outlet_cold, sink_cold.inlet);
  annotation(
    experiment(StartTime = 0, StopTime = 1, Interval = 1, Tolerance = 1e-6),
    __OpenModelica_commandLineOptions = "--matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection -d=initialization,NLSanalyticJacobian,aliasConflicts,bltdump",
    __OpenModelica_simulationFlags(lv = "LOG_DEBUG,LOG_EVENTS_V,LOG_INIT,LOG_LS,LOG_LS_V,LOG_NLS,LOG_NLS_V,LOG_NLS_JAC,LOG_NLS_RES,LOG_RES_INIT,LOG_STATS", outputFormat = "mat", s = "dassl", nls = "homotopy"));
end TestPCHECImpl;
