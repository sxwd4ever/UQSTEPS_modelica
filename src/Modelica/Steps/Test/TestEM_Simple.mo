within Steps.Test;

model TestEM_Simple
  "Simple Tests for ExternalMedia the library"  
  
  extends Modelica.Icons.Example;
    replaceable package Medium = ExternalMedia.Media.TestMedium;

    Medium.ThermodynamicState state;
  equation
    state = Medium.setState_ph(1e5, 1e5 + 1e5*time);
  annotation (
    Diagram(graphics),
    experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-2, Interval = 1),
    __OpenModelica_commandLineOptions = "--matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection -d=initialization,NLSanalyticJacobian,aliasConflicts,bltdump",    
    __OpenModelica_simulationFlags(lv = "LOG_DEBUG,LOG_NLS,LOG_NLS_V,LOG_STATS,LOG_INIT,LOG_STDOUT, -w", outputFormat = "mat", s = "dassl", nls = "homotopy"));
end TestEM_Simple;


