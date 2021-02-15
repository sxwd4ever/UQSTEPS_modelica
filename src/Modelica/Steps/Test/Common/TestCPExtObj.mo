within Steps.Test.Common;

model TestCPExtObj
  import Steps.Utilities.CoolProp;
  import Steps.Utilities.CoolPropExternalObject;
  
  model Context
    
    CoolPropExternalObject cp_ext_obj = CoolPropExternalObject("CO2", "test");  
  
  end Context;

  Context ctx;  
  
  Real prop1, prop2;
  
equation

  u = Modelica.Math.sin(time);
  u_pre = previous(u);
  u_pre_2 = previous(u_pre);  

  y = u > 0.5 or pre(y) and u >= -0.5;
  
annotation(
  experiment(StartTime = 0, StopTime = 10, Interval = 0.1, Tolerance = 1e-6),
    __OpenModelica_commandLineOptions = "--matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection -d=initialization,NLSanalyticJacobian,aliasConflicts");  
    
      
end TestCPExtObj;
