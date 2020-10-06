within Steps.Components;

model ScopeTest
  "Record to store the status parameter of a cell in heat exchanger"
  
  //Modelica.Blocks.Types.ExternalCombiTable1D table_th_inconel_750;
  
  inner Real var;
  
  model A
    outer Real var;
    
    model B
      outer Real var;
    end B;
  end A;
  
  A a[10];
 
end ScopeTest;
