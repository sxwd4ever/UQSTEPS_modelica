within Steps.Components;

record HXState
  "Record to store the status parameter of a cell in heat exchanger"
  
  //Modelica.Blocks.Types.ExternalCombiTable1D table_th_inconel_750;
  
  String name_material;
  
  Modelica.SIunits.Length t_wall "thickness of wall between two neighboring hot and cold";
  
  Modelica.SIunits.Area A_stack "surface area of all cells in a stack"; 
  
  Modelica.SIunits.Diameter d_h "Hydraulic Diameter";
  
  Real fit_const_a "fitting constant a in Eq[3] of [kim, 2011] ";
  
  Real fit_const_b "fitting constant b in Eq[3] of [kim, 2011] ";
  
  Real fit_const_c "fitting constant c in Eq[4] of [kim, 2011] ";
  
  Real fit_const_d "fitting constant d in Eq[4] of [kim, 2011] ";
  
  Modelica.SIunits.Area A_flow "Flow area of all channels";
 
end HXState;
