within Steps.Components;

record HXCellState
  "Record to store the status parameter of a cell in heat exchanger"
  
  import SI = Modelica.SIunits;
  
  //Local temperature
  Modelica.SIunits.Temperature T_c;
  Modelica.SIunits.Temperature T_h;
  
  //Local pressure
  Modelica.SIunits.Pressure p_c;
  Modelica.SIunits.Pressure p_h;
  
  // Local velocity of fluid
  Modelica.SIunits.Velocity u_c;
  Modelica.SIunits.Velocity u_h;
  
  // mass flow rate
  Modelica.SIunits.MassFlowRate mdot_c;
  Modelica.SIunits.MassFlowRate mdot_h;
  
  //local parameters of this cell listed as following
  //Local Dynamic Viscosity
  Modelica.SIunits.DynamicViscosity mu_c;
  Modelica.SIunits.DynamicViscosity mu_h;
      
  // Local Conductivity
  Modelica.SIunits.ThermalConductivity k_c;
  Modelica.SIunits.ThermalConductivity k_h;
  
  // Local Reynolds Number
  Modelica.SIunits.ReynoldsNumber Re_c;
  Modelica.SIunits.ReynoldsNumber Re_h;
  
  // Local Density
  Modelica.SIunits.Density rho_c;
  Modelica.SIunits.Density rho_h;
  
  // Local Nusselt Number
  Modelica.SIunits.NusseltNumber Nu_c;
  Modelica.SIunits.NusseltNumber Nu_h;
  
  // Local PrandtlNumber
  Modelica.SIunits.PrandtlNumber Pr_c;
  Modelica.SIunits.PrandtlNumber Pr_h;
  
  // local Thermal Conductance
  Modelica.SIunits.CoefficientOfHeatTransfer h_c;
  Modelica.SIunits.CoefficientOfHeatTransfer h_h;
  
  // Fanning Friction Factor - used to calculate pressure drop
  Real f_c; 
  Real f_h; 
  
  // Pressure drop
  Modelica.SIunits.PressureDifference dp_c;
  Modelica.SIunits.PressureDifference dp_h;
  
  // specific enthalpy to cal Heat flux
  Modelica.SIunits.SpecificEnthalpy h_mass_c;
  Modelica.SIunits.SpecificEnthalpy h_mass_h;
  
  String medium_name_c;
  String medium_name_h;
  
  // length of this cell
  Modelica.SIunits.Length length "unit m";   
  
  // Heat Flux
  Modelica.SIunits.HeatFlux q;  
  
  // wall thermal conductivity - determined by material of wall and local temperature
  Modelica.SIunits.ThermalConductivity k_wall;
  
  // overall Heat transfer coefficient
  Modelica.SIunits.CoefficientOfHeatTransfer U;
 
end HXCellState;
