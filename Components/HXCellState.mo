within Steps.Components;

record HXCellState
  "Record to store the status parameter of a cell in heat exchanger"
  
  import SI = Modelica.SIunits;
  
  // length of this cell
  Modelica.SIunits.Length length = 1.0e-3 "unit m"; 
  
  //Local temperature
  Modelica.SIunits.Temperature T_h = 0.0;
  Modelica.SIunits.Temperature T_c = 0.0;

  //Local pressure
  Modelica.SIunits.Pressure p_h = 0.0;
  Modelica.SIunits.Pressure p_c = 0.0;

  // Local velocity of fluid
  Modelica.SIunits.Velocity u_h = 0.0;
  Modelica.SIunits.Velocity u_c = 0.0;
  
  // Heat Flux
  Modelica.SIunits.HeatFlux q = 0.0;
  
  //local parameters of this cell listed as following
  //Local Dynamic Viscosity
  Modelica.SIunits.DynamicViscosity mu_h = 0;
  Modelica.SIunits.DynamicViscosity mu_c = 0;
  
  // Local Conductivity
  Modelica.SIunits.ThermalConductivity k_h = 0;
  Modelica.SIunits.ThermalConductivity k_c = 0;
  
  // Local Reynolds Number
  Modelica.SIunits.ReynoldsNumber Re_h = 0.0;
  Modelica.SIunits.ReynoldsNumber Re_c = 0.0;
  
  // Local Density
  Modelica.SIunits.Density rho_h = 0.0;
  Modelica.SIunits.Density rho_c = 0.0;
  
  // Local Nusselt Number
  Modelica.SIunits.NusseltNumber Nu_h = 0.0;
  Modelica.SIunits.NusseltNumber Nu_c = 0.0;
  
  // Local PrandtlNumber
  Modelica.SIunits.PrandtlNumber Pr_h = 0.0;
  Modelica.SIunits.PrandtlNumber Pr_c = 0.0;
  
  // wall thermal conductivity - determined by material of wall and local temperature
  Modelica.SIunits.ThermalConductivity k_wall = 0;
  
  // local Thermal Conductance
  Modelica.SIunits.CoefficientOfHeatTransfer h_h = 0;
  Modelica.SIunits.CoefficientOfHeatTransfer h_c = 0;
  
  // overall Heat transfer coefficient
  Modelica.SIunits.CoefficientOfHeatTransfer U = 0;
  
  // Fanning Friction Factor - used to calculate pressure drop
  Real f_h = 0.0;
  Real f_c = 0.0; 
  
  // Pressure drop
  Modelica.SIunits.PressureDifference dp_h = 0.0;
  Modelica.SIunits.PressureDifference dp_c = 0.0;
  
  // specific enthalpy to cal Heat flux
  Modelica.SIunits.SpecificEnthalpy h_mass_h = 0.0;
  Modelica.SIunits.SpecificEnthalpy h_mass_c = 0.0;
 
end HXCellState;
