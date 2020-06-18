within Steps.Components;

record HXCellState
  "Record to store the status parameter of a cell in heat exchanger"
  
  import SI = Modelica.SIunits;
  
  // length of this cell
  Modelica.SIunits.Length length "unit m"; 
  
  //Local temperature
  Modelica.SIunits.Temperature T_h;
  Modelica.SIunits.Temperature T_c;

  //Local pressure
  Modelica.SIunits.Pressure p_h;
  Modelica.SIunits.Pressure p_c;

  // Local velocity of fluid
  Modelica.SIunits.Velocity u_h;
  Modelica.SIunits.Velocity u_c;
  
  // Heat Flux
  Modelica.SIunits.HeatFlux q;
  
  //local parameters of this cell listed as following
  //Local Dynamic Viscosity
  Modelica.SIunits.DynamicViscosity mu_h;
  Modelica.SIunits.DynamicViscosity mu_c;
  
  // Local Conductivity
  Modelica.SIunits.ThermalConductivity k_h;
  Modelica.SIunits.ThermalConductivity k_c;
  
  // Local Reynolds Number
  Modelica.SIunits.ReynoldsNumber Re_h;
  Modelica.SIunits.ReynoldsNumber Re_c;
  
  // Local Density
  Modelica.SIunits.Density rho_h;
  Modelica.SIunits.Density rho_c;
  
  // Local Nusselt Number
  Modelica.SIunits.NusseltNumber Nu_h;
  Modelica.SIunits.NusseltNumber Nu_c;
  
  // Local PrandtlNumber
  Modelica.SIunits.PrandtlNumber Pr_h;
  Modelica.SIunits.PrandtlNumber Pr_c;
  
  // wall thermal conductivity - determined by material of wall and local temperature
  Modelica.SIunits.ThermalConductivity k_wall;
  
  // local Thermal Conductance
  Modelica.SIunits.CoefficientOfHeatTransfer h_h;
  Modelica.SIunits.CoefficientOfHeatTransfer h_c;
  
  // overall Heat transfer coefficient
  Modelica.SIunits.CoefficientOfHeatTransfer U;
  
  // Fanning Friction Factor - used to calculate pressure drop
  Real f_h;
  Real f_c; 
  
  // Pressure drop
  Modelica.SIunits.PressureDifference dp_h;
  Modelica.SIunits.PressureDifference dp_c;
  
  // specific enthalpy to cal Heat flux
  Modelica.SIunits.SpecificEnthalpy h_mass_h;
  Modelica.SIunits.SpecificEnthalpy h_mass_c;
 
end HXCellState;
