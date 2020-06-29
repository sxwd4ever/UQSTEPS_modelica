within Steps.Components;

record HXCellState
  "Record to store the status parameter of a cell in heat exchanger"
  
  import SI = Modelica.SIunits;
  
  record CellStatus
    "Cell Status for one cell in the hot-cell and cool-cell pair"
    
    //Local temperature
    Modelica.SIunits.Temperature T;
  
    //Local pressure
    Modelica.SIunits.Pressure p;
  
    // Local velocity of fluid
    Modelica.SIunits.Velocity u;
    
    // mass flow rate
    Modelica.SIunits.MassFlowRate mdot;
    
    //local parameters of this cell listed as following
    //Local Dynamic Viscosity
    Modelica.SIunits.DynamicViscosity mu;
        
    // Local Conductivity
    Modelica.SIunits.ThermalConductivity k;
    
    // Local Reynolds Number
    Modelica.SIunits.ReynoldsNumber Re;
    
    // Local Density
    Modelica.SIunits.Density rho;
    
    // Local Nusselt Number
    Modelica.SIunits.NusseltNumber Nu;
    
    // Local PrandtlNumber
    Modelica.SIunits.PrandtlNumber Pr;
    
    // local Thermal Conductance
    Modelica.SIunits.CoefficientOfHeatTransfer h;
    
    // Fanning Friction Factor - used to calculate pressure drop
    Real f; 
    
    // Pressure drop
    Modelica.SIunits.PressureDifference dp;
    
    // specific enthalpy to cal Heat flux
    Modelica.SIunits.SpecificEnthalpy h_mass;
    
    String medium_name;
    
  end CellStatus;
  
  // length of this cell
  Modelica.SIunits.Length length "unit m";   
  
  // Heat Flux
  Modelica.SIunits.HeatFlux q;  
  
  // wall thermal conductivity - determined by material of wall and local temperature
  Modelica.SIunits.ThermalConductivity k_wall;
  
  // overall Heat transfer coefficient
  Modelica.SIunits.CoefficientOfHeatTransfer U;
  
  // Status of the hot /cold cell in a pair
  
  CellStatus status[2];  
 
end HXCellState;
