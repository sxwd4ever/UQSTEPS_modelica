within Steps.Interfaces;

connector PBFluidPort_a
  extends Modelica.Fluid.Interfaces.FluidPort_a;
  Modelica.SIunits.ThermodynamicTemperature T "absolute temperature of the port";
  
  
  // Port type of the ports
  // To specify if the value(pressure, Temperature) of a port should be fixed by preset value or
  // determined by conncted port.
  parameter PortType PT = PortType.free;
  
  /*
    replaceable package Medium = Modelica.Media.Interfaces.PartialMedium
      "Medium model" annotation (choicesAllMatching=true);

    flow Medium.MassFlowRate m_flow
      "Mass flow rate from the connection point into the component";
    Medium.AbsolutePressure p "Thermodynamic pressure in the connection point";
    
    
    stream Medium.SpecificEnthalpy h_outflow
      "Specific thermodynamic enthalpy close to the connection point if m_flow < 0";
    stream Medium.MassFraction Xi_outflow[Medium.nXi]
      "Independent mixture mass fractions m_i/m close to the connection point if m_flow < 0";
    stream Medium.ExtraProperty C_outflow[Medium.nC]
      "Properties c_i/m close to the connection point if m_flow < 0";
    */
end PBFluidPort_a;
