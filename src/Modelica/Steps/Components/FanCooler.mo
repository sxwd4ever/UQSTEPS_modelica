within Steps.Components;

model FanCooler
  extends TwoPorts;
  
  import Steps.Interfaces.PortType;
  import Steps.Utilities.CoolProp.PropsSI;
  
  parameter Modelica.SIunits.TemperatureDifference delta_T = 18.0 "Incoming temperature difference";  
  
  parameter Modelica.SIunits.Temperature T_amb = Modelica.SIunits.Conversions.from_degC(15) "ambinent temperature";
  
  Modelica.SIunits.SpecificEnthalpy h_out;
/*
initial algorithm

    if init_inlet then
      inlet.p := p_init_inlet;
      inlet.h_outflow := PropsSI("H", "P", p_init_inlet, "T", T_init_inlet, PBMedia.mediumName); 
    end if;
    
    if init_outlet then
      outlet.p := p_init_outlet;
      outlet.h_outflow := PropsSI("H", "P", p_init_outlet, "T", T_init_outlet, PBMedia.mediumName);
    end if;
*/
equation

  // Note that, no matter what the CO2 inlet conditions are, the fan colder CO2 
  // exit temperature is always the same for fixed ITD and ambient air temperature.  
  h_out = PropsSI("H", "P", inlet.p, "T", T_amb + delta_T, PBMedia.mediumName);
  
  outlet.p = inlet.p;
  outlet.m_flow + inlet.m_flow = 0;
  outlet.h_outflow = h_out;
  
  // outlet.m_flow = - mdot_init;
  
  inlet.h_outflow = inStream(inlet.h_outflow); 

end FanCooler;
