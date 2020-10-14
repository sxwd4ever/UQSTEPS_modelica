within Steps.Components;

model Regulator
  extends Steps.Components.TwoPorts;
  
  import Steps.Utilities.CoolProp.PropsSI;
  /*    
  parameter Boolean init_inlet = false;
  
  parameter Modelica.SIunits.Temperature T_init_inlet = 500;
  
  parameter Modelica.SIunits.Pressure p_init_inlet = 8e6;
  
  parameter Boolean init_outlet = false;
  */
  
  parameter Modelica.SIunits.Temperature T_init = 600;
  
  parameter Modelica.SIunits.Pressure p_init = 9e6; 
  
  parameter Modelica.SIunits.MassFlowRate m_flow_init;   
  
  //Modelica.Blocks.Interfaces.RealInput T_input;
  
  Modelica.Blocks.Sources.BooleanStep bool_step(startValue = true, startTime = 0.01);  
      
equation  
  
  if bool_step.y == true then
    /*
    if init_inlet then
      inlet.p = p_init_inlet;
      inlet.h_outflow = PropsSI("H", "P", p_init_inlet, "T", T_init_inlet, PBMedia.mediumName); 
    end if;
    
    if init_outlet then

    end if;
    */
    outlet.p = p_init;
    outlet.h_outflow = PropsSI("H", "P", p_init, "T", T_init, PBMedia.mediumName);
    outlet.m_flow = - m_flow_init;
  else
  
    outlet.p = inlet.p;  
    outlet.m_flow + inlet.m_flow = 0;     
    outlet.h_outflow = inlet.h_outflow;    
  end if;
  
  inlet.h_outflow = inStream(inlet.h_outflow);
   
end Regulator;
