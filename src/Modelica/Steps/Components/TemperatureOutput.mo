within Steps.Components;

model TemperatureOutput
  extends Modelica.Blocks.Interfaces.SO;
  
  import SI = Modelica.SIunits;
  
  parameter SI.Temperature T_start "Start temperature";
  
  parameter SI.Temperature T_stop "Stop temperature";
  
  parameter SI.TemperatureDifference dT_step "detla Temperature between iterations";
  
  parameter Real t_sim_duration "duration of the simulation";
  
  SI.Time startTime = 0 "Output y = offset for time < startTime";
  
  Real period := t_sim_duration / max_count "interval of iteration";
  
  Real max_count =  ceil((T_stop - T_start) / dT_step + 1) "iteration number";
  
  Real count "iteration count";
  
  SI.Time t_start "Start of simulation";  
  
  SI.Temperature T_cur ;
  
initial algorithm

  count := integer((time - startTime)/period);

  t_start := startTime + count*period;  
  
  T_cur := T_start;
  
algorithm
  
  when ceil((time - startTime)/period) > pre(count) then
    count := pre(count) + 1;
    T_cur := if count <= max_count then T_start + (count - 1) * dT_step else T_stop;// else T_pcms[max_count];
    t_start := time;
  end when;
   
equation
  y = T_cur;
  
end TemperatureOutput;
