within Steps.Components;

record SimParam
  "Record for the simulation param, align with the struct definition in C"
  
    // error tolerance
    Real err;
    // initial delta_T between T_hot_in and T_cold_out
    Real delta_T_init;
    // maximum iteration number
    Integer N_iter;
    // relative step length between two trial values of T_cold_out
    // T_cold_out[i] += step_rel * (T_bc_cold_in  - T_cold_in[i-1])
    Real step_rel;   
    
    // log level for output contorl 
    // enum log_level {DEBUG = 0, INFO = 1, ERR = 2, SERVE = 3, OFF = 4};
    Integer log_level;
  
end SimParam;
