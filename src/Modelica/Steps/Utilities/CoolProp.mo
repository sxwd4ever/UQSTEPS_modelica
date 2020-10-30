within Steps.Utilities;

model CoolProp

  import Steps.Components.{ThermoState,PCHEGeoParam,KimCorrelations, SimParam, PCHEBoundaryCondition, SimulationResult, PCHECImplResult};
  import Modelica.SIunits.Conversions.{from_deg,from_bar};

  function PropsSI "Import the MyPropsLib's PropsSI Function"
    input String OutputType;
    input String Input1Type;
    input Real Input1Value;
    input String Input2Type;
    input Real Input2Value;
    input String FluidName;
    output Real OutputValue;
  
    external "C" OutputValue = MyPropsSI(OutputType, Input1Type, Input1Value, Input2Type, Input2Value, FluidName) annotation(
      Library = {"MyProps"},
      LibraryDirectory = "modelica://Steps/Resources/Library");
  end PropsSI;

  /*
    void MyPropsSI_pH(double p, double H, const std::string &FluidName, double &T, double &mu, double &k, double &rho);
    */

  function MyPropsSI "Import the CoolProp PropsSI Function"
    input Real p;
    input Real H;
    input String fluidName;
    output Real T;
    output Real mu;
    output Real k;
    output Real rho;
  
    external "C" T = MyPropsSI_pH(p, H, fluidName, mu, k, rho);
    annotation(
      Library = {"MyProps"},
      LibraryDirectory = "modelica://Steps/Resources/Library");
  end MyPropsSI;

  function MyPropsSI_pT "Import the CoolProp PropsSI Function"
    input Real p;
    input Real T;
    input String fluidName;
    output Real h;
    output Real rho;
  
    external "C" MyPropsSI_pT(p, T, fluidName, h, rho);
    annotation(
      Library = {"MyProps"},
      LibraryDirectory = "modelica://Steps/Resources/Library");
  end MyPropsSI_pT;

  /*
   * double EXPORT_MY_CODE PCHE_OFFD_Simulation_UQ_out(const char * name, const char * media_hot, const char * media_cold, PCHE_GEO_PARAM * geo, KIM_CORR_COE * cor, SIM_PARAM * sim_param, BoundaryCondtion * bc, PCHECImplResult * retOutput, double * Q, double * U)
   */

  function PCHE_OFFD_Simulation_UQ_out "test for transferring c struct as input/output parameter"
    input String pche_name;
    input String media_hot;
    input String media_cold;
    input PCHEGeoParam geo;
    input KimCorrCoe cor;
    input SimParam sim_param;
    input PCHEBoundaryCondition bc;
    output PCHECImplResult retOutput;
    output Real Q[200];
    output Real U[200];
  protected
    // mapping of array of structs or structs with array is not allowed.
    // so I have to use these 'ugly' approach
    // size of array should be specified as fixed value, could be larger than the real one
    // output Real st[10];
  
    external "C" PCHE_OFFD_Simulation_UQ_out(pche_name, media_hot, media_cold, geo, cor, sim_param, bc, retOutput, Q, U);
    annotation(
      Library = {"MyProps"},
      LibraryDirectory = "modelica://Steps/Resources/Library");
  end PCHE_OFFD_Simulation_UQ_out;
  
  /*
double EXPORT_MY_CODE PCHE_OFFD_Simulation(const char * name, const char * media_hot, const char * media_cold, PCHE_GEO_PARAM * geo, KIM_CORR_COE * cor, SIM_PARAM * sim_param, BoundaryCondtion * bc, PCHECImplResult * retOutput)
  */
  
  function PCHE_OFFD_Simulation "test for transferring c struct as input/output parameter"
    input String pche_name;
    input String media_hot;
    input String media_cold;
    input PCHEGeoParam geo;
    input KimCorrelations.KimCorrCoe cor;
    input Model.SimParam sim_param;
    input PCHEBoundaryCondition bc;
    output PCHECImplResult retOutput;
  protected
    // mapping of array of structs or structs with array is not allowed.
    // so I have to use these 'ugly' approach
    // size of array should be specified as fixed value, could be larger than the real one
    // output Real st[10];
  
    external "C" PCHE_OFFD_Simulation(pche_name, media_hot, media_cold, geo, cor, sim_param, bc, retOutput);
    annotation(
      Library = {"MyProps"},
      LibraryDirectory = "modelica://Steps/Resources/Library");
  end PCHE_OFFD_Simulation;
  
  /*
   * 
  void EXPORT_MY_CODE print_path_state(const char * name, const char * media, ThermoState * st, int log_level);
  */ 
  function PrintPathState
    input String name;
    input String medium;
    input Model.ThermoState st; 
    input Integer log_level = 1;    
    
    output Real result;
    external "C" result = print_path_state(name, medium, st, log_level);
    annotation(
      Library = {"MyProps"},
      LibraryDirectory = "modelica://Steps/Resources/Library");      
      //connect(from, to);
  end PrintPathState;
  
  /*
   * 
  double EXPORT_MY_CODE test_struct_param(SIM_PARAM * sim_para, PCHE_GEO_PARAM * geo, BoundaryCondtion * bc)
  */  
  function TestStructParam "test for transferring c struct as input/output parameter"
    input SimParam sim_param;
    input PCHEGeoParam geo;
    input PCHEBoundaryCondition bc;
    output Real h_hot[2];
    output Real h_cold[2];
    output Real p_hot[2];
    output Real p_cold[2];
  protected
    // mapping of array of structs or structs with array is not allowed.
    // so I have to use these 'ugly' approach
    // size of array should be specified as fixed value, could be larger than the real one
    // output Real st[10];
  
    external "C" test_struct_param(sim_param, geo, bc, h_hot, h_cold, p_hot, p_cold);
    annotation(
      Library = {"MyProps"},
      LibraryDirectory = "modelica://Steps/Resources/Library");
  end TestStructParam;

  function setState "Demo external function"
    input Real p;
    input Real M;
    output State state;
  
    external "C" setState_C_impl(p, M, state) annotation(
      Library = {"MyProps"},
      LibraryDirectory = "modelica://Steps/Resources/Library");
  end setState;
  
  function QueryProps "using exteranl object to query props"
  
    input CoolPropExternalObject cp_ex;
    input String input_pair "PT_INPUTS/HmassP_INPUTS/PSmass_INPUTS see DataStructures.h in CoolProp";
    input Real val1;
    input Real val2;
    input String output_name;
    output Real props;
    
    //double EXPORT_MY_CODE cp_query(void * wrapper, const char * input_pair,  double val1, double val2, const char * output_name);
    
    external "C" props = cp_query(cp_ex, input_pair, val1, val2, output_name);
    annotation(
      Library = {"MyProps"},
      LibraryDirectory = "modelica://Steps/Resources/Library");
  end QueryProps;
  
  /*
  PCHEGeoParam geo(pitch = 12e-3, phi = from_deg((180 - 108) / 2), length = 1000e-3, d_c = 1.5e-3, N_ch = 80000, N_seg = 200);
  SimParam sim_param(err = 1e-2, delta_T_init = 5, N_iter = 10000, step_rel = 0.2, log_level = 1);
  PCHEBoundaryCondition bc(st_hot_in(T = 730, p = from_bar(90), h = 932534, mdot = 10), st_cold_in(T = 500, p = from_bar(225), h = 627426, mdot = 10), st_hot_out(T = 576.69, p = from_bar(90), h = 754560, mdot = 10), st_cold_out(T = 639.15, p = from_bar(225), h = 805341, mdot = 10));
  PCHECImplResult retOutput;
  Real Q[200];
  Real U[200];
  KimCorrCoe cor(a = 0.37934, b = 0.82413, c = 0.03845, d = 0.73793);
  //PCHESimExternalObject pche_ext_obj = PCHESimExternalObject(geo = geo);
  Real p = 20;
  Real M = 30;
  State state;
  */
  CoolPropExternalObject cp_wrapper = CoolPropExternalObject("CO2", "test");
  Real prop1, prop2, prop3;
equation
// sr = TestStructParam(sim_param, geo, bc);
// (h_hot, h_cold, p_hot, p_cold) = TestStructParam(sim_param, geo, bc, geo.N_seg);
// (retOutput) = PCHE_OFFD_Simulation("CoolPCHE", "CO2", "CO2", geo, cor, sim_param, bc);
// (retOutput, Q, U) = PCHE_OFFD_Simulation_UQ_out("CoolPCHE", "CO2", "CO2", geo, cor, sim_param, bc);
// state = setState(p);
// state = setState(p, M);
// Test if this interface works
/*
  prop1 = QueryProps(cp_wrapper, "PT_INPUTS", 20e6, from_degC(700), "H");
  prop2 = QueryProps(cp_wrapper, "PT_INPUTS", 8e6, from_degC(700), "H");
  prop3 = QueryProps(cp_wrapper, "PT_INPUTS", 8e6, from_degC(35), "H");
*/

  prop1 = PropsSI("H", "P", 20e6, "T", from_degC(700), "CO2");
  prop2 = PropsSI("H", "P", 8e6, "T", from_degC(700), "CO2");
  prop3 = PropsSI("H", "P", 8e6, "T", from_degC(35), "CO2");
/*
  (T, mu, k , rho) = MyPropsSI(p, H, "CO2");
  //(H, rho) = MyPropsSI_pT(p, T, "CO2");

  T = PropsSI("T", "P", p, "H", H, "CO2");
  mu = PropsSI("V", "P", p, "H", H, "CO2");
  k = PropsSI("L", "P", p, "H", H, "CO2");
  rho = PropsSI("D", "P", p, "H", H, "CO2");
*/

/*
void EXPORT_MY_CODE print_path_state(const char * name, const char * media, ThermoState * st_in, ThermoState * st_out, int log_level);
*/



annotation(
  experiment(StartTime = 0, StopTime = 1, Interval = 1, Tolerance = 1e-6),
    __OpenModelica_commandLineOptions = "--matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection -d=initialization,NLSanalyticJacobian,aliasConflicts");  
end CoolProp;
