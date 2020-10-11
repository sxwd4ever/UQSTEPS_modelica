within Steps.Utilities;

model CoolProp

  import Steps.Components.{ThermoState, PCHEGeoParam, KimCorrCoe, SimParam, BoundaryCondition, SimulationResult};
  import Modelica.SIunits.Conversions.{from_deg, from_bar};
  constant Integer MAX_SEG_NUM = 500;
  
  function PropsSI "Import the MyPropsLib's PropsSI Function"
    
    input String OutputType;
    input String Input1Type;
    input Real   Input1Value;
    input String Input2Type;
    input Real   Input2Value;
    input String FluidName;
    output Real  OutputValue;
  
    external "C" OutputValue = MyPropsSI(OutputType, Input1Type, Input1Value, Input2Type, Input2Value, FluidName)
    annotation(Library={"MyProps"}, LibraryDirectory="modelica://Steps/Resources/Library"); 
    
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
    annotation(Library={"MyProps"}, LibraryDirectory="modelica://Steps/Resources/Library"); 

  end MyPropsSI;
  
  function MyPropsSI_pT "Import the CoolProp PropsSI Function"
    
    input Real p;
    input Real T;
    input String fluidName;
    output Real h;
    output Real rho;
  
    external "C" MyPropsSI_pT(p, T, fluidName, h, rho);
    annotation(Library={"MyProps"}, LibraryDirectory="modelica://Steps/Resources/Library"); 

  end MyPropsSI_pT;    
/*
 * PCHE_OFFD_Simulation(const char * media_hot, const char * media_cold, PCHE_GEO_PARAM * geo, KIM_CORR_COE * cor, SIM_PARAM * sim_param, BoundaryCondtion * bc, double * h_hot, double * h_cold, double * p_hot, double * p_cold, size_t N_seg)
 */
function PCHE_OFFD_Simulation "test for transferring c struct as input/output parameter"
    
    input String media_hot;
    
    input String media_cold;
    
    input PCHEGeoParam geo;
    
    input KimCorrCoe cor;
    
    input SimParam sim_param;    
    
    input BoundaryCondition bc;
    
    input Integer N_seg;       
    
    output Real h_hot[MAX_SEG_NUM];
    output Real h_cold[MAX_SEG_NUM];
    output Real p_hot[MAX_SEG_NUM];
    output Real p_cold[MAX_SEG_NUM];  
  protected     

    // mapping of array of structs or structs with array is not allowed.
    // so I have to use these 'ugly' approach
    // size of array should be specified as fixed value, could be larger than the real one

    // output Real st[10]; 
    
    external "C" PCHE_OFFD_Simulation(media_hot, media_cold, geo, cor, sim_param, bc, h_hot, h_cold, p_hot, p_cold, N_seg);
    
    annotation(Library={"MyProps"}, LibraryDirectory="modelica://Steps/Resources/Library");
    
  end PCHE_OFFD_Simulation;

/*
 * 
double EXPORT_MY_CODE test_struct_param(SIM_PARAM * sim_para, PCHE_GEO_PARAM * geo, BoundaryCondtion * bc)
*/
  function TestStructParam "test for transferring c struct as input/output parameter"
    
    input SimParam sim_param;
    
    input PCHEGeoParam geo;
    
    input BoundaryCondition bc;
    
    input Integer N_seg;       
    
    output Real h_hot[MAX_SEG_NUM];
    output Real h_cold[MAX_SEG_NUM];
    output Real p_hot[MAX_SEG_NUM];
    output Real p_cold[MAX_SEG_NUM];  
  protected     

    // mapping of array of structs or structs with array is not allowed.
    // so I have to use these 'ugly' approach
    // size of array should be specified as fixed value, could be larger than the real one

    // output Real st[10]; 
    
    external "C" test_struct_param(sim_param, geo, bc, h_hot, h_cold, p_hot, p_cold, N_seg);
    
    annotation(Library={"MyProps"}, LibraryDirectory="modelica://Steps/Resources/Library");
    
  end TestStructParam;
  
  function setState
    "Demo external function"
    
    input Real p;
    input Real M;
    output State state;
 
    external "C" setState_C_impl(p, M, state) annotation(Library={"MyProps"}, LibraryDirectory="modelica://Steps/Resources/Library");

  end setState;
  
  PCHEGeoParam geo(pitch = 12e-3, phi = from_deg((180-108)/2), length_cell = 5e-2, d_c = 1.5e-3,
  N_ch = 80000, N_seg = 200);
  
  SimParam sim_param(err=1e-2, delta_T_init = 5, N_iter = 10000, step_rel=0.2);
  
  BoundaryCondition bc(
  st_hot_in(T = 730, p = from_bar(90), h = 932534, mdot = 10), 
  st_cold_in(T = 500, p = from_bar(225), h = 627426, mdot = 10), 
  st_hot_out(T = 576.69, p = from_bar(90), h = 754560, mdot = 10), 
  st_cold_out(T = 639.15, p = from_bar(225), h = 805341, mdot = 10));   
  
  Real h_hot[MAX_SEG_NUM];
  Real h_cold[MAX_SEG_NUM];
  Real p_hot[MAX_SEG_NUM];
  Real p_cold[MAX_SEG_NUM];    
  
  KimCorrCoe cor(a=0.37934, b=0.82413, c=0.03845, d=0.73793);
  
  //PCHESimExternalObject pche_ext_obj = PCHESimExternalObject(geo = geo);

  Real p = 20;
  Real M = 30;
  State state;
equation

  // sr = TestStructParam(sim_param, geo, bc);
  // (h_hot, h_cold, p_hot, p_cold) = TestStructParam(sim_param, geo, bc, geo.N_seg);

  
  (h_hot, h_cold, p_hot, p_cold) = PCHE_OFFD_Simulation("CO2", "CO2", geo, cor, sim_param, bc, geo.N_seg);
  // state = setState(p);
  state = setState(p,M);
  
// Test if this interface works
/*
  (T, mu, k , rho) = MyPropsSI(p, H, "CO2");
  //(H, rho) = MyPropsSI_pT(p, T, "CO2");

  T = PropsSI("T", "P", p, "H", H, "CO2");
  mu = PropsSI("V", "P", p, "H", H, "CO2");
  k = PropsSI("L", "P", p, "H", H, "CO2");
  rho = PropsSI("D", "P", p, "H", H, "CO2");
*/
end CoolProp;
