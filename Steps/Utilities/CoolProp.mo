within Steps.Utilities;

model CoolProp
  function PropsSI "Import the CoolProp PropsSI Function"
    
    input String OutputType;
    input String Input1Type;
    input Real   Input1Value;
    input String Input2Type;
    input Real   Input2Value;
    input String FluidName;
    output Real  OutputValue;
  
    external "C" OutputValue = PropsSI(OutputType, Input1Type, Input1Value, Input2Type, Input2Value, FluidName)
    annotation(Library={"CoolProp"}, LibraryDirectory="modelica://Steps/Resources/Library"); 
    
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
    annotation(Library={"MyPropsLib"}, LibraryDirectory="modelica://Steps/Resources/Library"); 

  end MyPropsSI;
  
  function MyPropsSI_pT "Import the CoolProp PropsSI Function"
    
    input Real p;
    input Real T;
    input String fluidName;
    output Real h;
    output Real rho;
  
    external "C" MyPropsSI_pT(p, T, fluidName, h, rho);
    annotation(Library={"MyPropsLib"}, LibraryDirectory="modelica://Steps/Resources/Library"); 

  end MyPropsSI_pT;    

  parameter Real p = 101325;
  parameter Real H = 993771;
  //Real H;
  
  Real T;
  Real mu;
  Real k;
  Real rho;

equation
// Test if this interface works

  (T, mu, k , rho) = MyPropsSI(p, H, "CO2");
  //(H, rho) = MyPropsSI_pT(p, T, "CO2");
/*
  T = PropsSI("T", "P", p, "H", H, "CO2");
  mu = PropsSI("V", "P", p, "H", H, "CO2");
  k = PropsSI("L", "P", p, "H", H, "CO2");
  rho = PropsSI("D", "P", p, "H", H, "CO2");
*/
end CoolProp;
