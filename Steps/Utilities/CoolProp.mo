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
    annotation(Library={"CoolProp"}, LibraryDirectory="modelica://Steps/Resources/Library", IncludeDirectory="modelica://Steps/Resources/Include"); 
    
  end PropsSI;
end CoolProp;
