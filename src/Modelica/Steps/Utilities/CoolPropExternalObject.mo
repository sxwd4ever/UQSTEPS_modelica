within Steps.Utilities;

class CoolPropExternalObject
  "  CoolProp query function, Using low level api function of CoolProp to accelerate"
  extends ExternalObject;

  function constructor
    input String medium_name;
    input String id;
    output CoolPropExternalObject cp_wrapper;
    // void * EXPORT_MY_CODE init_cp_wrapper(const char * medium);
    external "C" cp_wrapper = init_cp_wrapper(medium_name, id)  annotation(Library={"MyProps"}, LibraryDirectory="modelica://Steps/Resources/Library");

  end constructor;  

  function destructor
    input CoolPropExternalObject cp_wrapper;
    // void EXPORT_MY_CODE close_cp_wrapper( void * state);
    external "C" close_cp_wrapper(cp_wrapper)  annotation(Library={"MyProps"}, LibraryDirectory="modelica://Steps/Resources/Library");
    
  end destructor;  
 
end CoolPropExternalObject;
