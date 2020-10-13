within Steps.Components;

class PCHESimExternalObject
  " I see no use of this External Object since I cannot transfer data back throught it"
  extends ExternalObject;

  function constructor
    input PCHEGeoParam geo;
    output PCHESimExternalObject pche_ext_obj;
    //void * EXPORT_MY_CODE init_PCHE_sim_ext_object(PCHE_GEO_PARAM * geo);
    external "C" pche_ext_obj = init_PCHE_sim_ext_object(geo)  annotation(Library={"MyProps"}, LibraryDirectory="modelica://Steps/Resources/Library");

  end constructor;  

  function destructor
    input PCHESimExternalObject pche_ext_obj;

    external "C" close_PCHE_sim_ext_object(pche_ext_obj)  annotation(Library={"MyProps"}, LibraryDirectory="modelica://Steps/Resources/Library");
    
  end destructor;  
 
end PCHESimExternalObject;
