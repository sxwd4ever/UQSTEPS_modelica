within Steps.Components;

record PCHEGeoParam    
  "Record for PCHEGeoParam, align with the struct definition in C"
  extends Model.EntityGeoParam;    
    
    // pitch length
    parameter Real pitch;
    // pitch angle
    parameter Real phi;    
    parameter Real length = L;    
    parameter Real d_c = d;
  
end PCHEGeoParam;
