within Steps.Components;

record PCHEGeoParam
  "Record for PCHEGeoParam, align with the struct definition in C"
  
    /* data */
    // pitch length
    parameter Real pitch;
    // pitch angle
    parameter Real phi;
    // length of PCHE
    parameter Real length;
    // Diameter of semi_circular
    parameter Real d_c;
    // number of channels
    parameter Integer N_ch;
    // number of segments
    parameter Integer N_seg; // [Kwon2019]'s maximum node number
  
end PCHEGeoParam;
