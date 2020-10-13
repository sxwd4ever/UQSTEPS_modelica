within Steps.Components;

record PCHEGeoParam
  "Record for PCHEGeoParam, align with the struct definition in C"
  
    /* data */
    // pitch length
    Real pitch;
    // pitch angle
    Real phi;
    // length of segment
    Real length_cell;
    // Diameter of semi_circular
    Real d_c;
    // number of channels
    Integer N_ch;
    // number of segments
    Integer N_seg; // [Kwon2019]'s maximum node number
  
end PCHEGeoParam;
