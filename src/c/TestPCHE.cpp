#include "MyPropsLib.h"
#include "PCHE.h"
#include <time.h>

int main()
{	
    std::string info = string_format("Ready to go %d\n", 1);

    time_t t1, t2;
    t1 = clock();
    PCHE_GEO_PARAM geo;
    // set the correlation coefficients mannually
    /* data */
    geo.pitch = 12.3e-3;
    geo.phi = from_deg((180-108)/2);
    geo.length_cell = 5e-3;
    geo.d_c = 1.5e-3;
    geo.N_ch = 80e3;
    geo.N_seg = 200; // [Kwon2019]'s maximum node number 

    // set the correlation coefficients mannually for demo purpose
    // should be determined by pitch and phi of PCHE
    KIM_CORR_COE cor;

    cor.a = 0.37934;
    cor.b = 0.82413;
    cor.c = 0.03845;
    cor.d = 0.73793;

    double mdot_hot = 10, mdot_cold = 10;
    const char * mname_hot = "CO2";
    const char * mname_cold = "CO2";

    BoundaryCondtion bc;

    bc.st_hot_in = * NewThermoState_pT(from_bar(90), 350, mdot_hot, mname_hot);
    bc.st_cold_in = * NewThermoState_pT(from_bar(225), 300, mdot_cold, mname_cold);
    bc.st_hot_out = * NewThermoState_pT(from_bar(90), 576.69, mdot_hot, mname_hot);
    bc.st_cold_out = * NewThermoState_pT(from_bar(225), 639.15, mdot_cold, mname_cold);

    SIM_PARAM sim_param;
    sim_param.N_iter = 1e4;
    sim_param.err = 1e-2; // 1 percent
    sim_param.step_rel = 0.2;
    sim_param.delta_T_init = 5;
    sim_param.log_level = 0;
       
    double * h_hot = new double[geo.N_seg];
    double * h_cold = new double[geo.N_seg];
    double * p_hot = new double[geo.N_seg];
    double * p_cold = new double[geo.N_seg];

    //double EXPORT_MY_CODE PCHE_OFFD_Simulation(const char * name, const char * media_hot, const char * media_cold, PCHE_GEO_PARAM * geo, KIM_CORR_COE * cor, SIM_PARAM * sim_param, BoundaryCondtion * bc, double * h_hot, double * h_cold, double * p_hot, double * p_cold);
    PCHE_OFFD_Simulation("PCHEDemo", "CO2", "CO2", &geo, &cor, &sim_param, &bc, h_hot, h_cold, p_hot, p_cold);     

    t2 = clock();
    printf("All done! Elapsed time for user defined tests: %g s",(double)(t2-t1)/CLOCKS_PER_SEC);
    return 1;
}