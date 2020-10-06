#ifndef MY_RPOPS_H__
#define MY_RPOPS_H__

#ifdef __cplusplus
// put it here to avoid error: template with C linkage
#include <string> 
#include "PCHE.h"

extern "C" {
#endif


// #define _GLIBCXX_USE_CXX11_ABI 0

#ifdef BUILDING_DLL
#define EXPORT_MY_CODE __declspec(dllexport)
#else
#define EXPORT_MY_CODE __declspec(dllimport)
#endif

void EXPORT_MY_CODE MyPropsSI_pT(double p, double T, const std::string &FluidName, double &h, double &rhomass);

double EXPORT_MY_CODE MyPropsSI_pH(double p, double H, const char * FluidName, double &mu, double &k, double &rho);

/**
 * off-design simulation for PCHE
 */
double EXPORT_MY_CODE PCHE_OFFD_Simulation(PCHE_GEO_PARAM * geo, KIM_CORR_COE * cor, SIM_PARAM * sim_param, ThermoState * st_hot_in, ThermoState * st_cold_in, ThermoState * st_hot_out, ThermoState * st_cold_out);

ThermoState * EXPORT_MY_CODE NewThermoState_pT(double p, double T, double mdot, std::string medium);

double EXPORT_MY_CODE from_deg(double deg);

/** 
 * convert preseaure in bar to Pa
 */
double EXPORT_MY_CODE from_bar(double p_bar);

double EXPORT_MY_CODE from_degC(double degC);

bool EXPORT_MY_CODE the_same(double x, double y, double eps, double & diff_per);

#ifdef __cplusplus
}
#endif

#endif