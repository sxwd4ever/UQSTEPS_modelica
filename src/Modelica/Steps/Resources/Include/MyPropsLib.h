#ifndef MY_RPOPS_H__
#define MY_RPOPS_H__

#ifdef __cplusplus
// put it here to avoid error: template with C linkage
#include <stdio.h>
#include <string> 
#include <memory>
#include <stdexcept>
#include "PCHE.h"

extern "C" {
#endif


// #define _GLIBCXX_USE_CXX11_ABI 0

#ifdef BUILDING_DLL
#define EXPORT_MY_CODE __declspec(dllexport)
#else
#define EXPORT_MY_CODE __declspec(dllimport)
#endif

typedef struct {

    //! Temperature
    double p;
    double M;

} State;

struct CoolPropWrapper
{
    long handle_cp_state;

    //optional name for this state
    const char * name;
    /* for CP's error handle */
    long err_code;
    char * cp_err_buf;
    long buffer_size;
};

void my_log(std::string content, log_level level);
/**
 * Wrapper of CoolProp's PropsSI function. 
 */
double EXPORT_MY_CODE MyPropsSI(const char *Output, const char *Name1, double Prop1, const char *Name2, double Prop2, const char *Ref);

/**
 * batch Props query function using (p,T): one thermo state in, multiple outputs out to save time of calling CoolProp
 */
void EXPORT_MY_CODE MyPropsSI_pT(double p, double T, const std::string &FluidName, double &h, double &rhomass);

/**
 * batch Props query function using (p,H), one thermo state in, multiple outputs out to save time of calling CoolProp
 */
double EXPORT_MY_CODE MyPropsSI_pH(double p, double H, const char * FluidName, double &mu, double &k, double &rho);

/**
 * off-design simulation for PCHE
 */
double EXPORT_MY_CODE PCHE_OFFD_Simulation(const char * name, const char * media_hot, const char * media_cold, PCHE_GEO_PARAM * geo, KIM_CORR_COE * cor, SIM_PARAM * sim_param, BoundaryCondtion * bc, PCHECImplResult * retOutput);

double EXPORT_MY_CODE PCHE_OFFD_Simulation_UQ_out(const char * name, const char * media_hot, const char * media_cold, PCHE_GEO_PARAM * geo, KIM_CORR_COE * cor, SIM_PARAM * sim_param, BoundaryCondtion * bc, PCHECImplResult * retOutput, double * Q, double * U);

/**
 * Create new ThermoState with given (p, T), mdot and medium name
 */
ThermoState * EXPORT_MY_CODE NewThermoState_pT(double p, double T, double mdot, const char * medium);

/******* UTILITY FUNCTIONS ********/

/**
 * convert angle in deg to rad
 */
double EXPORT_MY_CODE from_deg(double deg);

/** 
 * convert preseaure in bar to Pa
 */
double EXPORT_MY_CODE from_bar(double p_bar);

/**
 * Convert temperature in deg C to deg K
 */
double EXPORT_MY_CODE from_degC(double degC);

/**
 * test if x and y are almost same
 */
bool EXPORT_MY_CODE the_same(double x, double y, double eps, double & diff_per);

double EXPORT_MY_CODE material_conductivity(double T, bool extrapolate = false);

double EXPORT_MY_CODE print_path_state(const char * name, const char * media, ThermoState * st, int log_level);

void EXPORT_MY_CODE test_struct_param(SIM_PARAM * sim_param, PCHE_GEO_PARAM * geo, BoundaryCondtion * bc, double * h_hot, double * h_cold, double * p_hot, double * p_cold, size_t N_seg);

void EXPORT_MY_CODE setState_C_impl(double p, double M,  State *state);


/** 
 * CoolProp query function
 * Using modelica's external object and low level API of CoolProp to accelerate query speed.
 */
void * EXPORT_MY_CODE init_cp_wrapper(const char * medium, const char * name);

double EXPORT_MY_CODE cp_query(void * wrapper, const char * input_pair,  double val1, double val2, const char * output_name);

void EXPORT_MY_CODE close_cp_wrapper(void * state);

/** 
 * functions utilizing external objects - not very useful since I have to return data

struct PCHE_SIM_EXT_OBJ
{
    int handle_pche;

    ThermoState st_hot_in;
    ThermoState st_cold_in;
    ThermoState st_hot_out;
    ThermoState st_cold_out;
};

// test for transferring c struct as input/output parameter
void * EXPORT_MY_CODE init_PCHE_sim_ext_object(PCHE_GEO_PARAM * geo);

void EXPORT_MY_CODE PCHE_simulate(SIM_PARAM * sim_param, BoundaryCondtion * bc, void * PCHE_ext_obj);

void EXPORT_MY_CODE close_PCHE_sim_ext_object(void * ext_obj);

*/


#ifdef __cplusplus
}
#endif

// template function needs to be implemented in header file. 
template<typename ... Args>
std::string EXPORT_MY_CODE string_format( const std::string& format, Args ... args )
{
    try
    {
        size_t size = std::snprintf( nullptr, 0, format.c_str(), args ... ) + 1; // Extra space for '\0'
        if( size <= 0 ){ throw std::runtime_error( "Error during formatting." ); }
        std::unique_ptr<char[]> buf( new char[ size ] ); 
        std::snprintf( buf.get(), size, format.c_str(), args ... );
        return std::string( buf.get(), buf.get() + size - 1 ); // We don't want the '\0' inside
    }
    catch(const std::exception& e)
    {
        std::cerr << e.what() << '\n';
        return format;
    }   

}

#endif