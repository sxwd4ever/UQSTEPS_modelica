#include "CoolPropLib.h"
#include "AbstractState.h"
#include "DataStructures.h"
#include "MyPropsLib.h"
#include "PCHE.h"
#include <math.h>

// #include <string>
#include "crossplatform_shared_ptr.h"
// #include "example_dll.h"

using namespace CoolProp;

/**
 * convert angle from degree to rad
 */
double EXPORT_MY_CODE from_deg(double deg)
{
    return deg * M_PI / 180;
}

/** 
 * convert preseaure in bar to Pa
 */
double EXPORT_MY_CODE from_bar(double p_bar)
{
    return p_bar * 1e5;
}

double EXPORT_MY_CODE from_degC(double degC)
{
    return degC + 273.15;
}

bool EXPORT_MY_CODE the_same(double x, double y, double eps, double & diff_per)
{
    double err = 1e-4;
    if (abs(x - 0) < err)
    {
        diff_per = 1.0;
        return abs(y) <= err;
    }
        
    diff_per = abs((y-x) * 1.0 / x);
    return  diff_per<= eps;
}

// void EXPORT_MY_CODE MyCall()
// {
// 	hello("World");
// }

void EXPORT_MY_CODE MyPropsSI_pT(double p, double T, const std::string &FluidName, double &h, double &rhomass)
{

	// shared_ptr<AbstractState> medium(AbstractState::factory("HEOS", FluidName));
    // medium->update(PT_INPUTS, p, T); // SI units
	// h = medium->hmass();
	// rhomass = medium->rhomass();

	// Call functions declared in CoolPropLib.h
    const long buffersize = 500;
    long errcode = 0;
    char buffer[buffersize];
    long handle = AbstractState_factory("HEOS","CO2", &errcode, buffer, buffersize);
    long _PT = get_input_pair_index("PT_INPUTS");
    long _Hmass = get_param_index("HMASS");
	long _Dmass = get_param_index("DMASS");
	AbstractState_update(handle, _PT, p, T, &errcode, buffer, buffersize);
	h = AbstractState_keyed_output(handle, _Hmass, &errcode, buffer, buffersize);
	rhomass = AbstractState_keyed_output(handle, _Dmass, &errcode, buffer, buffersize);
	
    return;
}

double EXPORT_MY_CODE MyPropsSI_pH(double p, double H, const char * FluidName , double &mu, double &k, double &rho)
{
    const long buffersize = 500;
    long errcode = 0;
    char buffer[buffersize];
    long handle = AbstractState_factory("HEOS",FluidName, &errcode, buffer, buffersize);
    long _HP = get_input_pair_index("HmassP_INPUTS");
	long _T = get_param_index("T");
	long _Dmass = get_param_index("DMASS");
	long _MU = get_param_index("VISCOSITY");
	long _K = get_param_index("CONDUCTIVITY");
	double T = 0.0;

	AbstractState_update(handle, _HP, H , p, &errcode, buffer, buffersize);

	T = AbstractState_keyed_output(handle, _T, &errcode, buffer, buffersize);
	mu = AbstractState_keyed_output(handle, _MU, &errcode, buffer, buffersize);
	k = AbstractState_keyed_output(handle, _K, &errcode, buffer, buffersize);
	rho = AbstractState_keyed_output(handle, _Dmass, &errcode, buffer, buffersize);

	// mu = p + 2;
	// k = p + H + 1;
	// rho = H + 1;

	return T;
}

ThermoState * EXPORT_MY_CODE NewThermoState_pT(double p, double T, double mdot, std::string medium)
{
    ThermoState * st = new ThermoState();
    st->p = p;
    st->T = T;
    st->h = PropsSI("H", "P", p, "T", T, medium.c_str());
    st->mdot = mdot;
    st->medium = medium;
    return st;
}

/**
 * off-design simulation for PCHE
 */
double EXPORT_MY_CODE PCHE_OFFD_Simulation(PCHE_GEO_PARAM * geo, KIM_CORR_COE * cor, SIM_PARAM * sim_param, ThermoState * st_hot_in, ThermoState * st_cold_in, ThermoState * st_hot_out, ThermoState * st_cold_out)
{
    PCHE * pche = new PCHE(*geo);

    pche->set_kim_corr_coe(*cor);

    // *** boundary condition ***
    BoundaryCondtion bc;   

    bc.st_hot_in = st_hot_in;
    bc.st_cold_in = st_cold_in;
    bc.st_hot_out = st_hot_out;
    bc.st_cold_out = st_cold_out;

    pche->simulate(bc, * sim_param);

    delete pche;

	return 0.0;
}