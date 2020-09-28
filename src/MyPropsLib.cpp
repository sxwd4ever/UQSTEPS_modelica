#include "CoolPropLib.h"
#include "AbstractState.h"
#include "DataStructures.h"
#include "MyPropsLib.h"

// #include <string>
#include "crossplatform_shared_ptr.h"
// #include "example_dll.h"

using namespace CoolProp;

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