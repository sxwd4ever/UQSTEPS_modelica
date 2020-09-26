#include "CoolPropLib.h"
#include "AbstractState.h"
#include "DataStructures.h"
#include "MyPropsLib.h"

// #include <string>
#include "crossplatform_shared_ptr.h"
#include "example_dll.h"

using namespace CoolProp;

void EXPORT_MY_CODE MyCall()
{
	hello("World");
}

void EXPORT_MY_CODE MyPropsSI_pT(double p, double T, const std::string &FluidName, double &h, double &rhomass)
{

	/*
	T = CP.PropsSI("T", "P", p, "H", h, PBMedia.mediumName);
    
    mu = CP.PropsSI("V", "P", p, "T", T, PBMedia.mediumName); 

    k = CP.PropsSI("L", "P", p, "T", T, PBMedia.mediumName);   
    
    rho = CP.PropsSI("D", "P", p, "T", T, PBMedia.mediumName); 
	*/
	/*
	std::string str("HEOS");
	std::string fluidName("CO2");
	std::vector<std::string> fluid_names(1);
	fluid_names[0] = "CO2";
	
    AbstractState * pState = AbstractState::factory(str, fluid_names);

	shared_ptr<AbstractState> medium(pState);
    medium->update(PT_INPUTS, p, T); // SI units
	h = medium->hmass();
	rhomass = medium->rhomass();
	*/

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
	/*
    std::cout << "T: " << medium->T() << " K" << std::endl;
    std::cout << "rho': " << medium->rhomass() << " kg/m^3" << std::endl;
    std::cout << "rho': " << medium->rhomolar() << " mol/m^3" << std::endl;
    std::cout << "h': " << medium->hmass() << " J/kg" << std::endl;
    std::cout << "h': " << medium->hmolar() << " J/mol" << std::endl;
    std::cout << "s': " << medium->smass() << " J/kg/K" << std::endl;
    std::cout << "s': " << medium->smolar() << " J/mol/K" << std::endl;
	*/
}

double EXPORT_MY_CODE MyPropsSI_pH(double p, double H, const std::string &FluidName, double &mu, double &k, double &rho)
{
	/*
	shared_ptr<AbstractState> medium(AbstractState::factory("HEOS", FluidName));
    medium->update(HmassP_INPUTS, H, p); // SI units
	double T = medium->T();
	mu = medium->viscosity();
	k = medium->conductivity();
	rho = medium->rhomass();	
	return T;
	*/
	mu = p + 2;
	k = p + H + 1;
	rho = H + 1;
}