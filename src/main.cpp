//#define _GLIBCXX_USE_CXX11_ABI 0

// #include "example_dll.h"
#include <string>
#include <iostream>

#ifdef WITH_SHARED_LIB_WRAPPER
#include "MyPropsLib.h"
#endif

#ifdef WITH_SHARED_LIB_DIRECT
#include "CoolPropLib.h"
#include "AbstractState.h"
#include "DataStructures.h"
#include "crossplatform_shared_ptr.h"

//using namespace CoolProp;
#endif

#ifdef WITH_STATIC_LIB
#include "CoolProp.h"
#include "MyProps.h"
#include "AbstractState.h"
#include "DataStructures.h"
#include "crossplatform_shared_ptr.h"

#endif

using namespace CoolProp;

int main()

{	
    //hello("dll world");
    
	double p = 101325, T = 273.15 + 500;
	double h = 0, rhoh = 1, k = 0, mu = 0, rhomass = 0;

#ifdef WITH_SHARED_LIB_WRAPPER

    MyCall();	
	MyPropsSI_pT(p, T, "CO2", h, rhoh);
   
    std::cout << "rhoh': " << rhoh << " mol/m^3" << std::endl;
    std::cout << "h': " << h << " J/kg" << std::endl;
	
	// void MyPropsSI_pH(double p, double H, const std::string &FluidName, double &T, double &mu, double &k, double &rho);
	
	T = MyPropsSI_pH(p, h, "CO2",  mu, k, rhoh);

    std::cout << "T': " << T << " mol/m^3" << std::endl;
    std::cout << "mu': " << mu << " J/kg" << std::endl;
    std::cout << "k': " << k << " mol/m^3" << std::endl;
    std::cout << "rhoh: " << rhoh << " J/kg" << std::endl;	
#endif

#ifdef WITH_SHARED_LIB_DIRECT

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
    std::cout << "T: " << h << " K" << std::endl;
    std::cout << "rho': " << rhomass << " kg/m^3" << std::endl;

#endif

#ifdef WITH_STATIC_LIB

    shared_ptr<AbstractState> Water(AbstractState::factory("HEOS", "Water"));

    Water->update(PQ_INPUTS, 101325, 0); // SI units
    std::cout << "T: " << Water->T() << " K" << std::endl;
    std::cout << "rho': " << Water->rhomass() << " kg/m^3" << std::endl;
    std::cout << "rho': " << Water->rhomolar() << " mol/m^3" << std::endl;
    std::cout << "h': " << Water->hmass() << " J/kg" << std::endl;
    std::cout << "h': " << Water->hmolar() << " J/mol" << std::endl;
    std::cout << "s': " << Water->smass() << " J/kg/K" << std::endl;
    std::cout << "s': " << Water->smolar() << " J/mol/K" << std::endl;

#endif
    return 1;
}