#include "CoolProp.h"
#include "AbstractState.h"
#include "DataStructures.h"
#include "myProps.h"

#include <string>
#include <memory>

using namespace CoolProp;


void MyPropsSI_pT(double p, double T, const std::string &FluidName, double &h, double &rhomass)
{
	
	/*
	    T = CP.PropsSI("T", "P", p, "H", h, PBMedia.mediumName);
    
    mu = CP.PropsSI("V", "P", p, "T", T, PBMedia.mediumName); 

    k = CP.PropsSI("L", "P", p, "T", T, PBMedia.mediumName);   
    
    rho = CP.PropsSI("D", "P", p, "T", T, PBMedia.mediumName); 
	*/
	
	shared_ptr<AbstractState> medium(AbstractState::factory("HEOS", FluidName));
    medium->update(PT_INPUTS, p, T); // SI units
	h = medium->hmass();
	rhomass = medium->rhomass();
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

void MyPropsSI_pH(double p, double H, const std::string &FluidName, double &T, double &mu, double &k, double &rho)
{
	shared_ptr<AbstractState> medium(AbstractState::factory("HEOS", FluidName));
    medium->update(HmassP_INPUTS, H, p); // SI units
	T = medium->T();
	mu = medium->viscosity();
	k = medium->conductivity();
	rho = medium->rhomass();	
}