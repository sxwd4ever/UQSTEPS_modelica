#include "CoolProp.h"
#include <iostream>
#include <vector>
#include <string>
#include "AbstractState.h"
#include "DataStructures.h"
#include <memory>
#include "myProps.cpp"


using namespace CoolProp;

int main()
{	
	double p = 101325, T = 273.15 + 500;
	double h = 0, rhoh = 1, k = 0, mu = 0;
	
	MyPropsSI_pT(p, T, "CO2", h, rhoh);
    
    std::cout << "rhoh': " << rhoh << " mol/m^3" << std::endl;
    std::cout << "h': " << h << " J/kg" << std::endl;
	
	// void MyPropsSI_pH(double p, double H, const std::string &FluidName, double &T, double &mu, double &k, double &rho);
	
	MyPropsSI_pH(p, h, "CO2", T, mu, k, rhoh);

    std::cout << "T': " << T << " mol/m^3" << std::endl;
    std::cout << "mu': " << mu << " J/kg" << std::endl;
    std::cout << "k': " << k << " mol/m^3" << std::endl;
    std::cout << "rhoh: " << rhoh << " J/kg" << std::endl;	
	

	/*
    shared_ptr<AbstractState> Water(AbstractState::factory("HEOS", "Water"));
    Water->update(PQ_INPUTS, 101325, 0); // SI units
    std::cout << "T: " << Water->T() << " K" << std::endl;
    std::cout << "rho': " << Water->rhomass() << " kg/m^3" << std::endl;
    std::cout << "rho': " << Water->rhomolar() << " mol/m^3" << std::endl;
    std::cout << "h': " << Water->hmass() << " J/kg" << std::endl;
    std::cout << "h': " << Water->hmolar() << " J/mol" << std::endl;
    std::cout << "s': " << Water->smass() << " J/kg/K" << std::endl;
    std::cout << "s': " << Water->smolar() << " J/mol/K" << std::endl;
	*/
    return 1;
}