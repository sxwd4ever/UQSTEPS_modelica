#include "CoolProp.h"
#include "PCHE.h"
#include <iostream>
#include <math.h>

/**
 * convert angle from degree to rad
 */
double from_deg(double deg)
{
    return deg * M_PI / 180;
}

/** 
 * convert preseaure in bar to Pa
 */
double from_bar(double p_bar)
{
    return p_bar * 1e5;
}

double from_degC(double degC)
{
    return degC + 273.15;
}

using namespace CoolProp;

PCHE_CELL::PCHE_CELL(double p, double T, const char * medium)
{    
    this->p = p;
    this->T = T;
    this->h = PropsSI("H", "P", p, "T", T, medium);
    this->dp = 0;
    this->Q = 0;
    this->G = 0;
    this->h = 0;
    this->u = 0;
    this->mu = 0;
    this->k = 0;
    this->Re = 0;
    this->rho = 0;
    this->Nu = 0;
    this->hc = 0;

    return;

}

PCHE_CELL::~PCHE_CELL()
{
}

PCHE::PCHE(/* args */)
{
}

PCHE::~PCHE()
{

}


int main()
{	
    // // *** geometry parameter ***    
    // // pitch length
    // double pitch = 12e-3;
    // // angle
    // double phi = from_deg((180-108) / 2);
    // // Kim's correlation factors
    // double a = 0, b = 0, c = 0, d = 0;
    // // length of each segment
    // double length_cell = 12e-3;
    // // Diameter of semi_circular
    // double d_c = 2e-3;
    // // Area of semi-circular tube - area per channel
    // double A_c = M_PI * d_c * d_c / 8;
    // // perimeter of semi-circular
    // double peri_c = d_c * M_PI / 2 + d_c ;  
    // // Hydraulic Diameter
    // double d_h = 4 * A_c / peri_c;
    // // thickness of wall between two neighboring hot and cold
    // double t_wall = ( 2 - M_PI_4) * (d_c / 2);
    // // number of channels
    // int N_ch = 80e3;
    // // number of segments
    // int N_seg = 100;
    // // surface area of all cells in a stack 
    // double A_stack =  peri_c * length_cell * N_ch;
    // // Flow area for all channels
    // double A_flow = N_ch * A_c;
    // // length of one pipe in HeatExchanger unit m
    // double length_ch = length_cell * N_seg;

    // // *** boundary condition ***
    // // temperature and pressure
    // double T_cold_in = 500, p_cold_in = from_bar(225);
    // double T_cold_out = 639.15, p_cold_out = from_bar(225);
    // double T_hot_in = 730, p_hot_in = from_bar(90);
    // double T_hot_out = 576.69, p_hot_out = from_bar(90);
    // // mass flow rate of hot/cold side
    // double mdot_hot = 10, mdot_cold = 10;
    // double G_hot = mdot_hot / N_ch / A_c;
    // double G_cold = mdot_cold / N_ch / A_c;

    // // temperature difference between the hot inlet and cold inlet
    // double dt = T_hot_in - T_cold_in; // maximum possible difference

    // int N_iter = 1e4; // maximum iteration number
    // double err = 1e-3; // error for 'identical' test

    // PCHE_CELL cell_hot[N_seg], cell_cold[N_seg];

    // for (size_t i = 0; i < N_seg; i++)
    // {
    //     init_cell(cell_hot[i]);

    //     init_cell(cell_cold[i]);
    // }

    // // suppose (p, T) for hot/cold inlet are known
    // // try to find out (p, T) for cold/cold outlet by trial iteration

    // cell_hot[0].T = T_hot_in;
    // cell_hot[0].p = p_hot_in;
    // // start from the lowest possible temperature
    // cell_cold[0].T = T_cold_in; 
    // cell_cold[0].p = p_cold_in;

    // std::string mname_hot = "CO2", mname_cold = "CO2";
    
    // double h_hot, h_cold; // specific enthalpy for current hot/cold cell

    // for (size_t i = 0; i < N_iter; i++)
    // {
    //     h_hot = PropsSI("H", "P", cell_hot[0].p, "T", cell_hot[0].T, mname_hot);
    //     h_cold = PropsSI("H", "P", cell_cold[0].p, "T", cell_cold[0].T, mname_cold);



    //     if(abs(cell_cold[N_seg].T - T_cold_in) < err )
    //         break;
    // }
    


    std::cout << "Hello world!" << std::endl;

    //hello("dll world");
    
// 	double p = 101325, T = 273.15 + 500;
// 	double h = 0, rhoh = 1, k = 0, mu = 0, rhomass = 0;

// #ifdef WITH_SHARED_LIB_WRAPPER

//     //MyCall();	
// 	MyPropsSI_pT(p, T, "CO2", h, rhoh);
   
//     std::cout << "rhoh': " << rhoh << " mol/m^3" << std::endl;
//     std::cout << "h': " << h << " J/kg" << std::endl;
	
// 	// void MyPropsSI_pH(double p, double H, const std::string &FluidName, double &T, double &mu, double &k, double &rho);
	
// 	T = MyPropsSI_pH(p, h, "CO2",  mu, k, rhoh);

//     std::cout << "T': " << T << " mol/m^3" << std::endl;
//     std::cout << "mu': " << mu << " J/kg" << std::endl;
//     std::cout << "k': " << k << " mol/m^3" << std::endl;
//     std::cout << "rhoh: " << rhoh << " J/kg" << std::endl;	
// #endif

// #ifdef WITH_SHARED_LIB_DIRECT
//     MyCall();	
//     //MyCall("DLL");	
// 	const long buffersize = 500;
//     long errcode = 0;
//     char buffer[buffersize];
//     long handle = _AbstractState_factory("HEOS","CO2", &errcode, buffer, buffersize);
//     long _PT = get_input_pair_index("PT_INPUTS");
//     long _Hmass = get_param_index("HMASS");
// 	long _Dmass = get_param_index("DMASS");
// 	AbstractState_update(handle, _PT, p, T, &errcode, buffer, buffersize);
// 	h = AbstractState_keyed_output(handle, _Hmass, &errcode, buffer, buffersize);
// 	rhomass = AbstractState_keyed_output(handle, _Dmass, &errcode, buffer, buffersize);
//     std::cout << "T: " << h << " K" << std::endl;
//     std::cout << "rho': " << rhomass << " kg/m^3" << std::endl;
// #endif

// #ifdef WITH_STATIC_LIB

//     shared_ptr<AbstractState> Water(AbstractState::factory("HEOS", "Water"));

//     Water->update(PQ_INPUTS, 101325, 0); // SI units
//     std::cout << "T: " << Water->T() << " K" << std::endl;
//     std::cout << "rho': " << Water->rhomass() << " kg/m^3" << std::endl;
//     std::cout << "rho': " << Water->rhomolar() << " mol/m^3" << std::endl;
//     std::cout << "h': " << Water->hmass() << " J/kg" << std::endl;
//     std::cout << "h': " << Water->hmolar() << " J/mol" << std::endl;
//     std::cout << "s': " << Water->smass() << " J/kg/K" << std::endl;
//     std::cout << "s': " << Water->smolar() << " J/mol/K" << std::endl;

// #endif
    return 1;
}