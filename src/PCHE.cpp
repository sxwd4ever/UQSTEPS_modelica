#include "CoolProp.h"
#include "PCHE.h"
#include <iostream>
#include <math.h>

using namespace CoolProp;

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

ThermoState::ThermoState(double p, double T, double mdot, const char * medium)
{
    p = p;
    T = T;
    mdot = mdot;
    medium = medium;
}

ThermoState::~ThermoState()
{

}


BoundaryCondtion::BoundaryCondtion(/* args */)
{
    // mass flow rate of hot/cold side
    double mdot_hot = 10, mdot_cold = 10;

    st_hot_in = new ThermoState(730, from_bar(90), mdot_hot, "CO2");

    st_cold_in = new ThermoState(500, from_bar(225), mdot_cold, "CO2");

    st_hot_out = new ThermoState(576.69, from_bar(90), mdot_hot, "CO2");

    st_cold_out = new ThermoState(639.15, from_bar(225), mdot_cold, "CO2");

}

BoundaryCondtion::~BoundaryCondtion()
{
    delete st_hot_in;
    delete st_cold_in;
    delete st_hot_out;
    delete st_cold_out;
}

PCHE_CELL::PCHE_CELL()
{    
    _pche = NULL;
}

PCHE_CELL::~PCHE_CELL()
{
}

void PCHE_CELL::init(ThermoState & st, PCHE * pche, const char * medium)
{
    this->_pche = pche;
    this->p = st.p;
    this->T = st.T;
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

PCHE::PCHE(/* args */)
{
    pitch = 12e-3;
    phi = from_degC((180-108)/2);
    a = 0;
    b = 0;
    c = 0;
    d = 0;
    length_cell = 12e-3;
    d_c = 2e-3;
    N_ch = 80e3;
    N_seg = 100;
}

PCHE::~PCHE()
{

}

void PCHE::init()
{
    A_c = M_PI * d_c * d_c / 8;
    // perimeter of semi-circular
    peri_c = d_c * M_PI / 2 + d_c ;  
    // Hydraulic Diameter
    d_h = 4 * A_c / peri_c;
    // thickness of wall between two neighboring hot and cold
    t_wall = ( 2 - M_PI_4) * (d_c / 2);
    // surface area of all cells in a stack 
    A_stack =  peri_c * length_cell * N_ch;
    // Flow area for all channels
    A_flow = N_ch * A_c;
    // length of one pipe in HeatExchanger unit m
    length_ch = length_cell * N_seg;
}

double PCHE::avg_T(PCHE_CELL * cell_seq, int count)
{
    double sum = 0;
    for (size_t i = 0; i < count; i++)
        sum += cell_seq[i].T;

    return sum / count;    
}

void PCHE::calc_U(int idx, double & h_hot, double & h_cold, BoundaryCondtion & bc)
{
    PCHE_CELL * c_h = &cell_hot[idx];
    PCHE_CELL * c_c = &cell_cold[idx];

    double p_cold = bc.st_cold_in-> p;
    double p_hot = bc.st_hot_in->p;

    char * medium_hot = bc.st_hot_in->medium;
    char * medium_cold = bc.st_cold_in->medium;

    c_h->mu = PropsSI("V", "T", c_h->T, "P", p_hot, medium_hot);
    c_c->mu = PropsSI("V", "T", c_c->T, "P", p_cold, medium_cold);

    c_h->k = PropsSI("L", "T", c_h->T, "P", p_hot, medium_hot);
    c_c->k = PropsSI("L", "T", c_c->T, "P", p_cold, medium_cold);

    c_h->Re = this->G_hot * this->d_h / c_h->mu;
    c_c->Re = this->G_cold * this->d_h / c_c->mu;

    c_h->rho = PropsSI("D", "T", avg_T(this->cell_hot, i), "P", p_hot, medium_hot);
    c_c->rho = PropsSI("D", "T", avg_T(this->cell_cold, i), "P", p_cold, medium_cold);

    c_h->u = c_h.mdot / this->A_flow / c_h->rho;
    c_c->u = c_c.mdot / this->A_flow / c_c->rho;

    c_h->Nu = 4.089 + c * c_h->Re ^ d;
    c_c->Nu = 4.089 + c * c_c->Re ^ d;
    // thermal conductivity of the wall
    double kw = 1.0;
    c_h->hc = c_h->Nu * c_h->k / this->d_h;
    c_c->hc = c_c->Nu * c_c->k / this->d_h;

    U[idx] = 1 /( 1 / c_h->hc + 1 / c_c->hc + this->t_wall / kw);

    c_h->f = (15.78 + a * c_h->Re ^ b) / c_h->Re;
    c_c->f = (15.78 + a * c_h->Re ^ b) / c_c->Re;


    /*
    self.ReH[i]=self.Gh*self.dh/self.muh[i]
    self.ReC[i]=self.Gc*self.dh/self.muc[i]
    self.rhoH[i]=hotfluid.prop("D","T",np.average(self.Th[0:i+1]), "P", self.hot.P.p)
    self.rhoC[i]=coldfluid.prop("D","T",np.average(self.Tc[0:i+1]), "P", self.cold.P.p)
    uH=self.hot.P.mdot/self.flowarea()/self.rhoH[i]
    uC=self.cold.P.mdot/self.flowarea()/self.rhoC[i]
    area=self.peri*dL*self.N  # Surface area of the channel segment, m2
    if method=="Kim2012":
        self.NuH[i]=4.089+self.abcd[2]*(self.ReH[i])**self.abcd[3]
        self.NuC[i]=4.089+self.abcd[2]*(self.ReC[i])**self.abcd[3]
        kw=thconductivity(self.metal, (self.Th[i]+self.Tc[i])/2)
        self.hh[i]=self.NuH[i]*self.kh[i]/self.dh
        self.hc[i]=self.NuC[i]*self.kc[i]/self.dh
        self.U[i]=1/(1/self.hh[i]+1/self.hc[i]+self.tw/kw)
        self.V[3]=kw/self.tw
        self.fH[i]=(15.78+self.abcd[0]*self.ReH[i]**self.abcd[1])/self.ReH[i]
        self.fC[i]=(15.78+self.abcd[0]*self.ReC[i]**self.abcd[1])/self.ReC[i]  
        */  
}

void PCHE::simulate(BoundaryCondtion & bc, SIM_PARAM & sim_para)
{
    for (size_t i = 0; i < N_seg; i++)
    {
        ThermoState st(bc.st_hot_in->p, 0, bc.st_hot_in->mdot, "CO2");
        this->cell_hot[i].init(st, this);

        st = ThermoState(bc.st_cold_in->p, 0, bc.st_cold_in->mdot, "CO2");
        this->cell_cold[i].init(st, this);

        this->U[i] = 0;
    }

    this->G_hot = bc.st_hot_in->mdot / N_ch / A_c;
    this->G_cold = bc.st_cold_in->mdot / N_ch / A_c;

    // temperature difference between the hot inlet and cold inlet
    double dT = bc.st_hot_in->T - bc.st_cold_in->T; // maximum possible difference

    // suppose (p, T) for hot/cold inlet are known
    // try to find out (p, T) for cold/cold outlet by trial iteration

    cell_hot[0].T = bc.st_hot_in->T;
    cell_hot[0].p = bc.st_hot_in->p;
    // start from the lowest possible temperature
    cell_cold[0].T = bc.st_cold_in->T; 
    cell_cold[0].p = bc.st_cold_in->p;

    std::string mname_hot = "CO2", mname_cold = "CO2";
    
    double h_hot, h_cold; // specific enthalpy for current hot/cold cell

    for (size_t i = 0; i < sim_para.N_iter; i++)
    {
        h_hot = PropsSI("H", "P", cell_hot[0].p, "T", cell_hot[0].T, mname_hot);
        h_cold = PropsSI("H", "P", cell_cold[0].p, "T", cell_cold[0].T, mname_cold);

        for(size_t idx_cell = 0; idx_cell < this->N_seg; idx_cell ++)
            this->calc_U(idx_cell, h_hot, h_cold, bc);

        if(abs(cell_cold[N_seg].T - bc.st_cold_in->T) < sim_para.err )
            break;

        // update dt and T_cold_in accordingly
        dT = bc.st_cold_in->T - this->cell_cold[this->N_seg].T;
        this->cell_cold[0].T += 0.2 * dT;
    }
}

int main()
{	
    // // *** geometry parameter ***    

    // *** boundary condition ***

    PCHE pche;

    BoundaryCondtion bc;   

    SIM_PARAM sim_param;
    sim_param.N_iter = 1e4;
    sim_param.err = 1e-3;

    pche.simulate(bc, sim_param);

    std::cout << "All done" << std::endl;

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