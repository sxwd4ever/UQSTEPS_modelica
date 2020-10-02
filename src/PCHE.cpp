#include "CoolPropLib.h"
#include "PCHE.h"
#include <iostream>
#include <math.h>

// using namespace CoolProp;
using namespace std;

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

bool the_same(double x, double y, double eps, double & diff_per)
{
    double err = 1e-4;
    if (abs(x - 0) < err)
    {
        diff_per = 0.0;
        return abs(y) <= err;
    }
        
    diff_per = abs((y-x) * 1.0 / x);
    return  diff_per<= eps;
}

ThermoState::ThermoState(double p, double T, double mdot, std::string  medium_name)
{
    this->p = p;
    this->T = T;
    this->mdot = mdot;    
    this->medium = medium_name;

    this->h = PropsSI("H", "P", p, "T", T, medium_name.c_str());
}

ThermoState::~ThermoState()
{

}


BoundaryCondtion::BoundaryCondtion(/* args */)
{
    // mass flow rate of hot/cold side
    double mdot_hot = 10, mdot_cold = 10;
    std::string mname_hot = "CO2", mname_cold = "CO2";

    st_hot_in = new ThermoState(from_bar(90), 730, mdot_hot, mname_hot);

    st_cold_in = new ThermoState(from_bar(225), 500, mdot_cold, mname_cold);

    st_hot_out = new ThermoState(from_bar(90), 576.69, mdot_hot, mname_hot);

    st_cold_out = new ThermoState(from_bar(225), 639.15, mdot_cold, mname_cold);

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

void PCHE_CELL::init(ThermoState & st, PCHE * pche, std::string medium)
{
    this->_pche = pche;
    this->p = st.p;
    this->T = st.T;
    // this->h = PropsSI("H", "P", p, "T", T, medium.c_str());
    this->h = 0.0;
    this->dp = 0;
    this->Q = 0;
    this->G = 0;
    this->u = 0;
    this->mu = 0;
    this->k = 0;
    this->Re = 0;
    this->rho = 0;
    this->Nu = 0;
    this->hc = 0;

    return;
}

void PCHE_CELL::clone(PCHE_CELL * src)
{
    this->_pche = src->_pche;
    this->p = src->p;
    this->T = src->T;
    this->h = src->h;
    this->dp = src->dp;
    this->Q = src->Q;
    this->G = src->G;
    this->u = src->u;
    this->mu = src->mu;
    this->k = src->k;
    this->Re = src->Re;
    this->rho = src->rho;
    this->Nu = src->Nu;
    this->hc = src->hc;
}

PCHE::PCHE(/* args */)
{
    pitch = 12.3e-3;
    phi = from_deg((180-108)/2);
    // set the correlation coefficients mannually
    a = 0.37934;
    b = 0.82413;
    c = 0.03845;
    d = 0.73793;
    length_cell = 12e-3;
    d_c = 1.5e-3;
    N_ch = 80e3;
    N_seg = 200; // [Kwon2019]'s maximum node number

    this->init();

    // for coolprop query
    _buffer_size = 500;
    _err_code = 0;
    _cp_err_buf = new char[_buffer_size];

    _handle_HP_INPUT = get_input_pair_index("HmassP_INPUTS");

    _handle_T = get_param_index("T");
    _handle_Dmass = get_param_index("DMASS");
    _handle_MU = get_param_index("VISCOSITY");
    _handle_K = get_param_index("CONDUCTIVITY");       
}

PCHE::~PCHE()
{
    delete _cp_err_buf;

    if(_handle_cp_hot != 0)
        AbstractState_free(_handle_cp_hot, & _err_code, _cp_err_buf, _buffer_size);

    if(_handle_cp_cold != 0)    
        AbstractState_free(_handle_cp_cold, &_err_code, _cp_err_buf, _buffer_size);
}

void PCHE::init()
{
    A_c = M_PI * d_c * d_c / 8;
    // perimeter of semi-circular
    peri_c = d_c * M_PI / 2 + d_c ;  
    // Hydraulic Diameter
    // - 1.222 mm for d_c = 2e-3
    // - 0.922 mm for d_c = 1.5e-3
    // or for semi-circular pipe, we have 
    // d_c = ( 1 + 2 / pi) * d_h
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
    if (count == 0)
        return cell_seq[0].T;

    double sum = 0;   

    for (size_t i = 0; i < count; i++)
        sum += cell_seq[i].T;

    return sum / count;    
}

bool PCHE::calc_U(int idx, double & h_hot, double & h_cold, BoundaryCondtion & bc)
{
    bool heat_exchanged = true;

    PCHE_CELL * c_h = &cell_hot[idx];
    PCHE_CELL * c_c = &cell_cold[idx];

    AbstractState_update(_handle_cp_hot, _handle_HP_INPUT, h_hot, c_h->p, &_err_code, _cp_err_buf, _buffer_size);

    AbstractState_update(_handle_cp_cold, _handle_HP_INPUT, h_cold, c_c->p, &_err_code, _cp_err_buf, _buffer_size);

	// T = AbstractState_keyed_output(handle, _T, &errcode, buffer, buffersize);
	// mu = AbstractState_keyed_output(handle, _MU, &errcode, buffer, buffersize);
	// k = AbstractState_keyed_output(handle, _K, &errcode, buffer, buffersize);
	// rho = AbstractState_keyed_output(handle, _Dmass, &errcode, buffer, buffersize);     

    // double p_cold = bc.st_cold_in->p;
    // double p_hot = bc.st_hot_in->p;

    // const char * medium_hot = bc.st_hot_in->medium.c_str();
    // const char * medium_cold = bc.st_cold_in->medium.c_str();

    c_h->T = AbstractState_keyed_output(_handle_cp_hot, _handle_T, &_err_code, _cp_err_buf, _buffer_size);
    c_c->T = AbstractState_keyed_output(_handle_cp_cold, _handle_T, &_err_code, _cp_err_buf, _buffer_size);

    c_h->mu = AbstractState_keyed_output(_handle_cp_hot, _handle_MU, &_err_code, _cp_err_buf, _buffer_size);
    c_c->mu = AbstractState_keyed_output(_handle_cp_cold, _handle_MU, &_err_code, _cp_err_buf, _buffer_size);

    c_h->k = AbstractState_keyed_output(_handle_cp_hot, _handle_K, &_err_code, _cp_err_buf, _buffer_size);
    c_c->k = AbstractState_keyed_output(_handle_cp_cold, _handle_K, &_err_code, _cp_err_buf, _buffer_size);

    c_h->rho = AbstractState_keyed_output(_handle_cp_hot, _handle_Dmass, &_err_code, _cp_err_buf, _buffer_size);
    c_c->rho = AbstractState_keyed_output(_handle_cp_cold, _handle_Dmass, &_err_code, _cp_err_buf, _buffer_size);

    c_h->Re = this->G_hot * this->d_h / c_h->mu;
    c_c->Re = this->G_cold * this->d_h / c_c->mu;

    c_h->u = bc.st_hot_in->mdot / this->A_flow / c_h->rho;
    c_c->u = bc.st_cold_in->mdot / this->A_flow / c_c->rho;

    c_h->Nu = 4.089 + c * pow(c_h->Re, d);
    c_c->Nu = 4.089 + c * pow(c_c->Re, d);
    // thermal conductivity of the wall
    double kw = 27; // mock trial value
    c_h->hc = c_h->Nu * c_h->k / this->d_h;
    c_c->hc = c_c->Nu * c_c->k / this->d_h;

    c_h->f = (15.78 + a * pow(c_h->Re, b)) / c_h->Re;
    c_c->f = (15.78 + a * pow(c_h->Re, b)) / c_c->Re;

    U[idx] = 1 /( 1 / c_h->hc + 1 / c_c->hc + this->t_wall / kw);

    if(c_h->T > c_c->T)
        q[idx] = U[idx] * A_stack * (c_h->T - c_c->T);
    else
    {
        q[idx] = 0;
        heat_exchanged = false;
    }

    c_h->dp = c_h->f * length_cell * c_h->rho * pow(c_h->u, 2) / this->d_h;
    c_c->dp = c_c->f * length_cell * c_c->rho * pow(c_c->u, 2) / this->d_h;

    h_hot = (bc.st_hot_in->mdot * h_hot - q[idx]) / bc.st_hot_in->mdot;
    h_cold = (bc.st_cold_in->mdot * h_cold + q[idx] ) / bc.st_cold_in->mdot;

    // set temperature and pressure for next cell
    // if(q[idx] > 0)
    // {
    cell_hot[idx + 1].p = c_h->p - c_h->dp;

    cell_hot[idx + 1].h = h_hot; // PropsSI("T", "P", cell_hot[idx + 1].p, "H", h_hot, medium_hot);

    cell_cold[idx + 1].p = c_c->p + c_c->dp;

    cell_cold[idx + 1].h = h_cold; //PropsSI("T", "P", cell_cold[idx + 1].p, "H", h_cold, medium_cold);
    // }

    return heat_exchanged;
}

void PCHE::simulate(BoundaryCondtion & bc, SIM_PARAM & sim_para)
{
    ThermoState st_hot(bc.st_hot_in->p, 0, bc.st_hot_in->mdot, bc.st_hot_in->medium);

    ThermoState st_cold(bc.st_cold_in->p, 0, bc.st_cold_in->mdot, bc.st_cold_in->medium);

    for (size_t i = 0; i < N_seg; i++)
    {
        this->cell_hot[i].init(st_hot, this);
        this->cell_cold[i].init(st_cold, this);

        this->U[i] = 0;
        this->q[i] = 0;
    }

    _handle_cp_cold = AbstractState_factory("HEOS", st_cold.medium.c_str(), &this->_err_code, _cp_err_buf, _buffer_size);

    _handle_cp_hot = AbstractState_factory("HEOS", st_hot.medium.c_str(), &this->_err_code, _cp_err_buf, _buffer_size);

    this->G_hot = bc.st_hot_in->mdot / N_ch / A_c;
    this->G_cold = bc.st_cold_in->mdot / N_ch / A_c;

    // temperature difference between the hot inlet and cold inlet
    double dT = 0.0;

    // suppose (p, T) for hot/cold inlet are known
    // try to find out (p, T) for cold/cold outlet by trial iteration
    cell_hot[0].T = bc.st_hot_in->T;
    cell_hot[0].p = bc.st_hot_in->p;
    cell_cold[0].T = bc.st_hot_in->T - 5; 
    cell_cold[0].p = bc.st_cold_in->p;
    
    double h_hot, h_cold; // specific enthalpy for current hot/cold cell
    double diff_per = 0.0;
    bool heat_exchanged = true;

    for (size_t i = 0; i < sim_para.N_iter; i++)
    {
        h_hot = PropsSI("H", "P", cell_hot[0].p, "T", cell_hot[0].T, bc.st_hot_in->medium.c_str());
        h_cold = PropsSI("H", "P", cell_cold[0].p, "T", cell_cold[0].T, bc.st_cold_in->medium.c_str());

        heat_exchanged = true;
        diff_per = 0.0;

        for(size_t idx_cell = 0; idx_cell < this->N_seg; idx_cell ++)
            heat_exchanged = this->calc_U(idx_cell, h_hot, h_cold, bc);
            /*
            if(heat_exchanged)
                heat_exchanged = this->calc_U(idx_cell, h_hot, h_cold, bc);
            else if(idx_cell != 0) 
            // skip the rest of the cells if no heat exchanged between cold and hot side
                {
                    this->cell_hot[idx_cell].clone(&this->cell_hot[idx_cell - 1]);
                    this->cell_cold[idx_cell].clone(&this->cell_cold[idx_cell - 1]);
                }
                else
                {
                    cout<<"error initial trial value since no heat exchanged"<<endl;
                }
            */

        if(the_same(cell_cold[N_seg - 1].T, bc.st_cold_in->T, sim_para.err, diff_per) )
            break;

        cout<<i<<" th trial, diff(%)="<< diff_per<<endl;

        // update dt and T_cold_in accordingly
        dT = bc.st_cold_in->T - this->cell_cold[this->N_seg - 1].T;
        this->cell_cold[0].T += 0.5 * dT;
    }
}

int main()
{	
    PCHE * pche = new PCHE();

    // *** boundary condition ***
    BoundaryCondtion bc;   

    SIM_PARAM sim_param;
    sim_param.N_iter = 1e4;
    sim_param.err = 1e-2; // 1 percent

    pche->simulate(bc, sim_param);

    std::cout << "All done" << std::endl;
   
    delete pche;

    return 1;
}