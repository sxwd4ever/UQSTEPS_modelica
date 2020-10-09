#include "CoolPropLib.h"
#include "PCHE.h"
#include <iostream>
#include <math.h>
#include "MyPropsLib.h"

// using namespace CoolProp;
using namespace std;

PCHE_CELL::PCHE_CELL()
{    
    _pche = NULL;
    this->p = 0;
    this->T = 0;
    this->h = 0;
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
}

PCHE_CELL::~PCHE_CELL()
{
}

void PCHE_CELL::_init(ThermoState * st, PCHE * pche)
{
    this->_pche = pche;
    this->p = st->p;
    this->T = st->T;
    // this->h = PropsSI("H", "P", p, "T", T, medium.c_str());
    this->h = st->h;
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

PCHE::PCHE(PCHE_GEO_PARAM & geo)
{
    _cell_hot = NULL;
    _cell_cold = NULL;

    this->_init(geo);

    // for coolprop query
    _buffer_size = 500;
    _err_code = 0;
    _cp_err_buf = new char[_buffer_size];

    _handle_HP_INPUT = get_input_pair_index("HmassP_INPUTS");
    _handle_PT_INPUT = get_input_pair_index("PT_INPUTS");

    _handle_T = get_param_index("T");
    _handle_H = get_param_index("HMASS");
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

    delete [] _U;
    delete [] _Q;

    try
    {
        delete [] _cell_hot;
        delete [] _cell_cold;
    }
    catch(const std::exception& e)
    {
        std::cerr << e.what() << '\n';
    }    
}

void PCHE::set_kim_corr_coe(KIM_CORR_COE & coe)
{
    _corr_coe.a = coe.a;
    _corr_coe.b = coe.b;
    _corr_coe.c = coe.c;
    _corr_coe.d = coe.d;
}

void PCHE::_init(PCHE_GEO_PARAM & geo)
{
    _geo.pitch = geo.pitch;
    _geo.phi = geo.phi;
    _geo.length_cell = geo.length_cell;
    _geo.d_c = geo.d_c;
    _geo.N_ch = geo.N_ch;
    _geo.N_seg = geo.N_seg;

    _U = new double[_geo.N_seg];
    _Q = new double[_geo.N_seg];

    _cell_hot = new PCHE_CELL[_geo.N_seg]; 
    _cell_cold = new PCHE_CELL[_geo.N_seg];

    _A_c = M_PI * _geo.d_c * _geo.d_c / 8;
    _peri_c = _geo.d_c * M_PI / 2 + _geo.d_c ;  
    // Hydraulic Diameter
    // - 1.222 mm for d_c = 2e-3
    // - 0.922 mm for d_c = 1.5e-3
    // or for semi-circular pipe, we have 
    // d_c = ( 1 + 2 / pi) * _d_h
    _d_h = 4 * _A_c / _peri_c;
    _t_wall = ( 2 - M_PI_4) * (_geo.d_c / 2);    
    _A_stack =  _peri_c * _geo.length_cell * _geo.N_ch;    
    _A_flow = _geo.N_ch * _A_c;    
    _length_ch = _geo.length_cell * _geo.N_seg;
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

    PCHE_CELL * c_h = &_cell_hot[idx];
    PCHE_CELL * c_c = &_cell_cold[idx];

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

    c_h->T = _cp_props(_handle_cp_hot, _handle_T);
    c_c->T = _cp_props(_handle_cp_cold, _handle_T);

    c_h->mu = _cp_props(_handle_cp_hot, _handle_MU); 
    c_c->mu = _cp_props(_handle_cp_cold, _handle_MU); 

    c_h->k = _cp_props(_handle_cp_hot, _handle_K);  
    c_c->k = _cp_props(_handle_cp_cold, _handle_K);   

    c_h->rho = _cp_props(_handle_cp_hot, _handle_Dmass);  
    c_c->rho = _cp_props(_handle_cp_cold, _handle_Dmass); 

    c_h->Re = this->_G_hot * this->_d_h / c_h->mu;
    c_c->Re = this->_G_cold * this->_d_h / c_c->mu;

    c_h->u = bc.st_hot_in.mdot / _A_flow / c_h->rho;
    c_c->u = bc.st_cold_in.mdot / _A_flow / c_c->rho;

    c_h->Nu = 4.089 + _corr_coe.c * pow(c_h->Re, _corr_coe.d);
    c_c->Nu = 4.089 + _corr_coe.c * pow(c_c->Re, _corr_coe.d);
    // thermal conductivity of the wall
    double kw = 27; // mock trial value
    c_h->hc = c_h->Nu * c_h->k / this->_d_h;
    c_c->hc = c_c->Nu * c_c->k / this->_d_h;

    c_h->f = (15.78 + _corr_coe.a * pow(c_h->Re, _corr_coe.b)) / c_h->Re;
    c_c->f = (15.78 + _corr_coe.a * pow(c_h->Re, _corr_coe.b)) / c_c->Re;

    _U[idx] = 1 /( 1 / c_h->hc + 1 / c_c->hc + this->_t_wall / kw);

    if(c_h->T > c_c->T)
        _Q[idx] = _U[idx] * _A_stack * (c_h->T - c_c->T);
    else
    {
        _Q[idx] = 0;
        heat_exchanged = false;
        if(_idx_pinch < 0) // record the pinch point
            _idx_pinch = idx;
    }

    c_h->dp = c_h->f * _geo.length_cell * c_h->rho * pow(c_h->u, 2) / this->_d_h;
    c_c->dp = c_c->f * _geo.length_cell * c_c->rho * pow(c_c->u, 2) / this->_d_h;

    h_hot = (bc.st_hot_in.mdot * h_hot - _Q[idx]) / bc.st_hot_in.mdot;
    h_cold = (bc.st_cold_in.mdot * h_cold - _Q[idx] ) / bc.st_cold_in.mdot;

    // set temperature and pressure for next cell
    // if(_q[idx] > 0)
    // {
    _cell_hot[idx + 1].p = c_h->p - c_h->dp;

    _cell_hot[idx + 1].h = h_hot; // PropsSI("T", "P", _cell_hot[idx + 1].p, "H", h_hot, medium_hot);

    _cell_cold[idx + 1].p = c_c->p + c_c->dp;

    _cell_cold[idx + 1].h = h_cold; //PropsSI("T", "P", _cell_cold[idx + 1].p, "H", h_cold, medium_cold);
    // }

    return heat_exchanged;
}

double PCHE::_cp_props(long handle_stream, long handle_prop, double def_value)
{
    _err_code = 0;

    double prop = AbstractState_keyed_output(handle_stream, handle_prop, &_err_code, _cp_err_buf, _buffer_size);

    if(_err_code != 0)
        // handle the err
        prop = def_value;

    return prop;
}

void PCHE::_print_boundary_conditions(BoundaryCondtion & bc)
{
    cout<< "**** boundary condition ****"<<endl;
    cout<< "(p, T, h)_hot_in = ("<<bc.st_hot_in.p<<", "<<bc.st_hot_in.T<<", "<<bc.st_hot_in.h<<");"<< endl;
    cout<< "(p, T, h)_cold_in = ("<<bc.st_cold_in.p<<", "<<bc.st_cold_in.T<<", "<<bc.st_cold_in.h<<");"<< endl;
    cout<< "(p, T, h)_hot_out = ("<<bc.st_hot_out.p<<", "<<bc.st_hot_out.T<<", "<<bc.st_hot_out.h<<");"<< endl;
    cout<< "(p, T, h)_cold_out = ("<<bc.st_cold_out.p<<", "<<bc.st_cold_out.T<<", "<<bc.st_cold_out.h<<");"<< endl;
}

void PCHE::_print_geo_param(PCHE_GEO_PARAM & geo)
{
    cout<< "**** PCHE's Geometry params ****"<<endl;
    cout<< "pitch       = "<<geo.pitch<<endl;
    cout<< "phi         = "<<geo.phi<<endl;
    cout<< "length_cell = "<<geo.length_cell<<endl;
    cout<< "d_c         = "<<geo.d_c<<endl;
    cout<< "N_ch        = "<<geo.N_ch<<endl;
    cout<< "N_seg       = "<<geo.N_seg<<endl;
}

void PCHE::_print_sim_param(SIM_PARAM & sim_param)
{
    cout<< "**** simulation params ****"<<endl;
    cout<< "err       = "<<sim_param.err<<endl;
    cout<< "dT_init   = "<<sim_param.delta_T_init<<endl;
    cout<< "N_iter    = "<<sim_param.N_iter<<endl;
    cout<< "step_rel  = "<<sim_param.step_rel<<endl;
}

void PCHE::simulate(BoundaryCondtion & bc, SIM_PARAM & sim_para, SimulationResult & sr)
{
    _print_geo_param(_geo);

    _print_boundary_conditions(bc);

    _print_sim_param(sim_para);
    
    cout<<"ready to start simulation"<<endl;

    ThermoState * st_hot = NewThermoState_pT(bc.st_hot_in.p, 0, bc.st_hot_in.mdot, bc.media_hot);
    ThermoState * st_cold = NewThermoState_pT(bc.st_cold_in.p, 0, bc.st_cold_in.mdot, bc.media_cold);

    for (size_t i = 0; i < _geo.N_seg; i++)
    {
        this->_cell_hot[i]._init(st_hot, this);
        this->_cell_cold[i]._init(st_cold, this);

        this->_U[i] = 0;
        this->_Q[i] = 0;
    }

    _handle_cp_hot = AbstractState_factory("HEOS", bc.media_hot, &this->_err_code, _cp_err_buf, _buffer_size);
    _handle_cp_cold = AbstractState_factory("HEOS", bc.media_cold, &this->_err_code, _cp_err_buf, _buffer_size);

    this->_G_hot = bc.st_hot_in.mdot / _geo.N_ch / _A_c;
    this->_G_cold = bc.st_cold_in.mdot / _geo.N_ch / _A_c;

    // temperature difference between the hot inlet and cold inlet
    double dT = 0.0;

    // suppose (p, T) for hot/cold inlet are known
    // try to find out (p, T) for cold/cold outlet by trial iteration
    _cell_hot[0].T = bc.st_hot_in.T;
    _cell_hot[0].p = bc.st_hot_in.p;

    _cell_cold[0].T = bc.st_hot_in.T - sim_para.delta_T_init; // start from t_hot_in - 5
    _cell_cold[0].p = bc.st_cold_in.p;

    //dT = (_cell_cold[0].T - bc.st_cold_in->T) * sim_para.step_rel;
    
    double h_hot, h_cold; // specific enthalpy for current hot/cold cell
    double diff_per = 0.0;
    bool heat_exchanged = true;

    for (size_t i = 0; i < sim_para.N_iter; i++)
    {
        heat_exchanged = true;
        _idx_pinch = -1;
        diff_per = 0.0;

        AbstractState_update(_handle_cp_hot, _handle_PT_INPUT, _cell_hot[0].p, _cell_hot[0].T, &_err_code, _cp_err_buf, _buffer_size);
        AbstractState_update(_handle_cp_cold, _handle_PT_INPUT, _cell_cold[0].p, _cell_cold[0].T, &_err_code, _cp_err_buf, _buffer_size);    

        h_hot = _cp_props(_handle_cp_hot, _handle_H);
        h_cold = _cp_props(_handle_cp_cold, _handle_H);

        for(int idx_cell = 0; idx_cell < _geo.N_seg; idx_cell ++)
            // heat_exchanged = this->calc_U(idx_cell, h_hot, h_cold, bc);
            if(heat_exchanged)
                heat_exchanged = this->calc_U(idx_cell, h_hot, h_cold, bc);
            else 
            {
                _cell_hot[idx_cell].clone(&_cell_hot[idx_cell - 1]);
                _cell_cold[idx_cell].clone(&_cell_cold[idx_cell - 1]);
                //skip the rest of the cells if no heat exchanged between cold and hot side
                // T of rest cells remain 0 and the diff in T between bc.T_cold_in and last
                // cold cell's T will drive next iteration
                //break;
            }
        // calculate delta T for next iteration
        //dT = bc.st_cold_in->T - _cell_cold[_idx_pinch].T;

        dT = bc.st_cold_in.T - _cell_cold_in()->T;
        bool converged = the_same(_cell_cold_in()->T, bc.st_cold_in.T, sim_para.err, diff_per);

        cout<<i<<" th trial, diff(%)="<< diff_per<<"; ";
        cout<<"T_{cold,out}="<<_cell_cold[0].T <<" K; ";
        cout<<"dT="<<dT<<" K; ";
        if(_idx_pinch > 0)
            cout<<"T_{cold, pinch="<<_idx_pinch<<"}="<<_cell_cold[_idx_pinch].T<<" K"<<endl;
        else
            cout<<"T_{cold,in}="<< _cell_cold_in()->T<<" K"<<endl;

        if(converged)
            break;

        // update dt and T_cold_out accordingly, for next iteration
        _cell_cold[0].T += dT * sim_para.step_rel;

        // if(_cell_cold_in()->T > bc.st_cold_in->T)
        //     _cell_cold[0].T -= dT;
        // else
        //     _cell_cold[0].T += dT;
    }

    delete st_cold;
    delete st_hot;

    // assign the simulation result as output

    for (size_t i = 0; i < sr.N_seg; i++)
    {
        sr.p_hot[i] = _cell_hot[i].p;
        sr.p_cold[i] = _cell_cold[i].p;
        sr.h_hot[i] = _cell_hot[i].h;
        sr.h_cold[i] = _cell_cold[i].h;
    }

    return;
}

PCHE_CELL * PCHE::_cell_cold_in()
{
    return &_cell_cold[_geo.N_seg - 1];
}

