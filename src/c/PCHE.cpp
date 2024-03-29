#include "CoolPropLib.h"
#include "PCHE.h"
#include <iostream>
#include <sstream>
#include <exception>
#include <math.h>
#include "MyPropsLib.h"

// using namespace CoolProp;
using namespace std;

static int g_log_level = log_level::OFF;

PCHE_CELL::PCHE_CELL()
{    
    _pche = NULL;
    idx = 0;
    type_stream = 0;
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

void PCHE_CELL::init(int idx, int type, PCHE * pche)
{
    _pche = pche;
    this->idx = idx;
    this->type_stream = type;
}

std::string PCHE_CELL::to_string()
{
    return string_format("(p, T, h)_{%s, %d} = %f, %f, %f\n", type_name(),idx, p, T, h);
}

void PCHE_CELL::print_state()
{
    cout<<to_string();
}

bool PCHE_CELL::validate()
{
    bool valid = true;

    if(p < 7e6 || p > 30e6 )
        valid = false;
    
    if(T < from_degC(0) || T > from_degC(1200))
        valid = false;

    if(h < 0 || h > 1e8)
        valid = false;

    return valid;
}

const char * PCHE_CELL::type_name()
{
    const char * type = "hot";
    if(type_stream == 1)
        type = "cold";
    return type;
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

PCHE::PCHE(const char * name, PCHE_GEO_PARAM & geo)
{
    _name = name;
 
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
    if(_handle_cp_hot != 0)
        AbstractState_free(_handle_cp_hot, & _err_code, _cp_err_buf, _buffer_size);

    if(_handle_cp_cold != 0)    
        AbstractState_free(_handle_cp_cold, &_err_code, _cp_err_buf, _buffer_size);

    delete _cp_err_buf;

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
    _geo.length = geo.length;
    _geo.d_c = geo.d_c;
    _geo.N_ch = geo.N_ch;
    _geo.N_seg = geo.N_seg;

    int n_seg = _geo.N_seg;

    _U = new double[n_seg];
    _Q = new double[n_seg];

    _cell_hot = new PCHE_CELL[n_seg]; 
    _cell_cold = new PCHE_CELL[n_seg];

    for (size_t i = 0; i < n_seg; i++)
    {
        _cell_hot[i].init(i, 0, this);
        _cell_cold[i].init(i, 1, this);
        
        _U[i] = 0;
        _Q[i] = 0;
    }    

    _A_c = M_PI * _geo.d_c * _geo.d_c / 8;
    _peri_c = _geo.d_c * M_PI / 2 + _geo.d_c ;  
    // Hydraulic Diameter
    // - 1.222 mm for d_c = 2e-3
    // - 0.922 mm for d_c = 1.5e-3
    // or for semi-circular pipe, we have 
    // d_c = ( 1 + 2 / pi) * _d_h
    _d_h = 4 * _A_c / _peri_c;
    _t_wall = ( 2 - M_PI_4) * (_geo.d_c / 2);    
    _length_cell = _geo.length / _geo.N_seg;
    _A_stack =  _peri_c * _length_cell * _geo.N_ch;    
    _A_flow = _geo.N_ch * _A_c;    
}

double PCHE::_avg_T(PCHE_CELL * cell_seq, int count)
{
    if (count == 0)
        return cell_seq[0].T;

    double sum = 0;   

    for (size_t i = 0; i < count; i++)
        sum += cell_seq[i].T;

    return sum / count;    
}

bool PCHE::_calc_U(int idx, BoundaryCondtion & bc)
{
    bool heat_exchanged = true;
    bool cal_next = (idx != _geo.N_seg - 1);

    PCHE_CELL * c_h = &_cell_hot[idx];
    PCHE_CELL * c_c = &_cell_cold[idx];

    _cp_state_update(_handle_cp_hot, _handle_HP_INPUT, c_h->h, c_h->p);
    _cp_state_update(_handle_cp_cold, _handle_HP_INPUT, c_c->h, c_c->p);

    c_h->T = _cp_props(_handle_cp_hot, _handle_T);
    c_c->T = _cp_props(_handle_cp_cold, _handle_T);

    bool valid = c_h->validate();
    valid &= c_c->validate();

    if(valid)
    {
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
        double kw = material_conductivity((c_h->T + c_c->T) / 2); 
        //double kw = 27; 
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

        c_h->dp = c_h->f * _length_cell * c_h->rho * pow(c_h->u, 2) / this->_d_h;
        c_c->dp = c_c->f * _length_cell * c_c->rho * pow(c_c->u, 2) / this->_d_h;

        if(cal_next)
        {
            _cell_hot[idx + 1].p = c_h->p - c_h->dp;
            _cell_hot[idx + 1].h = (bc.st_hot_in.mdot * c_h->h - _Q[idx]) / bc.st_hot_in.mdot; // PropsSI("T", "P", _cell_hot[idx + 1].p, "H", h_hot, medium_hot);
            _cell_cold[idx + 1].p = c_c->p + c_c->dp;
            _cell_cold[idx + 1].h = (bc.st_cold_in.mdot * c_c->h - _Q[idx] ) / bc.st_cold_in.mdot; //PropsSI("T", "P", _cell_cold[idx + 1].p, "H", h_cold, medium_cold);
        }
    }
    else
    {        
        std::string err = string_format("%s !!! invalid cell state\n", _name);
        _log(err, log_level::ERR);      
        _log(c_h->to_string() + c_c->to_string(), log_level::ERR);

        throw std::out_of_range(err);
        /*
        heat_exchanged = false;
        c_h->dp = 0;
        c_c->dp = 0;

        if(cal_next)
        {
            _cell_hot[idx + 1].p = c_h->p;
            _cell_hot[idx + 1].h = c_h->h; 
            _cell_cold[idx + 1].p = c_c->p;
            _cell_cold[idx + 1].h = c_c->h;
        }
        */
    }

    return heat_exchanged;
}

bool PCHE::simulate(const char * media_hot, const char * media_cold, BoundaryCondtion & bc, SIM_PARAM & sim_param, SimulationResult & sr)
{
    std::ostringstream str_stream;
  
    std::string info;

    g_log_level = sim_param.log_level;

    _log(_geo_param_info(_geo), log_level::INFO);
    _log(_boundary_condition_info(bc), log_level::INFO);
    _log(_sim_param_info(sim_param), log_level::INFO);

    bool converged = false;
    
    // cout<<"ready to start simulation"<<endl;

    _handle_cp_hot = _cp_new_state(media_hot);

    _handle_cp_cold = _cp_new_state(media_cold);

    this->_G_hot = bc.st_hot_in.mdot / _geo.N_ch / _A_c;
    this->_G_cold = bc.st_cold_in.mdot / _geo.N_ch / _A_c;

    // temperature difference between the hot inlet and cold inlet
    double dT = 0.0;

    // suppose (p, T) for hot/cold inlet are known
    // try to find out (p, T) for cold/cold outlet by trial iteration
    _cell_hot[0].T = bc.st_hot_in.T;
    _cell_hot[0].p = bc.st_hot_in.p;

    _cell_cold[0].T = _cell_hot[0].T - sim_param.delta_T_init; // start from t_hot_in - 5
    _cell_cold[0].p = bc.st_cold_in.p;

    //double h_hot, h_cold; // specific enthalpy for current hot/cold cell
    double diff_per = 0.0;
    bool heat_exchanged = true;
    int idx_cell = 0;
    try
    {
        for (size_t i = 0; i < sim_param.N_iter; i++)
        {
            heat_exchanged = true;
            _idx_pinch = -1;
            diff_per = 0.0;

            _cp_state_update(_handle_cp_hot, _handle_PT_INPUT, _cell_hot[0].p, _cell_hot[0].T);
            _cp_state_update(_handle_cp_cold, _handle_PT_INPUT, _cell_cold[0].p, _cell_cold[0].T);

            // h_hot = _cp_props(_handle_cp_hot, _handle_H);
            // h_cold = _cp_props(_handle_cp_cold, _handle_H);
            _cell_hot[0].h = _cp_props(_handle_cp_hot, _handle_H);
            _cell_cold[0].h = _cp_props(_handle_cp_cold, _handle_H);

            str_stream<<endl<<_name<<": before "<< i<< " th trial:"<<endl;
            _print_PCHE_state(str_stream.str());

            for(idx_cell = 0; idx_cell < _geo.N_seg; idx_cell ++)
            {
                heat_exchanged = this->_calc_U(idx_cell, bc);
                if(!heat_exchanged)
                    throw std::out_of_range(string_format(" %s no heat exchanged", _name));

                // if(heat_exchanged)
                //     // heat_exchanged = this->calc_U(idx_cell, h_hot, h_cold, bc);
                //     heat_exchanged = this->_calc_U(idx_cell, bc);
                // else 
                // {

                //     //skip the rest of the cells if no heat exchanged between cold and hot side
                //     // T of rest cells remain 0 and the diff in T between bc.T_cold_in and last
                //     // cold cell's T will drive next iteration
                //     //break;
                // }
            }
            // calculate delta T for next iteration
            // dT = bc.st_cold_in->T - _cell_cold[_idx_pinch].T;
            dT = bc.st_cold_in.T - _cell_cold_in()->T;
            converged = the_same(_cell_cold_in()->T, bc.st_cold_in.T, sim_param.err, diff_per);

            str_stream.clear();
            str_stream<<endl<<_name<<": after "<< i<< " th trial:"<<endl;
            _print_PCHE_state(str_stream.str());

            std::string info;

            if(_idx_pinch > 0)
                info = string_format("-> T_{cold, out}=%f K; T_{cold,pinch=%d}=%f K; diff_rel=%f; dT_next= %f K\n", _cell_cold[0].T, _idx_pinch, _cell_cold[_idx_pinch].T, diff_per, dT);
            else
                info = string_format("-> T_{cold, out}=%f K; T_{cold,in}=%f K; diff_rel=%f; dT_next= %f K\n", _cell_cold[0].T, _cell_cold_in()->T, diff_per, dT);

            _log(info, log_level::DEBUG);

            if(converged)
                break;
            // update dt and T_cold_out accordingly, for next iteration
            _cell_cold[0].T += dT * sim_param.step_rel;
            _cell_cold[0].p = bc.st_cold_in.p;
        }

        std::string info = 
        string_format("%s simluation done with result:", _name) +
        _cell_hot[0].to_string() +
        _cell_hot[sr.N_seg - 1].to_string() + 
        _cell_cold[0].to_string() + 
        _cell_cold[sr.N_seg - 1].to_string();

        _log(info, log_level::DEBUG);        
    }
    catch(const std::exception& e)
    {
        // converged = true;

        for (size_t i = idx_cell; i < _geo.N_seg; i++)
        {
            _cell_hot[i].clone(&_cell_hot[i - 1]);
            _cell_cold[i].clone(&_cell_cold[i - 1]);
        }

        _log(string_format("%s clone cell state form %d th cell\n", _name, idx_cell), log_level::DEBUG);
        _print_PCHE_state(info);
        //std::cerr << e.what() << '\n';
    }
    
    // assign the simulation result as output
    for (size_t i = 0; i < sr.N_seg; i++)
    {
        sr.p_hot[i] = _cell_hot[i].p;
        sr.p_cold[i] = _cell_cold[i].p;
        sr.h_hot[i] = _cell_hot[i].h;
        sr.h_cold[i] = _cell_cold[i].h;
    }

    return converged;
}


void PCHE::_cp_state_update(long handle_stream, long handle_input, double val1, double val2)
{
    _err_code = 0;

    AbstractState_update(handle_stream, handle_input, val1, val2,&_err_code, _cp_err_buf, _buffer_size);

    if(_err_code != 0)
    {
        std::string err = string_format("%s error in state update\n", _name);
        _print_cp_err(err.c_str());
        throw std::out_of_range(err.c_str());
    }
}

double PCHE::_cp_props(long handle_stream, long handle_prop, double def_value, double max, double min)
{
    _err_code = 0;

    double prop = AbstractState_keyed_output(handle_stream, handle_prop, &_err_code, _cp_err_buf, _buffer_size);

    if(_err_code != 0)
    {
        // handle the err
        std::string err = string_format("%s error in CoolProp's PropsSI\n", _name);
        _print_cp_err(err);
        prop = def_value;
        throw std::out_of_range(err);
    }

    return prop;
}

long PCHE::_cp_new_state(const char * medium)
{
    _err_code = 0;

    long handle = AbstractState_factory("HEOS", medium, &_err_code, _cp_err_buf, _buffer_size);

    if(_err_code != 0)
        _print_cp_err(string_format("%s error in creating new state", _name));

    return handle;
}

PCHE_CELL * PCHE::_cell_cold_in()
{
    return &_cell_cold[_geo.N_seg - 1];
}

const char * PCHE::_get_prop_name(long handle_prop)
{
    const char * name;
    if(handle_prop == _handle_Dmass)
        name = "Density";
    else if(handle_prop == _handle_H)
        name = "Enthalpy";
    else if(handle_prop == _handle_K)
        name = "Conductivity";
    else if(handle_prop == _handle_MU)
        name = "Vicosity";
    else if(handle_prop == _handle_T)
        name = "Temperature";                        
    else
        name = "Unkown";    

    return name;
}

const char * PCHE::_get_stream_name(long handle_stream)
{
    const char * name;
    if(handle_stream == _handle_cp_cold)
        name = "Cold side";
    else if(handle_stream == _handle_cp_hot)
        name = "Hot side";
    else
        name = "Unkown";

    return name;
}

std::string new_state_string(ThermoState & st)
{
    const char * name;
    // hot_in=0, cold_in=1, hot_out=2, cold_out=3
    if(st.id == 0 )
        name = "hot_in";
    else if(st.id == 1)
        name = "cold_in";
    else if(st.id == 2)
        name = "hot_out";
    else if(st.id == 3)
        name = "cold_out";
    else
        name = "untitled";
    
    return new_state_string(name, st);
}

std::string new_state_string(const char * name, ThermoState & st)
{
    std::ostringstream str_stream;
    str_stream<< "(p, T, h, mdot)_"<<name<<" = ("<<st.p<<", "<<st.T<<", "<<st.h<<", "<<st.mdot<<");"<< endl;
    return str_stream.str();
}

std::string PCHE::_boundary_condition_info(BoundaryCondtion & bc)
{
    std::string info = string_format("**** %s 's boundary condition ****\n", _name);

    return info + new_state_string(bc.st_hot_in) + new_state_string(bc.st_cold_in) + new_state_string(bc.st_hot_out) + new_state_string(bc.st_cold_out);
}

std::string PCHE::_geo_param_info(PCHE_GEO_PARAM & geo)
{
    std::string info = 
    string_format("**** %s's Geometry params ****\n", _name) + 
    string_format("pitch        = %f\n", geo.pitch) +
    string_format("phi          = %f\n", geo.phi) + 
    string_format("length       = %f\n", geo.length) + 
    string_format("length_cell  = %f\n", _length_cell) + 
    string_format("d_c          = %f\n", geo.d_c) + 
    string_format("N_ch         = %d\n", geo.N_ch) + 
    string_format("N_seg        = %d\n", geo.N_seg);

    return info;
}

std::string PCHE::_sim_param_info(SIM_PARAM & sim_param)
{
    std::string info = 
    string_format("**** %s 's simulation params ****\n", _name) + 
    string_format("err          = %f\n", sim_param.err) +
    string_format("dT_init      = %f\n", sim_param.delta_T_init) + 
    string_format("N_iter       = %d\n", sim_param.N_iter) + 
    string_format("step_rel     = %f\n", sim_param.step_rel) + 
    string_format("log_level    = %d\n", sim_param.log_level);

    return info;
}

void PCHE::_print_PCHE_state(std::string title)
{
    int idx_last = _geo.N_seg - 1;

    std::string info = title;
    info += string_format("**** %s 's state of four ports (p, T, h) ****\n", _name);
    info += string_format("st_hot_in    = %s\n", _cell_hot[0].to_string());
    info += string_format("st_hot_out   = %s\n", _cell_hot[idx_last].to_string());
    info += string_format("st_cold_out  = %s\n", _cell_cold[0].to_string());
    info += string_format("st_cold_in   = %s\n", _cell_cold[idx_last].to_string());
    info += string_format("G_hot/cold   = %f/%f\n", _G_hot, _G_cold);

    _log(info, log_level::DEBUG);
}

void PCHE::_print_cp_err(std::string err)
{
    _print_cp_err(err.c_str());
}

void PCHE::_print_cp_err(const char * info)
{
    _log(string_format("%s(code=%d): %s\n", info, _err_code, _cp_err_buf));
}

void PCHE::_log(std::string content, log_level level)
{
    if(level < g_log_level) return;

    cout<<content;

}