#include "CoolPropLib.h"
#include "AbstractState.h"
#include "DataStructures.h"
#include "MyPropsLib.h"
#include "PCHE.h"
#include <math.h>
#include <exception>
#include <sstream>
#include <memory>
#include <string>
#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <vector>

// #include <string>
// #include "crossplatform_shared_ptr.h"
#include <tr1/memory>
// #include "example_dll.h"

using namespace CoolProp;

static int g_log_level = log_level::OFF;

class PCHELibrary{
private:
    std::map<std::size_t, std::tr1::shared_ptr<PCHE> > _PCHE_library;
    long _next_handle;
public:
    PCHELibrary(): _next_handle(0){};
    long add(std::tr1::shared_ptr<PCHE> pche){
        _PCHE_library.insert(std::pair<std::size_t, std::tr1::shared_ptr<PCHE> >(_next_handle,  pche));
        _next_handle++;
        return _next_handle-1;
    }
    void remove(long handle){
        std::size_t count_removed = _PCHE_library.erase(handle);
        if (count_removed != 1){
            cout<<"could not free handle="<<handle<<endl;
        }
    }
    std::tr1::shared_ptr<PCHE> & get(long handle){
        std::map<std::size_t, std::tr1::shared_ptr<PCHE> >::iterator it = _PCHE_library.find(handle);
        if (it != _PCHE_library.end()){
            return it->second;
        }
        else{
            throw exception();
        }
    }
};

static PCHELibrary handle_manager;

/**
 * convert angle from degree to rad
 */
double EXPORT_MY_CODE from_deg(double deg)
{
    return deg * M_PI / 180;
}

/** 
 * convert preseaure in bar to Pa
 */
double EXPORT_MY_CODE from_bar(double p_bar)
{
    return p_bar * 1e5;
}

double EXPORT_MY_CODE from_degC(double degC)
{
    return degC + 273.15;
}

bool EXPORT_MY_CODE the_same(double x, double y, double eps, double & diff_per)
{
    double err = 1e-4;
    if (abs(x - 0) < err)
    {
        diff_per = 1.0;
        return abs(y) <= err;
    }
        
    diff_per = abs((y-x) * 1.0 / x);
    return  diff_per<= eps;
}



double EXPORT_MY_CODE MyPropsSI(const char *Output, const char *Name1, double Prop1, const char *Name2, double Prop2, const char *Ref)
{
    return PropsSI(Output, Name1, Prop1, Name2, Prop2, Ref);
}

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

    if( T < 0 || T >= 3000|| mu < 0 || mu > 50e3|| k < 0 || k > 10e3 || rho < 0 || rho > 10e3)
    {
        // print error state values for debug
        cout<<string_format("(T, mu, k, rho)_(%d, %d) = (%d, %d, %d, %d)\n", p, H, T, mu, k, rho);
    }

	return T;
}

ThermoState * EXPORT_MY_CODE NewThermoState_pT(double p, double T, double mdot, const char * medium)
{
    ThermoState * st = new ThermoState();
    st->p = p;
    st->T = T;
    st->h = PropsSI("H", "P", p, "T", T, medium);
    st->mdot = mdot;
    return st;
}

void CompleteThermoState(ThermoState * st, const char * media)
{
    // st->p *= 1e6;
    // st->h *= 1e6;
    // st->h *= 1e6;

    if(st->T == 0)
        st->T = PropsSI("T", "P", st->p, "H", st->h, media);
    else if(st->h == 0)
        st->h = PropsSI("H", "P", st->p, "T", st->T, media);

    // modelica may transfer negetive mass flow rate - correct it here
    if(st->mdot < 0)
        st->mdot = -st->mdot;

    if(g_log_level == log_level::DEBUG)
        cout<<new_state_string(*st);
}

std::string new_state_string_ph(double p, double h, const char * media)
{
    std::ostringstream str_stream;
    str_stream<<"("<<p<<", "<< PropsSI("T", "P", p, "H", h, media) - 273.15<<", "<<h<<")";

    return str_stream.str();
}

double EXPORT_MY_CODE PCHE_OFFD_Simulation(const char * name, const char * media_hot, const char * media_cold, PCHE_GEO_PARAM * geo, KIM_CORR_COE * cor, SIM_PARAM * sim_param, BoundaryCondtion * bc, PCHECImplResult * retOutput)
{
    double * U = new double[geo->N_seg];
    double * Q = new double[geo->N_seg];

    double retVal = PCHE_OFFD_Simulation_UQ_out(name, media_hot, media_cold, geo, cor, sim_param, bc, retOutput, Q, U);
    
    delete [] U;
    delete [] Q;

    return retVal;
}


/**
 * off-design simulation for PCHE
 */
double EXPORT_MY_CODE PCHE_OFFD_Simulation_UQ_out(const char * name, const char * media_hot, const char * media_cold, PCHE_GEO_PARAM * geo, KIM_CORR_COE * cor, SIM_PARAM * sim_param, BoundaryCondtion * bc, PCHECImplResult * retOutput, double * Q, double * U)
{
    int N_seg = geo->N_seg;
    bool done = false;

    g_log_level = sim_param->log_level;

    if(N_seg <= 0)
    {
        if(g_log_level == log_level::DEBUG)
            cout<< "invalid N_seg="<<N_seg<<"; geo.N_seg="<< geo->N_seg<<endl;
        return -1;
    }
    else
    {
        if(g_log_level == log_level::DEBUG)
            cout<<name<<" in dll now N_seg="<<N_seg<<"; geo.N_seg="<< geo->N_seg<<endl;
    }

    SimulationResult sr;
    sr.h_hot = new double[N_seg];
    sr.h_cold = new double[N_seg];
    sr.p_hot = new double[N_seg];
    sr.p_cold = new double[N_seg];
    sr.N_seg = N_seg;

    if(g_log_level == log_level::DEBUG)
        cout<<name<<" initializing... "<<endl;

    PCHE * pche = new PCHE(name, *geo);    

    if(g_log_level == log_level::DEBUG)
        cout<<name<<" initialization done"<<endl;

    pche->set_kim_corr_coe(*cor);
    
    // bc->media_cold = "CO2";
    // bc->media_hot = "CO2";

    // fill the T in bc

    CompleteThermoState(& bc->st_hot_in, media_hot);
    CompleteThermoState(& bc->st_hot_out, media_hot);
    CompleteThermoState(& bc->st_cold_in, media_cold);
    CompleteThermoState(& bc->st_cold_out, media_cold);    

    if(g_log_level == log_level::DEBUG)
        cout<<"ready to simulate "<<name<<endl;

    done = pche->simulate(media_hot, media_cold, * bc, * sim_param, sr);

    int last = N_seg - 1;


    retOutput->p_hot = sr.p_hot[last];
    retOutput->h_hot = sr.h_hot[last];

    retOutput->p_cold = sr.p_cold[0];
    retOutput->h_cold = sr.h_cold[0];

    for (size_t i = 0; i < N_seg; i++)
    {
        Q[i] = pche->_Q[i];
        U[i] = pche->_U[i];
    }    

    // if (done)
    // {
    // // size of returned array is fixed at 2. 
    //     int last = N_seg - 1;

    //     h_hot[0] = sr.h_hot[0];
    //     h_hot[1] = sr.h_hot[last];
    //     h_cold[0] = sr.h_cold[0];
    //     h_cold[1] = sr.h_cold[last];

    //     p_hot[0] = sr.p_hot[0];
    //     p_hot[1] = sr.p_hot[last];
    //     p_cold[0] = sr.p_cold[0];
    //     p_cold[1] = sr.p_cold[last];
    // }
    // else 
    // {
    //     // use input boundary condition
    //     // treat pche as a direct pipe
    //     h_hot[0] = bc->st_hot_in.h;
    //     h_hot[1] = bc->st_hot_in.h;
    //     h_cold[0] = bc->st_cold_in.h;
    //     h_cold[1] = bc->st_cold_in.h; //sr.h_cold[last];

    //     p_hot[0] = bc->st_hot_in.p; //sr.p_hot[0];
    //     p_hot[1] = bc->st_hot_in.p; //sr.p_hot[last];
    //     p_cold[0] = bc->st_cold_in.p; //sr.p_cold[0];
    //     p_cold[1] = bc->st_cold_in.p; //sr.p_cold[last];

    // }

    if(g_log_level <= log_level::INFO)
    {
        const char * result = done ? "OK": "FAILED";
        cout<<name<<"("<<result<<"): [p, T, h] # MPa|oC|J/K @ ";
        cout<<"hi="<<new_state_string_ph(sr.p_hot[0], sr.h_hot[0], media_hot)<<"; ";
        cout<<"ci="<<new_state_string_ph(sr.p_cold[last], sr.h_cold[last], media_cold)<<"; ";
        cout<<"ho="<<new_state_string_ph(retOutput->p_hot, retOutput->h_hot, media_hot)<<"; ";
        cout<<"co="<<new_state_string_ph(retOutput->p_cold, retOutput->h_cold, media_cold)<<endl;
    }

    delete pche;
    delete [] sr.h_hot;
    delete [] sr.h_cold;
    delete [] sr.p_hot;
    delete [] sr.p_cold;
    /*
    // change magnitude 
    const int mag = 1e6;
    // scaling the magnitude
    retOutput->p_hot /= mag;
    retOutput->h_hot /= mag;

    retOutput->p_cold /= mag;
    retOutput->h_cold /= mag;
    */    
	return 0.0;
}

double EXPORT_MY_CODE print_path_state(const char * name, const char * media, ThermoState * st, int log_level)
{
    if(g_log_level <= log_level)
    {
        std::string info = 
        string_format("at %s: (p, T, h)=%s\n", name, 
        new_state_string_ph(st->p, st->h, media).c_str());
        cout<<info;
    }

    return 0.0;
}

#include "CoolPropLib.h"

boolean validate_cp_result(CoolPropWrapper * cp_wrapper)
{
    if(cp_wrapper->err_code != 0)
    {
        // my_log(string_format("error in creating new state for %s: %s", cp_wrapper->name, cp_wrapper->cp_err_buf), log_level::ERR);
        cout<<string_format("error in creating new state for %s: %s\n", cp_wrapper->name, cp_wrapper->cp_err_buf);
        return false;
    }

    return true;
}

void * EXPORT_MY_CODE init_cp_wrapper(const char * medium, const char * name)
{
    CoolPropWrapper * cp_wrapper = new CoolPropWrapper();

    // for coolprop query
    cp_wrapper->handle_cp_state = 0;
    cp_wrapper->buffer_size = 500;
    cp_wrapper->err_code = 0;
    cp_wrapper->cp_err_buf = new char[cp_wrapper->buffer_size];
    cp_wrapper->name = name;

    cp_wrapper->handle_cp_state = AbstractState_factory("HEOS", medium, &cp_wrapper->err_code, cp_wrapper->cp_err_buf, cp_wrapper->buffer_size);

    if(!validate_cp_result(cp_wrapper)) return NULL;
    
    return (void *)cp_wrapper;
}

double EXPORT_MY_CODE cp_query(void * wrapper, const char * input_pair,  double val1, double val2, const char * output_name)
{
    CoolPropWrapper * cp_wrapper = (CoolPropWrapper *)wrapper;

    cp_wrapper->err_code = 0;
    long handle_input =  get_input_pair_index(input_pair);   

    long handel_output = get_param_index(output_name);

    AbstractState_update(cp_wrapper->handle_cp_state, handle_input, val1, val2, &cp_wrapper->err_code, cp_wrapper->cp_err_buf, cp_wrapper->buffer_size);

    if(!validate_cp_result(cp_wrapper)) return -1;

    double prop = AbstractState_keyed_output(cp_wrapper->handle_cp_state, handel_output, &cp_wrapper->err_code, cp_wrapper->cp_err_buf, cp_wrapper->buffer_size);

    if(!validate_cp_result(cp_wrapper)) return -1;

    return prop;
}

void EXPORT_MY_CODE close_cp_wrapper(void * state)
{
    CoolPropWrapper * cp_wrapper = (CoolPropWrapper *) state;

    cp_wrapper->err_code = 0;

    if(cp_wrapper->handle_cp_state != 0)
        AbstractState_free(cp_wrapper->handle_cp_state, & cp_wrapper->err_code, cp_wrapper->cp_err_buf, cp_wrapper->buffer_size);

    validate_cp_result(cp_wrapper);

    delete [] cp_wrapper->cp_err_buf;

    delete cp_wrapper;
}

/**
 * test for transferring c struct as input/output parameter


void * EXPORT_MY_CODE init_PCHE_sim_ext_object(PCHE_GEO_PARAM * geo)
{
    std::tr1::shared_ptr<PCHE> ptr_pche (new PCHE(*geo)); 

    PCHE_SIM_EXT_OBJ * ext_obj = new PCHE_SIM_EXT_OBJ();
    ext_obj->handle_pche = handle_manager.add(ptr_pche);

    return (void *) ext_obj;
}

void EXPORT_MY_CODE PCHE_simulate(SIM_PARAM * sim_param, BoundaryCondtion * bc, void * ext_obj)
{
    PCHE_SIM_EXT_OBJ * pche_ext_obj = (PCHE_SIM_EXT_OBJ *)ext_obj;

    std::tr1::shared_ptr<PCHE> & pche = handle_manager.get(pche_ext_obj->handle_pche);
    
    // pche->simulate(*bc,*sim_param);

    pche_ext_obj->st_cold_in.T = bc->st_cold_in.T + 123;
}

void EXPORT_MY_CODE close_PCHE_sim_ext_object(void * ext_obj)
{
    PCHE_SIM_EXT_OBJ * pche_ext_obj = (PCHE_SIM_EXT_OBJ *)ext_obj;

    std::tr1::shared_ptr<PCHE> & pche = handle_manager.get(pche_ext_obj->handle_pche);
    
    handle_manager.remove(pche_ext_obj->handle_pche);
}
*/

double EXPORT_MY_CODE material_conductivity(double T, bool extrapolate)
{
    // Returns interpolated value at x from parallel arrays ( xData, yData )
    //   Assumes that xData has at least two elements, is sorted and is strictly monotonic increasing
    //   boolean argument extrapolate determines behaviour beyond ends of array (if needed)

    // Original data
    std::vector<double> xData = { 149, 316, 538, 649, 760, 871};
    std::vector<double> yData = { 16.9, 20.5, 26.5, 28.7, 31.4, 35.3};

    int size = xData.size();
    double x = T;

    int i = 0;                                                                  // find left end of interval for interpolation
    if ( x >= xData[size - 2] )                                                 // special case: beyond right end
    {
    i = size - 2;
    }
    else
    {
    while ( x > xData[i+1] ) i++;
    }
    double xL = xData[i], yL = yData[i], xR = xData[i+1], yR = yData[i+1];      // points on either side (unless beyond ends)
    if ( !extrapolate )                                                         // if beyond ends of array and not extrapolating
    {
    if ( x < xL ) yR = yL;
    if ( x > xR ) yL = yR;
    }

    double dydx = ( yR - yL ) / ( xR - xL );                                    // gradient

    return yL + dydx * ( x - xL );                                              // linear interpolation    
}

void EXPORT_MY_CODE test_struct_param(SIM_PARAM * sim_param, PCHE_GEO_PARAM * geo, BoundaryCondtion * bc, double * h_hot, double * h_cold, double * p_hot, double * p_cold, size_t N_seg)
{
    // sim_param = new SIM_PARAM();

    sim_param->err = geo->pitch;
    sim_param->step_rel = geo->phi + geo->length;
    sim_param->delta_T_init = 10.0;

    // st_result[1] = 11222;

    for (size_t i = 0; i < N_seg; i++)
    {
        /* code */
        //st_result[i].T = i * 100 + 1;
        //st_result[i].p = i * 1000 + 2;
        // st_result[i] = i * 100 + 1;
        h_hot[i] = i * 100 + 1;
        h_cold[i] = i * 100 + 3;
        p_hot[i] = i * 100 + 5;
        p_cold[i] = i * 100 + 7;
    }    

    // return sim_param;
}

void EXPORT_MY_CODE setState_C_impl(double p, double M,  State *state)
{
    state->p = p;
    state->M = M;
}

void my_log(std::string content, log_level level)
{
    if(level < g_log_level) return;

    cout<<content;

}
