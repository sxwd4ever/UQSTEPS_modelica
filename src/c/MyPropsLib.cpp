#include "CoolPropLib.h"
#include "AbstractState.h"
#include "DataStructures.h"
#include "MyPropsLib.h"
#include "PCHE.h"
#include <math.h>
#include <exception>

// #include <string>
// #include "crossplatform_shared_ptr.h"
#include <tr1/memory>
// #include "example_dll.h"

using namespace CoolProp;

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

/**
 * off-design simulation for PCHE
 */
double EXPORT_MY_CODE PCHE_OFFD_Simulation(const char * media_hot, const char * media_cold, PCHE_GEO_PARAM * geo, KIM_CORR_COE * cor, SIM_PARAM * sim_param, BoundaryCondtion * bc, double * h_hot, double * h_cold, double * p_hot, double * p_cold, size_t N_seg)
{
    SimulationResult sr;
    sr.h_hot = h_hot;
    sr.h_cold = h_cold;
    sr.p_hot = p_hot;    
    sr.p_cold = p_cold;
    sr.N_seg = N_seg;

    PCHE * pche = new PCHE(*geo);

    pche->set_kim_corr_coe(*cor);

    // bc->media_cold = "CO2";
    // bc->media_hot = "CO2";

    pche->simulate(media_hot, media_cold, * bc, * sim_param, sr);


    delete pche;

	return 0.0;
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

void EXPORT_MY_CODE PCHE_simulate(SIM_PARAM * sim_para, BoundaryCondtion * bc, void * ext_obj)
{
    PCHE_SIM_EXT_OBJ * pche_ext_obj = (PCHE_SIM_EXT_OBJ *)ext_obj;

    std::tr1::shared_ptr<PCHE> & pche = handle_manager.get(pche_ext_obj->handle_pche);
    
    // pche->simulate(*bc,*sim_para);

    pche_ext_obj->st_cold_in.T = bc->st_cold_in.T + 123;
}

void EXPORT_MY_CODE close_PCHE_sim_ext_object(void * ext_obj)
{
    PCHE_SIM_EXT_OBJ * pche_ext_obj = (PCHE_SIM_EXT_OBJ *)ext_obj;

    std::tr1::shared_ptr<PCHE> & pche = handle_manager.get(pche_ext_obj->handle_pche);
    
    handle_manager.remove(pche_ext_obj->handle_pche);
}
*/

void EXPORT_MY_CODE test_struct_param(SIM_PARAM * sim_para, PCHE_GEO_PARAM * geo, BoundaryCondtion * bc, double * h_hot, double * h_cold, double * p_hot, double * p_cold, size_t N_seg)
{
    // sim_para = new SIM_PARAM();

    sim_para->err = geo->pitch;
    sim_para->step_rel = geo->phi + geo->length_cell;
    sim_para->delta_T_init = 10.0;

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

    // return sim_para;
}

void EXPORT_MY_CODE setState_C_impl(double p, double M,  State *state)
{
    state->p = p;
    state->M = M;
}
