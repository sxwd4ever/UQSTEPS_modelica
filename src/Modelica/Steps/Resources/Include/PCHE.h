#ifndef __PCHE_H
#define __PCHE_H

#include <iostream>

#define MAX_CELL_NUM 1000

using namespace std;

enum log_level {DEBUG = 0, INFO = 1, ERR = 2, SERVE = 3, OFF = 4};

struct ThermoState
{
    double T;
    double p;
    double h;
    double mdot;
    // optional id, to specify which point this state is assigned to 
    // because of the mapping constraintns of Modelica and external functions
    // integer id is used instead of string(char *) name
    // more detailed naming specificaiton should be provided
    // hot_in=0, cold_in=1, hot_out=2, cold_out=3
    int id = -1; 
};

/**
 * basic geometric parameters of PCHE
 */
struct PCHE_GEO_PARAM
{
    /* data */
    // pitch length
    double pitch;
    // pitch angle
    double phi;
    // length of segment
    double length;
    // Diameter of semi_circular
    double d_c;
    // number of channels
    int N_ch;
    // number of segments
    int N_seg; // [Kwon2019]'s maximum node number   
};

/** 
 * Kim's Correlation coeffients
 */
struct KIM_CORR_COE
{
    /* data */
    double a;
    double b;
    double c;
    double d;
};


struct SIM_PARAM
{
    // error tolerance
    double err;
    // initial delta_T between T_hot_in and T_cold_out
    double delta_T_init;
    // maximum iteration number
    int N_iter;
    // relative step length between two trial values of T_cold_out
    // T_cold_out[i] += step_rel * (T_bc_cold_in  - T_cold_in[i-1])
    double step_rel; 

    // log level for output contorl
    int log_level = log_level::OFF;
};

typedef struct
{
    ThermoState st_hot_in;
    ThermoState st_cold_in;
    ThermoState st_hot_out;
    ThermoState st_cold_out;
} BoundaryCondtion;

// typedef struct
// {
//     ThermoState st_hot_in;
//     ThermoState st_cold_in;
//     ThermoState st_hot_out;
//     ThermoState st_cold_out;
// }SimulationResult;

typedef struct
{
  double * h_hot;
  double * h_cold;
  
  double * p_hot;
  double * p_cold;

  /* segment number */
  int N_seg;

} SimulationResult;

class PCHE;

class PCHE_CELL
{
private:
    /* data */
    PCHE * _pche;
    const char * type_name();
public:
    PCHE_CELL();
    ~PCHE_CELL();

    // type: 0=hot, 1=cold
    void init(int idx, int type, PCHE * pche);
    void clone(PCHE_CELL * src);

    void print_state();

    bool validate();

    // idx in sequence
    int idx; 
    // hot/cold side
    int type_stream;    
    // heat Flux
    double Q;
    // heat flux
    double G;
    // local temperature
    double T;
    // local pressure
    double p;
    // local specific enthalpy
    double h;
    // local velocity of fluid
    double u;
    // local dynamic vsicosity
    double mu;
    // local conductivity
    double k;
    // local renorlds number
    double Re;
    // local density
    double rho;
    // local Nusselt Number
    double Nu;
    // local Thermal Conductance
    double hc;
    // Fanning friction factor
    double f;
    // local pressure drop
    double dp;    
};

class PCHE
{
private:
    /* data */
    void _init(PCHE_GEO_PARAM & geo);

    PCHE_CELL * _cell_cold_in();

    double _cp_props(long handle_stream, long handle_prop, double def_value = 0, double max = 1e6, double min = 0);

    long _cp_new_state(const char * medium);

    void _cp_state_update(long handle_stream, long handle_input, double val1, double val2);

    bool _calc_U(int idx, BoundaryCondtion & bc);

    double _avg_T(PCHE_CELL * cell_seq, int count);

    /*** function outputs parameters for debug ***/
    bool _can_log(int level);

    void _print_boundary_conditions(BoundaryCondtion & bc);

    void _print_geo_param(PCHE_GEO_PARAM & geo);

    void _print_sim_param(SIM_PARAM & sim_param);

    void _print_PCHE_state(std::string title);

    void _print_cp_err(const char * info);

    const char * _get_prop_name(long handle_prop);

    const char * _get_stream_name(long handle_stream);
    /*** end of debug utility functions ***/

    // Kim's correlation factors
    KIM_CORR_COE _corr_coe;

    /* variables for cool prop*/
    // handles for querry 
    long _handle_cp_hot;
    long _handle_cp_cold;

    long _handle_HP_INPUT;
    long _handle_PT_INPUT;

    long _handle_T;
    long _handle_H;
    long _handle_Dmass;
    long _handle_MU;
    long _handle_K;

    /* for CP's error handle */
    long _err_code;
    char * _cp_err_buf;
    long _buffer_size;

    // *** geometry parameter ***    
    PCHE_GEO_PARAM _geo;

    // Area of semi-circular tube - area per channel
    double _A_c;
    // perimeter of semi-circular
    double _peri_c;  
    // Hydraulic Diameter
    double _d_h;
    // thickness of wall between two neighboring hot and cold
    double _t_wall;
    // surface area of all cells in a stack 
    double _A_stack;
    // Flow area for all channels
    double _A_flow;
    // length of one segment in a tube, m
    double _length_cell;

    // mass flux for cold/hot side    
    double _G_hot, _G_cold;

    PCHE_CELL * _cell_hot;
    PCHE_CELL * _cell_cold;
    
    // index of the pinch point cell
    int _idx_pinch;

    const char * _name;
    
public:
    PCHE(const char * name, PCHE_GEO_PARAM & geo);
    ~PCHE();

    /**
     * set kim's correlation coefficients
     */
    void set_kim_corr_coe(KIM_CORR_COE & coe);

    bool simulate(const char * media_hot, const char * media_cold, BoundaryCondtion & bc, SIM_PARAM & sim_param, SimulationResult & sr);

    // overall heat transfer coefficients
    double * _U;
    // array of local transferred heat 
    double * _Q;    
};

std::string new_state_string(ThermoState & st);

std::string new_state_string(const char * name, ThermoState & st);

#endif