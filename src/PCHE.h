#ifndef __PCHE_H
#define __PCHE_H

#include <iostream>

using namespace std;

struct ThermoState
{
    double T;
    double p;
    double h;
    double mdot;
    std::string medium;
};

/**
 * basic geometric parameters of PCHE
 */
struct PCHE_GEO_PARAM
{
    /* data */
    // pitch length
    double pitch;
    // angle
    double phi;
    // length of each segment
    double length_cell;
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
    double err;
    int N_iter;
    // relative step length between two trial values of T_cold_out
    // T_cold_out[i] += step_rel * (T_bc_cold_in  - T_cold_in[i-1])
    double step_rel; 
};

class BoundaryCondtion
{
private:
    /* data */
public:
    BoundaryCondtion(/* args */);
    ~BoundaryCondtion();

    ThermoState * st_hot_in;
    ThermoState * st_cold_in;
    ThermoState * st_hot_out;
    ThermoState * st_cold_out;
};


class PCHE;

class PCHE_CELL
{
private:
    /* data */
    PCHE * _pche;

public:
    PCHE_CELL();
    ~PCHE_CELL();

    void _init(ThermoState * st, PCHE * pche, std::string medium = "CO2");

    void clone(PCHE_CELL * src);

    /* data */
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
    // length of one pipe in HeatExchanger unit m
    double _length_ch;

    // mass flux for cold/hot side    
    double _G_hot, _G_cold;

    PCHE_CELL * _cell_hot;
    PCHE_CELL * _cell_cold;
    
    // overall heat transfer coefficients
    double * _U;
    // local transferred heat 
    double * _q;

public:
    PCHE(PCHE_GEO_PARAM & geo);
    ~PCHE();

    /**
     * set kim's correlation coefficients
     */
    void set_kim_corr_coe(KIM_CORR_COE & coe);

    void simulate(BoundaryCondtion & bc, SIM_PARAM & sim_para);

    bool calc_U(int idx, double & h_hot, double & h_cold, BoundaryCondtion & bc);

    double avg_T(PCHE_CELL * cell_seq, int count);
};


ThermoState * NewThermoState_pT(double p, double T, double mdot, std::string medium);

#endif