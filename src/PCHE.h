#ifndef __PCHE_H
#define __PCHE_H
#define MAX_SEG_LEN 10000

class ThermoState
{
private:
    /* data */
public:
    ThermoState(double p, double T, double mdot, std::string medium);
    ~ThermoState();

    double T;
    double h;
    double p;
    double mdot;

    // medium name
    std::string medium;
};

struct SIM_PARAM
{
    double err;
    int N_iter;
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

    void init(ThermoState & st, PCHE * pche, std::string medium = "CO2");

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
    void init();
    /* variables for cool prop*/

    char * _cp_err_buf = NULL;
    long _handle_cp_hot;
    long _handle_cp_cold;
    long _handle_HP_INPUT;
    long _handle_T;
    long _handle_Dmass;
    long _handle_MU;
    long _handle_K;
    long _err_code;
public:
    PCHE(/* args */);
    ~PCHE();

    void simulate(BoundaryCondtion & bc, SIM_PARAM & sim_para);

    bool calc_U(int idx, double & h_hot, double & h_cold, BoundaryCondtion & bc);

    double avg_T(PCHE_CELL * cell_seq, int count);
    // *** geometry parameter ***    
    // pitch length
    double pitch;
    // angle
    double phi;
    // Kim's correlation factors
    double a, b, c, d;
    // length of each segment
    double length_cell;
    // Diameter of semi_circular
    double d_c;
    // Area of semi-circular tube - area per channel
    double A_c;
    // perimeter of semi-circular
    double peri_c;  
    // Hydraulic Diameter
    double d_h;
    // thickness of wall between two neighboring hot and cold
    double t_wall;
    // number of channels
    int N_ch;
    // number of segments
    int N_seg;
    // surface area of all cells in a stack 
    double A_stack;
    // Flow area for all channels
    double A_flow;
    // length of one pipe in HeatExchanger unit m
    double length_ch;   

    double G_hot, G_cold;

    PCHE_CELL cell_hot[MAX_SEG_LEN];
    PCHE_CELL cell_cold[MAX_SEG_LEN];
    
    // overall heat transfer coefficients
    double U[MAX_SEG_LEN];
    // local transferred heat 
    double q[MAX_SEG_LEN];

};



#endif