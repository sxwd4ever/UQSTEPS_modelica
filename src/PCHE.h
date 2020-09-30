#ifndef __PCHE_H
#define __PCHE_H
#define MAX_SEG_LEN 10000

class PCHE_CELL
{
private:
    /* data */
public:
    PCHE_CELL(double p = 0, double T = 0, const char * medium="CO2");
    ~PCHE_CELL();

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
public:
    PCHE(/* args */);
    ~PCHE();

    // *** geometry parameter ***    
    // pitch length
    double pitch;
    // angle
    double phi;
    // Kim's correlation factors
    double a = 0, b = 0, c = 0, d = 0;
    // length of each segment
    double length_cell = 12e-3;
    // Diameter of semi_circular
    double d_c = 2e-3;
    // Area of semi-circular tube - area per channel
    double A_c;
    // perimeter of semi-circular
    double peri_c;  
    // Hydraulic Diameter
    double d_h;
    // thickness of wall between two neighboring hot and cold
    double t_wall;
    // number of channels
    int N_ch = 80e3;
    // number of segments
    int N_seg = 100;
    // surface area of all cells in a stack 
    double A_stack;
    // Flow area for all channels
    double A_flow;
    // length of one pipe in HeatExchanger unit m
    double length_ch;   

    PCHE_CELL cell_hot[MAX_SEG_LEN];
    PCHE_CELL cell_cold[MAX_SEG_LEN];
};



#endif