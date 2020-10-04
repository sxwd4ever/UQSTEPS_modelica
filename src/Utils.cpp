#include "Utils.h"
#include <math.h>

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