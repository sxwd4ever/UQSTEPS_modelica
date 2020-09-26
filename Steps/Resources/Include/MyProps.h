#ifndef MY_RPOPS_H__
#define MY_RPOPS_H__

#include <string>

void MyPropsSI_pT(double p, double T, const std::string &FluidName, double &h, double &rhomass);

double MyPropsSI_pH(double p, double H, const std::string &FluidName, double &mu, double &k, double &rho);

#endif