#ifndef MY_RPOPS_H__
#define MY_RPOPS_H__

#include <string>

void MyPropsSI_pT(double p, double T, const std::string &FluidName, double &h, double &rhomass);

void MyPropsSI_pH(double p, double H, const std::string &FluidName, double &T, double &mu, double &k, double &rho);

#endif