#ifndef MY_RPOPS_H__
#define MY_RPOPS_H__

#ifdef __cplusplus
extern "C" {
#endif

// #define _GLIBCXX_USE_CXX11_ABI 0

void MyPropsSI_pT(double p, double T, const std::string &FluidName, double &h, double &rhomass);

double MyPropsSI_pH(double p, double H, const std::string &FluidName, double &mu, double &k, double &rho);

#ifdef __cplusplus
}
#endif

#endif