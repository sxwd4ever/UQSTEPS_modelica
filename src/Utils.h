#ifndef __UTILS_H
#define __UTILS_H

double from_deg(double deg);

/** 
 * convert preseaure in bar to Pa
 */
double from_bar(double p_bar);

double from_degC(double degC);

bool the_same(double x, double y, double eps, double & diff_per);

#endif
