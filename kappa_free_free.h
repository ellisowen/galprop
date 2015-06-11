#ifndef KAPPA_FREE_FREE_h
#define KAPPA_FREE_FREE_h


#include <valarray>



inline double gff(const double nu, const double Te);
inline double kff(const double nu, const double Ne, const double Te);
inline double eff(const double nu, const double Ne, const double Te);

double kappa_free_free(double nu,double Ne, double Te,  double &emiss_free_free, int options, int debug); //AWS20110704
double kappa_free_free_2D(double nu, double R,           double z,  double Te, double clumping_factor, double &emiss_free_free, int options, int debug); //AWS20110704
double kappa_free_free_3D(double nu, double x, double y, double z,  double Te, double clumping_factor, double &emiss_free_free, int options, int debug); //AWS20110704



#endif
