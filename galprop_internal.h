//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * galprop.h *                                   galprop package * 08/16/2001 
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|
#ifndef _galprop_internal_h_
#define _galprop_internal_h_

#include <constants.h>

#include <cmath>
#include <iostream>
#include <string>
#include <valarray>
#include "GasCubeData.h"

using namespace std;                                             //AWS20050912

//#include "fort_interface.h"

// function prototypes

// only those functions which are not in Galprop class
// and hence stand-alone.

double eprop(double,double,double,double,double,double,double,double,double,double);//IMOS20061030
void Kcapture_cs(double,int,int,double*,double*);               // IMOS20010816
void nucleon_cs(int,double,int,int,int,double*,double*,double*,double*,double*,double*);// IMOS20010511
double isotope_cs(double,int,int,int,int,int,int*);
void read_nucdata(const string& path);
void cleanup_nucdata();                                         //IMOS20060420
double nucdata(int,int,int,int,int,int,int*,int*,double*);      // IMOS20010816
double nucleon_loss(int,int,double,double,double,double,        // nucleon
		    double*, double*);                                    // energy losses
double electron_loss(double,double,double,double,double,double, // electron
		     double*,double*,double*,double*,double*,double*);     // energy losses
double blattnig_gamma(double,double,int,int,int);               // gammas from pi0-decay, Blattnig etal. formalism

double kamae(double,double,int,int,int, valarray<double>, valarray<double>, valarray<double>, valarray<double>);//double**);              // gammas from pi0-decay, Kamae etal. formalism IMOS20061114
 // functions for component parameter calculation IMOS20061114
valarray<double> kamae_gamma_param_nd(double);//, double*);
valarray<double> kamae_gamma_param_diff(double);//, double*);
valarray<double> kamae_gamma_param_delta(double);//, double*);
valarray<double> kamae_gamma_param_res(double);//, double*);

double sigma_boron_dec_heinbach_simon(int, int, int, int, double);

//double sim(double,double,double,double,double,double(*)(double)); //integration
int kinematic(int, int, char*, double&, double&, double&, double&, double&, double&, int);

int tridag(float*, float*, float*, float*, float*, int);
int tridag(double*, double*, double*, double*, double*, int); //IMOS20030217
int tridag_sym(float*, float*, float*, float*, float*, int);
int tridag_sym(double*, double*, double*, double*, double*, int); //IMOS20030217
double sim(double, double, double, double, double, double(*)(double)); //integration

int He_to_H_CS(double, int, int, int, int, double*, double*);

double aic_cc(int, int, double, double, double, double, double, double, double, double, double);//IMOS20060420
double fjones_cc(double, double, double);          //IMOS20060420

double ionization_bethe(int Z, double beta);

double nHI(double, double);
//double nHI_3D(double, double,double);//ECC20141024
double nH2(double, double);
double nH2(double, double, double);//AWS20090616

double nH2_3D(double, double, double, double);  //ECC20141024
double nHI_3D(double, double, double);          //ECC20141024

//double nH2_3D(double, double, double, double);//ECC20141024
double nHII(double, double);
double nHII3D(const double x, const double y, const double z);//AWS20131004

double interp_HI_gas_cube(double x, double y, double z, GasCubeData gcd);
double interp_HI_gas_cube_rlb(double r, double l, double b, GasCubeData gcd);
double interp_CO_gas_cube(double x, double y, double z, double H2toCO, GasCubeData gcd);
double interp_CO_gas_cube_rlb(double r, double l, double b, double H2toCO, GasCubeData gcd);

int get_gas_cube_coords(int, int, int, int, int, int);//ECC20140824

int    nH_set_model(int nHI_model_, int nH2_model_, int nHII_model_, int debug_); //AWS20090814

//double nHI_av (double,double,double,double,double); //IMOS20080114
//double nH2_av (double,double,double,double,double); //IMOS20080114
//double nHII_av(double,double,double,double,double); //IMOS20080114

double nH_av(double, double, double, double, double, double(*)(double, double)); //IMOS20080114
double nH_av(double, double, double, double, double, double, double(*)(double, double, double));                      //AWS20090616
double nH_av_lb(double, double, double, double, double, double, double, double(*)(double, double)); //IMOS20080114
double nH_av_lb(double, double, double, double, double, double, double, double,double(*)(double, double, double)); //AWS20090616

/* ECC20141024:  3d versions with/without GasCubeData and with/without H2toCO */
double nH_av_3D(double, double, double, double, double, double(*)(double, double, double));                         //ECC20141024
double nH_av_3D(double, double, double, double, double,double, double(*)(double, double, double, double));          //ECC20141024
double nH_av_3D(double, double, double, double, double,GasCubeData, double(*)(double, double, double,GasCubeData)); //ECC20141024
double nH_av_3D(double, double, double, double, double,double, GasCubeData, double(*)(double, double, double, double, GasCubeData)); //ECC20141024

double nH_av_lb_3D(double, double, double, double, double, double, double, double(*)(double, double, double)); //ECC20141024
double nH_av_lb_3D(double, double, double, double, double, double, double, double, double(*)(double, double, double, double)); //ECC20141024
double nH_av_lb_3D(double, double, double, double, double, double, double, GasCubeData, double(*)(double, double, double, GasCubeData)); //ECC20141024
double nH_av_lb_3D(double, double, double, double, double, double, double, double, GasCubeData, double(*)(double, double, double, double, GasCubeData)); //ECC20141024
double nH_av_lb_3D_rlb(double, double, double, double, double, double, double, double,GasCubeData,double(*)(double,double,double,double,GasCubeData)); //ECC20141024
double nH_av_lb_3D_rlb(double, double, double, double, double, double, double, GasCubeData, double(*)(double, double, double, GasCubeData)); //ECC20141024

/* ECC20141024:  Precomputed angles for full 3D F07 model */
const double SIN_THETA_C = 0.939693, COS_THETA_C = 0.34202;
const double SIN_ALPHA = 0.233445  , COS_ALPHA = 0.97237;
const double SIN_BETA  = 0.34202   , COS_BETA = 0.939693;
const double SIN_THETA_D = 0.748956, COS_THETA_D = 0.66262;
const double PHI_00 = -3.84229569333e+14; // gravitational potential at zero


double B_field_model(double, double, int);
double B_field_model(double, double, double, int);

double gauss(double mean, double sigma);

int test_sigma_boron_dec_heinbach_simon();
int test_kinematic();
int test_He_to_H_CS();
int test_nH(); 
int test_Distribution();
int test_float_accuracy();

#endif
