#include "B_field_3D_model.h"
#include "synchrotron_emissivity.h"
#include "synchrotron_emissivity_B_field.h"

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <string>

using namespace std;
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double synchrotron_emissivity_aws(double gamma, double nu, double x, double y, double z,  const std::string &name, const std::vector<double> &parameters, int debug=0)
{

 int options;
  double x0,y0,z0;

 double synch_emissivity_total;
 double synch_emissivity_reg;
 double synch_emissivity_par;
 double synch_emissivity_perp;
 double synch_emissivity_random;


  x0=8.5; // solar position
  y0=0.;
  z0=0.;

  options=0;


 synch_emissivity_total  =
 synchrotron_emissivity_B_field(  gamma,  nu, 
				  name, parameters, x, y, z, options,
                                  x0, y0, z0,
                                  synch_emissivity_reg,synch_emissivity_par,synch_emissivity_perp,synch_emissivity_random,
			      	  debug );

 if(debug==1)
 {
 cout<<"-------------------  synchrotron_emissivity_aws: synchrotron_emissivity_B_field at ";
 cout<<"(x, y, z) = ("<<x<<", "<<y<<", "<<z<<")  "<<endl;
 cout<<"B field name="<<name<<"  ";

 cout
  <<" gamma= "<<gamma<<" nu="<<nu << " Hz" <<endl
  <<" synch emissivity random        =  "  << synch_emissivity_random <<endl
  <<" synch emissivity regular       =  "  << synch_emissivity_reg    <<endl
  <<" synch emissivity  pol parallel =   " << synch_emissivity_par    <<endl
  <<" synch emissivity  pol perp     =   " << synch_emissivity_perp   <<endl
  <<" synch emissivity total         =  "  << synch_emissivity_total  <<endl
  <<endl;

}



  return synch_emissivity_total;
}


//////////////////////////////////////////////////////////////////////////////////////////////////

double synchrotron_emissivity_aws(double gamma, double nu, double R, double z, const std::string &name, const std::vector<double> &parameters, int debug=0)

{

 double synch_emissivity_total;
 double x,y;
 double phi;
 const double pi=acos(-1.0);
 phi=pi/4.;             // arbitrary angle, should really take azimuthal average
 x=R*cos(phi); 
 y=R*sin(phi);
 synch_emissivity_total  =  synchrotron_emissivity_aws(gamma,nu,x,y,z,  name, parameters, debug);

 if(debug==1)
 {
 cout<<"-------------------  synchrotron_emissivity_aws: synchrotron_emissivity_B_field at ";
 cout<<"(R, x, y, z) = ("<<R<<","<<x<<", "<<y<<", "<<z<<")  "
     <<" synch emissivity total         =  "  << synch_emissivity_total  <<endl;
 }

 return synch_emissivity_total;
}



// Stokes parameters version AWS20100707
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double synchrotron_emissivity_aws(double gamma, double nu, double x, double y, double z,  const std::string &name, const std::vector<double> &parameters, double &I, double &Q, double &U, int debug=0)
{

 int options;
  double x0,y0,z0;

 double synch_emissivity_total;
 double synch_emissivity_reg;
 double synch_emissivity_par;
 double synch_emissivity_perp;
 double synch_emissivity_random;


  x0=8.5; // solar position
  y0=0.;
  z0=0.;

  options=0;


 synch_emissivity_total  =
 synchrotron_emissivity_B_field(  gamma,  nu, 
				  name, parameters, x, y, z, options,
                                  x0, y0, z0,
                                  synch_emissivity_reg,synch_emissivity_par,synch_emissivity_perp,synch_emissivity_random,
                                  I, Q, U,
			      	  debug );

 if(debug==1)
 {
 cout<<"-------------------  synchrotron_emissivity_aws: synchrotron_emissivity_B_field at ";
 cout<<"(x, y, z) = ("<<x<<", "<<y<<", "<<z<<")  "<<endl;
 cout<<"B field name="<<name<<"  ";

 cout
  <<" gamma= "<<gamma<<" nu="<<nu << " Hz" <<endl
  <<" synch emissivity random        =  "  << synch_emissivity_random <<endl
  <<" synch emissivity regular       =  "  << synch_emissivity_reg    <<endl
  <<" synch emissivity  pol parallel =   " << synch_emissivity_par    <<endl
  <<" synch emissivity  pol perp     =   " << synch_emissivity_perp   <<endl
  <<" synch emissivity Stokes I      =  "  << I                       <<endl
  <<" synch emissivity Stokes Q      =  "  << Q                       <<endl
  <<" synch emissivity Stokes U      =  "  << U                       <<endl
  <<endl;

}



  return synch_emissivity_total;
}

//////////////////////////////////////////////////////////////////////////////////////////////////

double synchrotron_emissivity_aws(double gamma, double nu, double R, double z, const std::string &name, const std::vector<double> &parameters, double &I, double &Q, double &U, int debug=0)

{

 double synch_emissivity_total;
 double x,y;
 double phi;
 double pi=acos(-1.0);
 phi=pi/4.;             // arbitrary angle, should really take azimuthal average
 x=R*cos(phi); 
 y=R*sin(phi);
 synch_emissivity_total  =  synchrotron_emissivity_aws(gamma,nu,x,y,z,  name, parameters, I, Q, U,  debug);

 if(debug==1)
 {
 cout<<"-------------------  synchrotron_emissivity_aws: synchrotron_emissivity_B_field at ";
 cout<<"(R, x, y, z) = ("<<R<<","<<x<<", "<<y<<", "<<z<<")  "
     <<" synch emissivity total         =  "  << synch_emissivity_total  <<endl;
 }

 return synch_emissivity_total;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// following routines using B_field_model were obsolete so removed AWS2011215
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


