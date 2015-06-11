#include "galprop_internal.h"
#include "kappa_free_free.h"



//Convenience functions to split up the emission and absorption calculation
inline double gff(const double nu, const double Te) {
  return 10.6 + 1.9*log10(Te) - 1.26*log10(nu);           // Allen Astrophysical Quantities 4th Edition Section 5.9
}

inline double kff(const double nu, const double Ne, const double Te) {
  return 0.0178*gff(nu,Te)*Ne*Ne/(nu*nu*Te*sqrt(Te));// Allen Astrophysical Quantities 4th Edition Section 5.9
}

inline double eff(const double nu, const double Ne, const double Te) {
  return 5.444e-39*gff(nu,Te)*Ne*Ne/sqrt(Te); // Allen Astrophysical Quantities 4th Edition Section 5.9
}


///////////////////////////////////////////////////////////
double kappa_free_free(double nu,double Ne, double Te, double &emiss_free_free, int options, int debug)
{ 
  // AWS20110630
  // radio free-free absorption coefficient, in cm^-1
  // nu = frequency in Hz, Ne=electron density in cm^-3, Te=electron temperature in K
  // formula from Peterson and Webber, ApJ 575,217 (2002), equations 3 and 4.
  // quoting Allen's Astrophysical Quantities. Now use 4th Edition.
  // checked against Spitzer Physical Processes in the Interstellar Medium eq 3-55 and 3-57
  // where Gaunt factor gff is given explicitly using ln and log10, and agrees reasonably
  // Allen formula for Gaunt factor includes Z multiplying nu, explains why tabulated values smaller (Z>1 presumably).
  // and kappa~Z^2, but included in Ne^2 presumably. Here using Z=1 but what about He etc ?
  // Webber and Higbie JGR 113, A11106 (2008) quote Lang Astrophysical Formulae (1999) for WIM so this is an option also.
  // The comparison below gives Lang < Allen so needs following up.
  
  // The free-free emissivity (erg s-1 sr-1 Hz-1) from Allen is also return since it is closely related to absorption
  // and uses the same Gaunt factor. Since it is used for radio frequencies the Planck term is omitted.
  // Again Z=1 is used.

  const double gff_ = gff(nu,Te);//10.6 + 1.9*log10(Te) - 1.26*log10(nu);           // Allen Astrophysical Quantities 4th Edition Section 5.9
  

  const double kappa_free_free_= kff(nu,Ne,Te);//   0.0178*gff*Ne*Ne/(nu*nu*Te*sqrt(Te));// Allen Astrophysical Quantities 4th Edition Section 5.9
               emiss_free_free = eff(nu,Ne,Te);//5.444e-39*gff*Ne*Ne/          sqrt(Te); // Allen Astrophysical Quantities 4th Edition Section 5.9

  // Lang Astrophysical Formulae 3rd Edition equation 1.223 gives this approximation from Altenhoff etal 1960:
  // Lang formula is for tau (kappa*pathlength) with  EM in cm-6 pc so convert to cm, and nu in GHz
  const double kappa_free_free_lang= 8.235e-2 * pow(Te,-1.35) * Ne * Ne* pow(nu*1.0e9,-2.1)*3.08e18;




  if(debug==1)
  {
   cout<<" kappa_free_free:"<<endl;
   cout<<" nu="<<nu<<endl;
   cout<<" Ne="<<Ne<<endl;
   cout<<" Te="<<Te<<endl;
   cout<<" gff="<<gff_<<endl;
   cout<<" kappa_free_free="<<kappa_free_free_<<endl;
   cout<<" kappa_free_free_lang="<<kappa_free_free_lang<<endl;
   cout<<" emiss_free_free="<<emiss_free_free<<endl;
  }


  return  kappa_free_free_;
}
///////////////////////////////////////////////////////////////

double kappa_free_free_2D(double nu, double R, double z,  double Te, double clumping_factor, double &emiss_free_free, int options, int debug)
{
  // AWS20110630
  // free-free absorption coefficient in cm^-1 at position in Galaxy
  // based on model for free-electron distribution

 
  double Ne = nHII  (R  ,z);

  int debug_=0;

  double  kappa_free_free_ =   kappa_free_free( nu, Ne,  Te, emiss_free_free, 0,debug_);
          kappa_free_free_ *= clumping_factor;//AWS20110704
          emiss_free_free  *= clumping_factor;//AWS20110704

  if(debug==1)
  {
   cout<<" kappa_free_free_2D for free electron distribution at:"<<endl;
   cout<<" R="<<R<<endl;
   cout<<" z="<<z<<endl;
   cout<<" nu="<<nu<<endl;
   cout<<" Ne="<<Ne<<endl;
   cout<<" Te="<<Te<<endl;
   cout<<" kappa_free_free at R,z ="<<kappa_free_free_<<endl;
   cout<<" emiss_free_free at R,z ="<<emiss_free_free <<endl;
  }

return kappa_free_free_;
}
///////////////////////////////////////////////////////////////

double kappa_free_free_3D(double nu, double x, double y, double z,  double Te,  double clumping_factor, double &emiss_free_free, int options, int debug)
{
  // AWS20110630
  // free-free absorption coefficient in cm^-1 at position in Galaxy
  // based on model for free-electron distribution
  

  double Ne = nHII3D(x,y,z);

  int debug_=0;

  double  kappa_free_free_ =   kappa_free_free( nu, Ne,  Te, emiss_free_free, 0,debug_);
          kappa_free_free_ *= clumping_factor;//AWS20110704
          emiss_free_free  *= clumping_factor;//AWS20110704

  if(debug==1)
  {
   cout<<" kappa_free_free_3D for free electron distribution at:"<<endl;
   cout<<" x="<<x<<endl;
   cout<<" y="<<y<<endl;
   cout<<" z="<<z<<endl;
   cout<<" nu="<<nu<<endl;
   cout<<" Ne="<<Ne<<endl;
   cout<<" Te="<<Te<<endl;
   cout<<" kappa_free_free at xyz ="<<kappa_free_free_<<endl;
   cout<<" emiss_free_free at xyz ="<<emiss_free_free <<endl;
  }

return kappa_free_free_;
}

