


//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|
//  AWS20140507
// source function for  p+p -> d + pi+ fusion reaction
// this source function is added to existing deuterium source function (e.g. spallation)
//
// Started with copy of gen_secondary_positron_source.cc for template.
// 
// I thank Nicolas Picot-Clemente for suggesting this process, and help in implementation and testing.
//
// CR density gcr.cr_density is in c/4pi * n(p) [cm s^-1 sr^-1 cm^-3 MeV^-1]
//
// 
//
// The CR source function [cm^-2 s^-2 sr^-1 MeV^-1] as used in
// galprop is defined as  (c/4pi * q)  [q = cm^-3 s^-1 MeV^-1]:
//                 ___
//         c       \    /                   c   d sigma_ij(p,p')
// q(p) * --- = c  /__  \ dp' beta n_j(p') ---  --------------- ,
//        4pi           /                  4pi        dp
// 
// where n_i is the gas density, d sigma_ij(p,p')/dp is
// the production cross section, n_j(p') is the CR species density, 
// and p' is the total momentum of a nucleus.
// Substitution of dp' with d(log Ekin) gives:
//                                     ___
//       c                /            \               c  d sigma_ij(p,Ekin)
// q(p)*--- = c      d(log Ekin)  Ekin /__  n_j(Ekin) --- -----------------
//      4pi                /                    j             4pi       dp  
//                         ___     ___      ___
//                         \       \        \              c   d sigma_ij(p,Ekin)
//      = c A /\(log Ekin) /__ n_i /__ Ekin /__ n_j(Ekin) ---  -----------------,
//                        i=H,He   Ekin      j            4pi        dp             
// 
// where /\=Delta, and we used dp'=1/beta A Ekin d(log Ekin).
// 
// 
// 
//=="====!===="====!===="====!===="====!===="====!===="====!===="====!===="====!
using namespace std;

#include"galprop_classes.h"
#include"galprop_internal.h"

#include <fort_interface.h>

#include <string>
#include <cstring>
#include <sstream>

#include <ErrorLogger.h>

double cs_ppd_fusion(double Ekin_d, double Ekin_p,int cs_init, int debug); // local to this file
double cs_ppd_total (               double Ekin_p,int cs_init, int debug); // local to this file

/////////////////////////////////////////////////////////////////////////////////

int Galprop::gen_secondary_deuterium_fusion_source(Particle& particle) 
{

  INFO("Entry");

  cout<<">>gen_secondary_deuterium_fusion_source"<<endl;
  
  
  ostringstream buf;
  buf << "Generating " << particle.name << " source function for spatial dimensions " << gcr[0].n_spatial_dimensions;
  INFO(buf.str());

  

  const string particleName = particle.name;

  int debug=0;
  if(galdef.verbose==-5000) debug=1;// selectable debug

  int stat=0, iprotons=-1,  i;
  double cs_ppd_fusion_; // cross section pp->d + pi+
  Distribution protons;              
  
  // identify CR protons                 
  if (2 == galdef.n_spatial_dimensions) 
    protons.init(gcr[0].n_rgrid, gcr[0].n_zgrid, gcr[0].n_pgrid);
  
  if (3 == galdef.n_spatial_dimensions) 
    protons.init(gcr[0].n_xgrid, gcr[0].n_ygrid, gcr[0].n_zgrid, gcr[0].n_pgrid);

  protons=0.;

  // sum protons primary and secondary etc 

  for(i=0; i<n_species; i++)  
    if(101==100*gcr[i].Z+gcr[i].A)
      {
	iprotons=i;
	protons+=gcr[iprotons].cr_density;
	buf.str("");
	buf<<"  CR protons found as species #"<<iprotons;
	INFO(buf.str());
      }
  if(iprotons==-1) { cout<<"  CR protons not found!"<<endl;  return 1; }
  
  int cs_init=1;
  cs_ppd_fusion(0., 0., cs_init, debug ); // initialize data table
  cs_init=0;

  // need to apply factor before additing to existing source function

  double factor=1.e-24 *1.e-3 *C *log(galdef.Ekin_factor); // transformation to cm^2/MeV and logarithmic integration factor
  
  for(int ip_sec=0; ip_sec<particle.n_pgrid; ip_sec++) 
  {

   for(int ip=0; ip<gcr[iprotons].n_pgrid; ip++) 
   {

      //ostringstream buf;
      //buf << "Generating for secondary energy " << ip_sec << " " << particle.Etot[ip_sec] << " primary energy " << ip << " " << gcr[iprotons].p[ip];
      //INFO(buf.str());

    
    
    
     cs_ppd_fusion_=cs_ppd_fusion(particle.Ekin[ip_sec]*1.e-3, gcr[iprotons].Ekin[ip]*1.e-3,cs_init, debug ); // barn / GeV
      
      if(galaxy.n_spatial_dimensions==2)
	for(int ir=0; ir<gcr[iprotons].n_rgrid; ir++)
	for(int iz=0; iz<gcr[iprotons].n_zgrid; iz++)
	{
	      particle.secondary_source_function.d2[ir][iz].s[ip_sec ]+=                  
		(galaxy.n_HI.d2[ir][iz].s[0] + 2.0*galaxy.n_H2.d2[ir][iz].s[0] + galaxy.n_HII.d2[ir][iz].s[0])
		 *cs_ppd_fusion_  
                 *protons.d2[ir][iz].s[ip] *gcr[iprotons].Ekin[ip]
                 *factor;
    
        }  //  iz //  ir //  particle.n_spatial_dimensions==2
 


      if(galaxy.n_spatial_dimensions==3)
            for(int ix=0; ix<gcr[iprotons].n_xgrid; ix++)
            for(int iy=0; iy<gcr[iprotons].n_ygrid; iy++)
            for(int iz=0; iz<gcr[iprotons].n_zgrid; iz++)
	    {
                  particle.secondary_source_function.d3[ix][iy][iz].s[ip_sec ]+=         
                  (galaxy.n_HI.d3[ix][iy][iz].s[0] + 2.0*galaxy.n_H2.d3[ix][iy][iz].s[0] + galaxy.n_HII.d3[ix][iy][iz].s[0])
                   *cs_ppd_fusion_   
                   *protons.d3[ix][iy][iz].s[ip] *gcr[iprotons].Ekin[ip]
		   *factor;

             }  //  iz //  iy  //  ix  //  particle.n_spatial_dimensions==3


      }  //  ip
   }  //  ip_sec
 
 
  

   protons.delete_array();                        

   if(galdef.verbose>=2)
   {
      cout<<"   particle.secondary_source_function for "<<particle.name<<endl;
      particle.secondary_source_function.print();
   }
   
   cout<<" << gen_secondary_deuterium_fusion_source"<<endl;
   
   INFO("Exit");
   return stat;
}

///////////////////////////////////////////////////////////////
double cs_ppd_fusion(double Ekin_d, double Ekin_p, int cs_init, int debug)
{
  // differential cross-section for fusion reaction  pp->d + pi+, barn/GeV

  // based on Meyer, J.P., A&A Supp 7, 417 (1972) http://adsabs.harvard.edu/abs/1972A%26AS....7..417M
  // and Coste, B., etal,  A&A    539, A88 (2012) http://adsabs.harvard.edu/abs/2012A%26A...539A..88C

  // cs_ppd_fusion = differential cross-section per momentum of deuterium, barn/GeV
  // Ekin_d = deuterium kinetic energy per nucleon in GeV
  // Ekin_p = proton    kinetic energy             in GeV
  // cs_init  1= initialize data, 0= assume data initialized

  // using GeV and barn/GeV, and order of arguments here for consistency with other gen_secondary.. cross-section routines in GALPROP (pp_meson, antiprotons etc)

  double result=0;

  double cs_ppd_total_=0.;
   
  cs_ppd_total_=cs_ppd_total(Ekin_p,cs_init,debug); // barn
  

  // kinematics from Meyer (1972), Section 3.1

  // particle rest masses from constants.h 
  double M1=Mp;   // proton 
  double M2=Mp;   // proton
  double M3=Md;   // deuteron
  double M4=Mpi1; // pi+

  
  double gamma1 = (M1+Ekin_p)/M1;              // gamma of proton in LAB system
  double  beta1 = sqrt(1.-1/(gamma1*gamma1));  // beta  of proton in LAB system

  double  beta_CM = gamma1*beta1/(gamma1 +1. );      //  beta of CM system (M1=M2)
  double gamma_CM = sqrt(1./(1. -beta_CM*beta_CM));  // gamma of CM system

  double Wt       = M1 + M2 + Ekin_p;   // total energy in LAB system
  double Wt_prime = Wt/gamma_CM;        // total energy in  CM system

  // if Wt_prime < M3 + M4 the reaction is energetically excluded
  if(Wt_prime<M3+M4 && debug==0){ result=0.;return result;}
    //{ cout<< " NB reaction energetically excluded!"; cout<<endl; result=0.;return result;}

  double gamma3_prime = (Wt_prime*Wt_prime + M3*M3 - M4*M4)/(M3 * 2.*Wt_prime); // gamma of deuteron in CM system
  double  beta3_prime = sqrt(1.-1/(gamma3_prime*gamma3_prime));                 //  beta of deuteron in CM system

  double A = M3*(gamma_CM * gamma3_prime - 1.0); 
  double B = M3* gamma_CM * gamma3_prime * beta_CM * beta3_prime;

 
  double E3_min = A - B; // deuteron minimum Ekin for given Ekin_p
  double E3_max = A + B; // deuteron maximum Ekin for given Ekin_p
 
  double E3    = 2 * Ekin_d;        // total kinetic energy of deuteron (two nucleons)
  double ET3   = E3 + M3;           //       total   energy of deuteron
  double p3    = sqrt(ET3*ET3-M3*M3); //             momentum of deuteron
  double beta3 = p3/ET3;             //             beta     of deuteron

  double costheta3_prime = (E3 - A)/B;

 
  double dsigmadOmegaprime=0;

  if(costheta3_prime>-1. && costheta3_prime<+1.)
  { 
  // dsigma/dE3 = (2pi/B) dsigma/dOmega', where E3 is total kinetic energy of deuteron
  // dsigma/dEn    4pi/B   dsigma/dOmega' for              (kinetic energy / nucleon)^-1, En = E3 / 2
  // dsigma/dp = (2pi/B) dsigma/dOmega' dE3/dp3 =  (2pi/B) dsigma/dOmega' p3/E3 =  (2pi/B) dsigma/dOmega' beta3
  // dsigma/dOmega'               ~ 0.22 + cos^2 (theta')  approximately (Meyer 1972), need to normalize integral to total sigma.
  // dsigma/dOmega' = sigma_total ( 0.22 + cos^2 (theta')) /( 4pi (0.22 + 1/3))
  //               using integral  (0.22 + cos^2 (theta') dOmega= integral  ( 0.22 + cos^2 (theta') 2pi sin theta dtheta = 4pi (0.22 + 1/3)
  // dsigma/dEn = sigma_total (0.22 + costheta3_prime * costheta3_prime)/(0.22+1/3) /B
  // dsigma/dp  = sigma_total (0.22 + costheta3_prime * costheta3_prime)/(0.22+1/3) /2B *  beta3

    if(debug==1)
    {
     dsigmadOmegaprime  = 0.22 + costheta3_prime * costheta3_prime; // used only for debug
     dsigmadOmegaprime /= 4.*Pi*(0.22+1./3);
    }

   
    //     result = cs_ppd_total_ *  (0.22 + costheta3_prime * costheta3_prime) / (0.22+1./3.) /B ;              //AWS20140508 this is for dsigma/dEkin

           result = cs_ppd_total_ *  (0.22 + costheta3_prime * costheta3_prime) / (0.22+1./3.) /(2.*B) * beta3 ; //AWS20140508 this is for dsigma/dp as required
   }

  if(debug==1)
  {
   cout<<endl<<"cs_ppd_fusion: debug"<<endl;
   cout<<"Ekin_p="<<Ekin_p<<" GeV Ekin_d="<<Ekin_d<<" GeV/nucleon"<<endl;
   cout<<"proton in LAB: gamma1 ="<<gamma1 <<" beta1 ="<<beta1 <<endl;
   cout<<"CM: gamma_CM ="<<gamma_CM <<" beta_CM ="<<beta_CM <<endl;
   cout<<"Wt ="<<Wt <<" Wt_prime ="<<Wt_prime <<endl;
   cout<<" M1+M2="<< M1+M2 <<" M3+M4 ="<<M3+M4;
   if(Wt_prime<M3+M4)
    cout<< " reaction energetically excluded!"; cout<<endl;
   cout<<"deuteron in CM: gamma3_prime ="<<gamma3_prime <<" beta3_prime="<<beta3_prime <<endl;
   cout<<"A ="<<A <<" B="<< B <<endl;
   cout<<"deuteron total Ekin range E3_min ="<<E3_min <<" E3_max="<< E3_max <<endl;
   cout<<" E3="<<E3 <<" costheta3_prime="<< costheta3_prime;
   if(!(costheta3_prime>-1. && costheta3_prime<+1.))
    cout<< " cos theta out of range!";cout<<endl; 
   cout<<" ET3="<<ET3 <<" p3="<<p3 <<" beta3="<< beta3;
   cout<<" dsigmadOmegaprime="<<dsigmadOmegaprime <<"  cs_ppd_total_ ="<< cs_ppd_total_ <<endl;
   cout<<"cs_ppd_fusion:"<< " Ekin_p="<<Ekin_p<<" GeV Ekin_d="<<Ekin_d<<" result="<< result<<" barn / GeV" <<endl<<endl;
  }

  return result;
}
/////////////////////////////////////////////
double cs_ppd_total(double Ekin_p, int cs_init, int debug)
{

  //  pp-> d + pi+      total cross-section in barn. Ekin_p in GeV. 

  /*




  */
  // static to keep after initialization
  static int n;
  static valarray<double> E,cs;


  int i;

  if (cs_init==1)
  {
    cout<<" gen_secondary_deuterium_fusion_source: initializing total cross-section data"<<endl;
   n=29;
   E .resize(n);
   cs.resize(n);

  // from Nicolas Picot-Clemente, digitalized from  Coste etal 2012 Fig B.3

   i=0;


  //     proton Ekin, MeV         mbarn
   E[i]=     313.073   ;cs[i]=    0.0164899 ; i++;
   E[i]=     313.186   ;cs[i]=    0.0264723 ; i++;
   E[i]=     315.180   ;cs[i]=    0.0655127 ; i++;  
   E[i]=     322.109   ;cs[i]=    0.118792  ; i++;       
   E[i]=     331.693   ;cs[i]=    0.171161  ; i++;       
   E[i]=     339.002   ;cs[i]=    0.290065  ; i++;      
   E[i]=     354.283   ;cs[i]=    0.435262  ; i++;      
   E[i]=     378.575   ;cs[i]=    0.627206  ; i++;      
   E[i]=     416.661   ;cs[i]=    0.980264  ; i++;       
   E[i]=     479.372   ;cs[i]=    1.63943   ; i++;      
   E[i]=     543.459   ;cs[i]=    2.49409   ; i++;      
   E[i]=     598.225   ;cs[i]=    3.18231   ; i++;      
   E[i]=     644.211   ;cs[i]=    2.93483   ; i++;       
   E[i]=     698.922   ;cs[i]=    2.52967   ; i++;     
   E[i]=     752.743   ;cs[i]=    1.98346   ; i++;     
   E[i]=     810.742   ;cs[i]=    1.47329   ; i++;       
   E[i]=     873.246   ;cs[i]=    1.03671   ; i++;   
   E[i]=     947.584   ;cs[i]=    0.700512  ; i++;      
   E[i]=    1028.23    ;cs[i]=    0.486318  ; i++;      
   E[i]=    1115.67    ;cs[i]=    0.366156  ; i++;      
   E[i]=    1219.57    ;cs[i]=    0.264728  ; i++;   
   E[i]=    1373.30    ;cs[i]=    0.174121  ; i++;    
   E[i]=    1592.92    ;cs[i]=    0.109981  ; i++;      
   E[i]=    1861.32    ;cs[i]=    0.0723451 ; i++;     
   E[i]=    2111.31    ;cs[i]=    0.0537441 ; i++; 
   E[i]=    2466.99    ;cs[i]=    0.0368167 ; i++; 
   E[i]=    2777.70    ;cs[i]=    0.0273500 ; i++;
   E[i]=    3245.59    ;cs[i]=    0.0189909 ; i++; 
   E[i]=    4054.18    ;cs[i]=    0.0104799 ;  
 
   E /=1000.; // MeV -> GeV
   cs/=1000.; // mbarn - barn

   if (debug==1) cout<<"last i="<<i<<endl;
   if(i!=n-1){cout<<"error in cs_ppd_fusion_total data"<<endl; exit(1);}

  } // cs_init==1

   double result=0;

   if(Ekin_p>E[0] && Ekin_p<E[n-2])
   {    
    for(i=1;i<n-1;i++)   if(Ekin_p < E[i])  break;

    result=  cs[i] +  (cs[i+1] - cs[i]) * (Ekin_p - E[i])/ (E[i+1] - E[i]); // linear interpolation
   
    if (debug==1){cout<<"cs_ppd_fusion_total: "<<" Ekin_p="<<Ekin_p<<" i="<<i<< " E[i]="<<E[i]<<" E[i+1]="<<E[i+1]
			<<" cs[i]="<<cs[i]<<" cs[i+1]="<<cs[i+1] <<" result="<<result <<endl;}
   
   }




  if(debug==1) cout<<"cs_ppd_total: Ekin_p="<<Ekin_p<<" GeV  result (barn)="<<result<<endl;

  return result;
}
/////////////////////////////////////
/*
int main()
{
  cout<<" gen_secondary_deuterium_fusion_source test programme"<<endl;
  double Ekin_d, Ekin_p; int cs_init,debug;
 
  cs_init=1;
  debug=1;
  cs_ppd_fusion(0.,0., cs_init, debug);

  cs_init=0;
  
  for(Ekin_p=0.2;Ekin_p<5.0;Ekin_p*=1.2)
  for(Ekin_d=0.2;Ekin_d<5.0;Ekin_d*=1.2)
  { 

    cout<<"Ekin_p="<<Ekin_p<<" GeV Ekin_d="<<Ekin_d<<" GeV/nucleon"
       <<" cs_ppd_fusion="<< cs_ppd_fusion(Ekin_d, Ekin_p, cs_init, debug)<<" barn / GeV"<<endl;

  }

  debug=0;
  for(Ekin_p=0.2;Ekin_p<5.0;Ekin_p*=1.2)
  for(Ekin_d=0.2;Ekin_d<5.0;Ekin_d*=1.2)
  { 

    cout<<"Ekin_p="<<Ekin_p<<" GeV Ekin_d="<<Ekin_d<<" GeV/nucleon"
       <<" cs_ppd_fusion="<< cs_ppd_fusion(Ekin_d, Ekin_p, cs_init, debug)<<" barn / GeV"<<endl;

  }

    return 0;
}
*/ 


/////////////////////////////////////////////////////////
/*
sample compilation of standalone test programme 

 g++ -DHAVE_CONFIG_H -I. -I.. -I.. -I../include -I/afs/ipp-garching.mpg.de/home/a/aws/propagate/c/cfitsio/3.26/gcc_sles11_olga2/cfitsio/include -I/afs/ipp-garching.mpg.de/home/a/aws/propagate/c/CCfits/gcc_olga2_sles11_cfitsio3.26/usr/local/include -I/afs/ipp-garching.mpg.de/home/a/aws/Healpix/Healpix_3.11/src/cxx/generic_gcc/include -I/afs/ipp-garching.mpg.de/home/a/aws/propagate/c/CLHEP/2.0.4.3/installed_here/include -I/afs/ipp-garching.mpg.de/home/a/aws/gsl/gsl-1.10/olga/include -DGALDEF_PATH=\"/afs/ipp/u/aws/propagate/c/GALDEF\" -DFITSDATA_PATH=\"/afs/ipp/u/aws/propagate/c/FITS\" -DDATA_PATH=\"/afs/ipp/u/aws/propagate/c/galprop/fortran_data\" -g -pipe -O3 -Wall -MT gen_secondary_deuterium_fusion_source.o -MD -MP -MF .deps/gen_secondary_deuterium_fusion_source.Tpo gen_secondary_deuterium_fusion_source.cc -c

g++ -g -pipe -O3 -Wall gen_secondary_deuterium_fusion_source.o libgalprop.a libskymap.a -L/usr/lib64/gcc/x86_64-suse-linux/4.3 -L/usr/lib64/gcc/x86_64-suse-linux/4.3/../../../../lib64 -L/lib/../lib64 -L/usr/lib/../lib64 -L/usr/lib64/gcc/x86_64-suse-linux/4.3/../../../../x86_64-suse-linux/lib -L/usr/lib64/gcc/x86_64-suse-linux/4.3/../../.. -L/afs/ipp-garching.mpg.de/home/a/aws/propagate/c/cfitsio/3.26/gcc_sles11_olga2/cfitsio/lib -L/afs/ipp-garching.mpg.de/home/a/aws/propagate/c/CCfits/gcc_olga2_sles11_cfitsio3.26/usr/local/lib /afs/ipp-garching.mpg.de/home/a/aws/propagate/c/CCfits/gcc_olga2_sles11_cfitsio3.26/usr/local/lib/libCCfits.so -lcfitsio -L/afs/ipp-garching.mpg.de/home/a/aws/Healpix/Healpix_3.11/src/cxx/generic_gcc/lib -lhealpix_cxx -lcxxsupport -lfftpack -lgfortranbegin -lgfortran -L/afs/ipp-garching.mpg.de/home/a/aws/gsl/gsl-1.10/olga/lib /afs/ipp-garching.mpg.de/home/a/aws/gsl/gsl-1.10/olga/lib/libgsl.so /afs/ipp-garching.mpg.de/home/a/aws/gsl/gsl-1.10/olga/lib/libgslcblas.so -lm -L/afs/ipp-garching.mpg.de/home/a/aws/propagate/c/CLHEP/2.0.4.3/installed_here/lib -lCLHEP -Wl,-rpath -Wl,/afs/ipp-garching.mpg.de/home/a/aws/propagate/c/CCfits/gcc_olga2_sles11_cfitsio3.26/usr/local/lib -Wl,-rpath -Wl,/afs/ipp-garching.mpg.de/home/a/aws/gsl/gsl-1.10/olga/lib -Wl,-rpath -Wl,/afs/ipp-garching.mpg.de/home/a/aws/propagate/c/CCfits/gcc_olga2_sles11_cfitsio3.26/usr/local/lib -Wl,-rpath -Wl,/afs/ipp-garching.mpg.de/home/a/aws/gsl/gsl-1.10/olga/lib

./a.out


*/
