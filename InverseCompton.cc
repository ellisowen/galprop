/* general-purpose Inverse Compton emissivity routines.
For isotropic and anisotropic radiation fields.
For usage see InverseCompton::test() examples, and the main program at end of this file.
A. Strong, 20131204
History:
20120104: first version for Stellarics package.
20131204: improved test case
*/

using namespace std;
#include<iostream>
#include<cmath>

#include"InverseCompton.h"

/////////////////////////////////////////////////////////
void InverseCompton:: init(int debug_)
{

  cout<<">> InverseCompton:: init()"<<endl;

#include"Constants.h"



  constant=pi*re*re*me*me;
  initialized=1;
  debug=debug_;

  cout<<"<< InverseCompton:: init()"<<endl;
}
/////////////////////////////////////////////////////////

double InverseCompton::IC_sigma_KN(double Ee, double Eph, double eta, double Egamma)
{
  // Inverse Compton differential cross-section using Klein-Nishina cross-section.

  // units: cm^2 MeV^-1
  // uses equations (6)+(7) of    Orlando & Strong A&A 480, 847 (2008)
  // and  equation  (9)     of Moskalenko & Strong ApJ 528, 357 (2000)
  // stand-alone form, all constants defined here
  // energies in MeV
  // Ee     = electron energy
  // Eph    = target photon energy
  // eta    = angle between electron and target photon in radians
  // Egamma = gamma ray energy

#include"Constants.h"

  double result;
  if(Egamma>=Ee) {result=0.; return result;}
  if(initialized!=1) init(0); // just in case


  double g    = Ee/me;                      // electron Lorentz factor
  double beta = sqrt(1.0-1.0/(g*g));        // electron velocity/c
  double v    = Egamma/Ee;                  // fractional energy transfer
  double Ephd = g*Eph*(1.0+beta*cos(eta));  // energy of target photon in electron frame

  double Egamma_max = 2.0*g*Ephd/(1.0+2.0*Ephd/me); //NB MS2000 eq (9) Ephd has units of me
  if(Egamma>=Egamma_max) {result=0.; return result;}

  
 
  result= constant/(Eph*Ee*Ee)
         *(
           (me/Ephd)*(me/Ephd)*v*v/((1.0-v)*(1.0-v))
                -2.0*(me/Ephd)  *v/ (1.0-v)
                                   +(1.0-v) + 1.0/(1.0-v)
          );
  

 
  if(debug==1)
  cout<<"IC_sigma_KN: Ee="<<Ee<<" Eph="<<Eph<<" eta="<<eta<<" Egamma="<<Egamma
      <<" g="<<g<<" beta="<<beta
      <<" Ephd="<< Ephd <<" me/Ephd="<< me/Ephd<<" Ephd/me="<< Ephd/me
      <<" v="<<v<<" Egamma_max="<<Egamma_max
      << " result="<<result<<endl;


  return result;
};

/////////////////////////////////////////////////////////

double InverseCompton::IC_sigma_KN_isotropic(double Ee, double Eph,  double Egamma)
{
  // Inverse Compton differential cross-section using Klein-Nishina cross-section for isotropic target photon distribution

  // units: cm^2 MeV^-1
  // using equation  (12)     of Moskalenko & Strong ApJ 528, 357 (2000)
  // stand-alone form, all constants defined here
  // energies in MeV
  // Ee     = electron energy
  // Eph    = target photon energy
  // eta    = angle between electron and target photon in radians
  // Egamma = gamma ray energy

#include"Constants.h"

  double result;
  if(Egamma>=Ee) {result=0.; return result;}

  if(initialized!=1) init(0); // just in case


  double g    = Ee/me;                      // electron Lorentz factor
  
  double Egamma_max = 4.0*g*g*Eph/(1.0 + 4.0*Eph/me); // MS2000 eq (10)
  if(Egamma>=Egamma_max) {result=0.; return result;}

  double q_prime= Egamma/me / (4.0* Eph/me *g*g*(1.0- Egamma/me /g)); // Eph and Egamma units of me in equation (12) of reference

  if(q_prime>1.0 || q_prime<=1.0/(4.0*g*g)) {result=0.; return result;} // condition below equation (12) of reference
  

  result=2.0*pi*re*re/ (Eph*g*g)
       *(2.0*q_prime*log(q_prime) + (1.0 + 2.0*q_prime)*(1.0 - q_prime)
	 +0.5*pow(4. *Eph/me *g *q_prime,2) * (1.0-q_prime) / (1.0 + 4.0*Eph/me *g *q_prime)
	 );// NB eq (12) has Eph in units of me on RHS

  
 
  if(debug==1)
  cout<<"IC_sigma_KN_isotropic: Ee="<<Ee<<" Eph="<<Eph<<" Egamma="<<Egamma
      <<" g="<<g      <<" Egamma_max="<<Egamma_max<<" q_prime="<<q_prime
      << " result="<<result<<endl;
  

  return result;
};

///////////////////////////////////////////////////////////////////////////

  void  InverseCompton:: IC_emiss(valarray<double> Ee,   valarray<double> Eph, double eta,valarray<double> Egamma,
                                  valarray<double> fluxe,valarray<double> nph,            valarray<double> &emiss,
                                  int method )

{

  if(debug==1) cout<<">> InverseCompton::IC_emiss()"<< " method="<<method<<endl;

  // assuming both electron and target photon energy grids are logarithmic

  ///  emiss.resize(Egamma.size());

  emiss=0.0;

 if(method==1)
 {

  double IC_sigmaKN;

  for (int iEe    =0;iEe    <Ee    .size();iEe    ++)
  for (int iEph   =0;iEph   <Eph   .size();iEph   ++)
  for (int iEgamma=0;iEgamma<Egamma.size();iEgamma++)
  {
    IC_sigmaKN      = IC_sigma_KN   (Ee[iEe],  Eph[iEph], eta, Egamma[iEgamma]);

    emiss[iEgamma] += IC_sigmaKN* fluxe[iEe] * nph[iEph]
                                    *Ee[iEe] * Eph[iEph];

    if(debug==1) cout<<"IC_emiss: method="<<method<<" Egamma="<< Egamma[iEgamma]<<" Ee="<<Ee[iEe]<<" Eph ="<<Eph[iEph]
                <<"  IC_sigmaKN= "<<IC_sigmaKN<<" emiss="<<emiss[iEgamma]<<endl;
  }


  emiss *= log(Ee[1]/Ee[0])*log(Eph[1]/Eph[0]);// NB assumes at least 2 values in spectra  !

 }// method==1

 ////////////// optimized inline method

 if(method==2)
 {
#include"Constants.h"
 if(initialized!=1) init(0); // just in case

  double IC_sigmaKN;

  for (int iEe    =0;iEe    <Ee    .size();iEe    ++)
  {
  double g    = Ee[iEe]/me;                      // electron Lorentz factor
  double beta = sqrt(1.0-1.0/(g*g));        // electron velocity/c

  for (int iEph   =0;iEph   <Eph   .size();iEph   ++)
  {
 
   double Ephd       = g*Eph[iEph]*(1.0+beta*cos(eta));  // energy of target photon in electron frame
   double Egamma_max = 2.0*g*Ephd/(1.0+2.0*Ephd/me);

   double term1=  constant/(Eph[iEph]*Ee[iEe]   *Ee[iEe])
                                 * fluxe[iEe] * nph[iEph]
                                     *Ee[iEe] * Eph[iEph];

  for (int iEgamma=0;iEgamma<Egamma.size();iEgamma++)
  {
   

  if(Egamma[iEgamma] < Ee[iEe]) 
  {
    
   if(Egamma[iEgamma]< Egamma_max) 
   {

  
    double v    = Egamma[iEgamma]/Ee[iEe];                  // fractional energy transfer
    double v1   = 1.0-v;
    double v2   = v/v1;
    double v3   = v1 + 1.0/v1;

    double me_Ephd_v2= me/Ephd*v2;

 

    double term2= me_Ephd_v2*(me_Ephd_v2 -2.0) + v3;


    emiss[iEgamma] +=   term1 * term2;

     /* for reference to follow the formula above

       (me/Ephd)*(me/Ephd)*v*v/((1.0-v)*(1.0-v))
            -2.0*(me/Ephd)  *v/ (1.0-v)
                               +(1.0-v) + 1.0/(1.0-v);
    
    */

  if(debug==1)
  cout<<"IC_emiss : method="<<method
      <<" Ee="<<Ee[iEe]<<" Eph="<<Eph[iEph]<<" eta="<<eta<<" Egamma="<<Egamma[iEgamma]
      <<" g="<<g<<" beta="<<beta
      <<" Ephd="<< Ephd <<" me/Ephd="<< me/Ephd<<" v="<<v<<" Egamma_max="<<Egamma_max
      << " IC_sigmaKN="<<IC_sigmaKN<<endl;

 
     }// Egamma < Egamma_max
  }//   Egamma < Ee

    if(debug==1) cout<<"IC_emiss : method="<<method<<" Egamma="<< Egamma[iEgamma]<<" Ee="<<Ee[iEe]<<" Eph ="<<Eph[iEph]
                <<"  IC_sigmaKN= "<<IC_sigmaKN<<" emiss="<<emiss[iEgamma]<<endl;

 

    }//Egamma
   }//Eph
  }//Ee

  
  emiss *= log(Ee[1]/Ee[0])*log(Eph[1]/Eph[0]);// NB assumes at least 2 values in spectra  !

 }// method==2


  if(debug==1)  cout<<"<< InverseCompton::IC_emiss()"<<endl;
  return;
};


///////////////////////////////////////////////////////////////////////

 void InverseCompton::  IC_loss (          double  Ee,  
                  valarray<double> Eph,  valarray<double> nph,    
		  double &dEdt, double &mean_energy_loss,
                  int method)

{
#include"Constants.h"

  if(debug==1) cout<<">> InverseCompton::IC_loss()"<<endl;
  

  //  int method=1; // 1=integrate over angle, 2=use angle-integrated isotropic cross-section
       

  double factor_Egamma= 1.01; 
  double deta          =0.01;// radian

  dEdt=0.;
  mean_energy_loss=0.;
  double energy_loss_weight=0.;

  for(int iEph   =0;iEph<Eph.size();   iEph++)  
  {
 
    double     E_gamma_sum=0;//AWS20120208
    double IC_sigma_KN_sum=0;//AWS20120208

  

  for(double Egamma=1e-6; Egamma<=1e8; Egamma *=factor_Egamma) // low energy required for correct result for low Ee
  {
   if(Egamma<Ee)
   {

  if(method==1) //AWS20120303
   {

  for(double eta   =0;    eta   <= pi;     eta+=deta  )
  {
   double IC_sigmaKN=IC_sigma_KN( Ee,  Eph[iEph], eta,  Egamma);

   if(debug==1)
   if(IC_sigmaKN>0.)
    cout<<"InverseCompton::IC_loss : Ee="<<Ee<< " MeV  Eph="<<Eph[iEph]*1e6<<" eV  Egamma="<<Egamma<<" MeV  eta="<<eta
       << " rad  IC_sigma_KN="<<IC_sigmaKN<<endl;

                                                                 
   double solid_angle_factor=0.5*sin(eta)*deta; //normalized to 1 //AWS20120313

   // Egamma is log scale so need to include factor in integration
   double weight = Egamma*IC_sigmaKN*solid_angle_factor;

   E_gamma_sum     += Egamma*weight;
   IC_sigma_KN_sum +=        weight;

  }//eta

   }// method==1
 
  if(method==2) //AWS20120303
  {
  double IC_sigmaKN=IC_sigma_KN_isotropic( Ee,  Eph[iEph],   Egamma);

   double weight = Egamma*IC_sigmaKN;

   E_gamma_sum     += Egamma*weight;
   IC_sigma_KN_sum +=        weight;


  if(debug==1)
   if(IC_sigmaKN>0.)
    cout<<"InverseCompton::IC_loss : isotropic case Ee="<<Ee<< " MeV  Eph="
        <<Eph[iEph]*1e6<<" eV  Egamma="<<Egamma<<" MeV  IC_sigma_KN="<<IC_sigmaKN<<endl;
  }


   }// Egamma<Ee

  }//Egamma

  
   double IC_sigma_total=IC_sigma_KN_sum         *log(factor_Egamma);// since logarithmic integration  //AWS20120313
   double IC_sigma_Thompson=8./3.*pi*re*re;

   double E_gamma_mean=E_gamma_sum/IC_sigma_KN_sum;
   double mean_energy_transfer=E_gamma_mean/Ee;
   double IC_sigma_total_v=IC_sigma_total*mean_energy_transfer;
   double dvdt=IC_sigma_total_v*c;
   double t_loss_for_1eVcm3=1./dvdt*Eph[iEph]*1e6 /(86400.*365); // seconds->years

   double analytical_mean_energy_transfer_low = 4./3. * pow(Ee/me,2)*Eph[iEph] / Ee; // Thompson limit for comparison, only valid for Ephd<<me

   double Ephd=4./3.*(Ee/me)*Eph[iEph];// average photon energy in electron system
   double analytical_mean_energy_transfer_high = (log(2.*Ephd/me)-5./6.)/(log(2.0*Ephd/me)+1./2.);// ultra-high-energy limit, Allcock and Wdowczyk 1972, http://adsabs.harvard.edu/abs/1972NCimB...9..315A  only valid for Ephd>>me


   dEdt += IC_sigma_total * E_gamma_mean * nph[iEph]*Eph[iEph] *c; // Eph for logarithmic integration

   mean_energy_loss       +=  IC_sigma_total * E_gamma_mean/Ee * nph[iEph]*Eph[iEph];
        energy_loss_weight+=  IC_sigma_total                   * nph[iEph]*Eph[iEph];


  if(debug==1)
  {
    cout<<"InverseCompton::IC_loss: method="<<method<<" Ee="<<Ee<< " MeV Eph="<<Eph[iEph]*1e6<<" eV"<<" nph="<<nph[iEph]
       <<" total sigma="   << IC_sigma_total
       <<" Thompson sigma="<< IC_sigma_Thompson
       <<" E_gamma_mean="<<E_gamma_mean
       <<" dEdt="<<dEdt
       <<" mean energy loss="<<mean_energy_loss
       <<" energy loss weight="<<energy_loss_weight
       <<" mean fractional energy transfer="<<mean_energy_transfer;
       if(Ephd<=me) cout<<" analytical approx low energy="<<analytical_mean_energy_transfer_low;
       if(Ephd> me) cout<<" analytical ultra-high energy limit="     <<analytical_mean_energy_transfer_high;
       cout<<" t_loss for 1 eV cm-3="<<t_loss_for_1eVcm3<<" yr";
      
   cout<<endl;
  }// debug==1


  }//Eph

  dEdt *= log(Eph[1]/Eph[0]); // for logarithmic integration

  mean_energy_loss/=energy_loss_weight   ;

  return;
};


//////////////////////////////
// other functions foreseen for arrays of energies at constant eta, etc, for more efficient computation
//////////////////////////////






/////////////////////////////////////////////////////////
void InverseCompton::test(int options)
{
#include"Constants.h"
  cout<<">> InverseCompton::test()"<<endl;




  int debug_=1;
  init(debug_);


  int test_IC_sigma_KN=0;
  int test_IC_emiss   =0;
  int test_IC_loss    =0;

  if(options==1) test_IC_sigma_KN=1;
  if(options==2) test_IC_emiss   =1;
  if(options==3) test_IC_loss    =1;

  if( test_IC_sigma_KN==1)
  {

  cout<<   " testing IC_sigma_KN"<<endl;

  double Ee,  Eph, eta,  Egamma;
  Ee=1e3;
  Eph=1e-3*1.e-6;// NB requires value in MeV 
  Egamma=1e2;
  eta=0.5;

  
  
  double factor_Ee    =10.0;       //AWS20120208
  double factor_Eph   = 2.0;       //AWS20120208
  double factor_Egamma= 2.0;//1.05;//AWS20120208 AWS20131204
  double deta          =0.2;//0.01;// radian     AWS20131204


  for(Eph   =1e-10;Eph   <=1e-5;    Eph*=factor_Eph)    //AWS20120208
  {
  for(Ee    =1e2;  Ee    <=1e8;     Ee *=factor_Ee)     //AWS20120208
  {
    double     E_gamma_sum=0;//AWS20120208
    double IC_sigma_KN_sum=0;//AWS20120208
   
  

  for(Egamma=1e-3; Egamma<=1e8; Egamma *=factor_Egamma) //AWS20120208
  {

   double IC_sigma_KN_sum_eta=0;//AWS20120208

  for(eta   =0;    eta   <= pi;     eta+=deta  )
  {
   double IC_sigmaKN=IC_sigma_KN( Ee,  Eph, eta,  Egamma);

   if(IC_sigmaKN>0.)
    cout<<"InverseCompton::test : Ee="<<Ee<< " MeV  Eph="<<Eph*1e6<<" eV  Egamma="<<Egamma<<" MeV  eta="<<eta
       << " rad  IC_sigma_KN="<<IC_sigmaKN<<endl;

   double solid_angle_factor=2.0*pi*sin(eta)*deta;            //AWS20120209
 

   // Egamma is log scale so need to include factor in integration
   double weight = Egamma*IC_sigmaKN*solid_angle_factor;

   E_gamma_sum         += Egamma*weight;//AWS20120208
   IC_sigma_KN_sum     +=        weight;//AWS20120208

   IC_sigma_KN_sum_eta += IC_sigmaKN*solid_angle_factor;//AWS20120224

  }//eta
   cout<<endl;

   double IC_sigmaKN_isotropic=IC_sigma_KN_isotropic( Ee,  Eph,  Egamma);
   double IC_sigma_total_eta=IC_sigma_KN_sum_eta/(4.0*pi);

   if(IC_sigmaKN_isotropic>0.)
    cout<<"InverseCompton::test isotropic : Ee="<<Ee<< " MeV  Eph="<<Eph*1e6<<" eV  Egamma="<<Egamma<<" MeV "
       << " IC_sigma from angle integration = "<<IC_sigma_total_eta
       << " IC_sigma analytical = "            <<IC_sigmaKN_isotropic
       << " ratio ="<<IC_sigma_total_eta/ IC_sigmaKN_isotropic	
       <<endl;

  }//Egamma

   double IC_sigma_total=IC_sigma_KN_sum/(4.0*pi)*log(factor_Egamma);// since logarithmic integration
   double IC_sigma_Thompson=8./3.*pi*re*re;

   double E_gamma_mean=E_gamma_sum/IC_sigma_KN_sum;
   double mean_energy_transfer=E_gamma_mean/Ee;
   double IC_sigma_total_v=IC_sigma_total*mean_energy_transfer;
   double dvdt=IC_sigma_total_v*c;
   double t_loss_for_1eVcm3=1./dvdt*Eph*1e6 /(86400.*365); // seconds->years

   double analytical_mean_energy_transfer_low = 4./3. * pow(Ee/me,2)*Eph / Ee; // Thompson limit for comparison, only valid for Ephd<<me

   double Ephd=4./3.*(Ee/me)*Eph;// average photon energy in electron system
   double analytical_mean_energy_transfer_high = (log(2.*Ephd/me)-5./6.)/(log(2.0*Ephd/me)+1./2.);// ultra-high-energy limit, Allcock and Wdowczyk 1972, http://adsabs.harvard.edu/abs/1972NCimB...9..315A  only valid for Ephd>>me

   cout<<"InverseCompton::test IC_sigma_KN:  Ee="<<Ee<< " MeV Eph="<<Eph*1e6<<" eV"
       <<" total sigma="   << IC_sigma_total
       <<" Thompson sigma="<< IC_sigma_Thompson
       <<" E_gamma_mean="<<E_gamma_mean
       <<" mean fractional energy transfer="<<mean_energy_transfer;
       if(Ephd<=me) cout<<" analytical approx low energy="<<analytical_mean_energy_transfer_low;
       if(Ephd> me) cout<<" analytical ultra-high energy limit="     <<analytical_mean_energy_transfer_high;
       cout<<" t_loss for 1 eV cm-3="<<t_loss_for_1eVcm3<<" yr";
      
   cout<<endl;
  }//Ee
   cout<<endl;
  }//Eph

  }//if  test_IC_sigma_KN==1


//----------------------------------- test emissivity

  if(test_IC_emiss==1)
  {

  cout<<" testing emissivity"<<endl;

  int nEe    =400;   // NB must be >1 for correct run
  int nEph   =200   ; // NB must be >1 for correct run
  int nEgamma=100 ;

  double  Eemin       =1.0;
  double  Eefactor    =1.05;
  double  Ephmin      =1.0e-4 * 1e-6;//NB MeV
  double  Ephfactor   =1.05;
  double  Egammamin   =1.0;
  double  Egammafactor=pow(10.,0.1);

  double ge=3.0; // electron spectral index
 
  valarray<double> Ee_grid     (nEe);
  valarray<double> Eph_grid    (nEph);
  valarray<double> Egamma_grid (nEgamma);
  valarray<double> fluxe       (nEe);
  valarray<double> nph         (nEph);
  valarray<double> emiss       (nEgamma);

  int method = 2;
  debug=0;
  

  for (int iEgamma=0;iEgamma<Egamma_grid.size();iEgamma++) Egamma_grid [iEgamma]=Egammamin*pow(Egammafactor,iEgamma);
 
  for (int iEph   =0;iEph   <Eph_grid   .size();iEph   ++) Eph_grid    [iEph]   =Ephmin   *pow(Ephfactor,   iEph   );
  for (int iEph   =0;iEph   <Eph_grid   .size();iEph   ++) nph         [iEph]   = 1.0e6;// per MeV

  for (int iEe    =0;iEe    <Ee_grid    .size();iEe    ++) Ee_grid     [iEe]    =Eemin    *pow(Eefactor,    iEe    );
  for (int iEe    =0;iEe    <Ee_grid    .size();iEe    ++) fluxe       [iEe]    =          pow(Ee_grid[iEe], -ge   );

  // Ackermann etal Fermi-LAT arXiv:1109.052 electrons fit, converting from GeV m^2 to MeV cm^2 
  ge=3.19;
  for (int iEe    =0;iEe    <Ee_grid    .size();iEe    ++) fluxe       [iEe]    = 2.07e-9* pow(Ee_grid[iEe]/2.e4, -ge   );
   for (int iEe    =0;iEe    <Ee_grid    .size();iEe++)     cout<<"Ee= "    << Ee_grid    [iEe]<<" fluxe="<<fluxe[iEe]<<endl;
   for (int iEph   =0;iEph   <Eph_grid   .size();iEph   ++) cout<<"Eph= "   << Eph_grid   [iEph]<<" nph=" <<  nph[iEph]<<endl;

  double eta=0.1; // 0= head on

  for(method=1;method<=2;method++)
  {
   IC_emiss(Ee_grid, Eph_grid, eta, Egamma_grid,
	   fluxe,   nph,           emiss,
           method );

   

   for (int iEgamma=0;iEgamma<Egamma_grid.size();iEgamma++)
       cout<<" method="<<method<<" Egamma= "<< Egamma_grid[iEgamma]<<" emiss="<<emiss[iEgamma]<<endl;

  }//method
  

  }// if test_IC_emiss==1





  //-------------------------------- test IC_loss
  if(test_IC_loss==1)
  {
  cout<<"testing IC_loss"<<endl;

  debug=0;
 
  int nEe    =20;   // NB must be >1 for correct run
  int nEph   =200   ; // NB must be >1 for correct run
 

  double  Eemin       =1.0e0;
  double  Eefactor    =pow(10,0.5);
  double  Ephmin      =1.0e-4 * 1e-6;//NB MeV
  double  Ephfactor   =1.05;

  valarray<double> Ee_grid     (nEe);
  valarray<double> Eph_grid    (nEph);
  valarray<double> nph         (nEph);

  for (int iEph   =0;iEph   <Eph_grid   .size();iEph   ++) Eph_grid    [iEph]   =Ephmin   *pow(Ephfactor,   iEph   );
  for (int iEph   =0;iEph   <Eph_grid   .size();iEph   ++) nph         [iEph]   = 1.0e6;// per MeV
  for (int iEe    =0;iEe    <Ee_grid    .size();iEe    ++) Ee_grid     [iEe]    =Eemin    *pow(Eefactor,    iEe    );

  for (int iEph   =0;iEph   <Eph_grid   .size();iEph   ++) 
     cout<<"Radiation field: Eph="<<Eph_grid[iEph]*1.0e6<<" eV, nph="<<nph[iEph]*1e-6<<" cm-3 eV-1"<<endl;

  double total_energy_density=0.;
  for (int iEph   =0;iEph   <Eph_grid   .size();iEph   ++)  total_energy_density+=nph[iEph]*Eph_grid[iEph]*Eph_grid[iEph];
  total_energy_density*=log(Eph_grid[1]/Eph_grid[0])*1.0e6; // since log scale, MeV -> eV
  cout<<"Radiation field: total energy density = "<<total_energy_density<< " eV cm-3"<<endl;
  

  // normalize to 1 eV cm-3
  nph /= total_energy_density;

  for (int iEph   =0;iEph   <Eph_grid   .size();iEph   ++) 
     cout<<"Radiation field after normalization: Eph="<<Eph_grid[iEph]*1.0e6<<" eV, nph="<<nph[iEph]*1e-6<<" cm-3 eV-1"<<endl;

   total_energy_density=0.;
  for (int iEph   =0;iEph   <Eph_grid   .size();iEph   ++)  total_energy_density+=nph[iEph]*Eph_grid[iEph]*Eph_grid[iEph];
  total_energy_density*=log(Eph_grid[1]/Eph_grid[0])*1.0e6; // since log scale, MeV -> eV
  cout<<"Radiation field: total energy density after normalization= "<<total_energy_density<< " eV cm-3"<<endl;

  double dEdt,mean_energy_loss;

  for (int iEe    =0;iEe    <Ee_grid    .size();iEe    ++)
  { 
    int method=1;
    IC_loss(Ee_grid[iEe], Eph_grid, nph, dEdt, mean_energy_loss, method);

    double loss_time=Ee_grid[iEe]/dEdt/(86400.*365); // in years

    cout<<"InverseCompton::test IC_loss: method="<<method<<" Ee="
        <<    Ee_grid[iEe]<<" MeV  dEdt="<<dEdt<< " MeV s-1  mean fractional energy loss="
        <<mean_energy_loss<<" loss time="<<loss_time<< " yr"<<endl;

       method=2;
    IC_loss(Ee_grid[iEe], Eph_grid, nph, dEdt, mean_energy_loss, method);

           loss_time=Ee_grid[iEe]/dEdt/(86400.*365); // in years

    cout<<"InverseCompton::test IC_loss: method="<<method<<" Ee="
        <<    Ee_grid[iEe]<<" MeV  dEdt="<<dEdt<< " MeV s-1  mean fractional energy loss="
        <<mean_energy_loss<<" loss time="<<loss_time<< " yr"<<endl;

  }

  }// test_IC_loss==1

  cout<<"<< InverseCompton::test()"<<endl;
  return;
}

///////////////////////////// Test program: to run,  remove /* ... */ and then do:  g++ InverseCompton.cc ; ./a.out
/*
int main()
{
cout<<">>InverseCompton::test"<<endl;

  InverseCompton inversecompton;

  int options;                     //AWS20131204

  options=1;                       //AWS20131204
  inversecompton.test(options);    //AWS20131204

  options=2;                       //AWS20131204
  inversecompton.test(options);    //AWS20131204

  options=3;                       //AWS20131204
  inversecompton.test(options);    //AWS20131204


  cout<<"<<InverseCompton::test"<<endl;
  return 0;
}

*/
